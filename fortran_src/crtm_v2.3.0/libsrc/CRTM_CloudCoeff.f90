!
! CRTM_CloudCoeff
!
! Module containing the shared CRTM scattering coefficient data
! (CloudCoeff) and their load/destruction routines. 
!
! PUBLIC DATA:
!       CloudC:  Data structure containing the cloud bulk optical
!                properties data
! ----- Added by Yinghui Lu ------------
!       CloudC_ice, CloudC_water, CloudC_rain, CloudC_snow, CloudC_graup,
!       CloudC_hail :  
!                Data structure containing the cloud bulk optical
!                properties data, for ice/water/rain/snow/graupel/hail
!                respectively, for MP-scheme-specific LUTs
! ----- END Added by Yinghui Lu ------------
!
! SIDE EFFECTS:
!       Routines in this module modify the contents of the public
!       data structure CloudC.
!
! RESTRICTIONS:
!       Routines in this module should only be called during the
!       CRTM initialisation.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 24-Jun-2004
!                       paul.vandelst@noaa.gov
!

MODULE CRTM_CloudCoeff

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use
  USE Message_Handler,      ONLY: SUCCESS, FAILURE, Display_Message
  USE CloudCoeff_Define,    ONLY: CloudCoeff_type, &
                                  CloudCoeff_Associated, &
                                  CloudCoeff_Destroy                  
  USE CloudCoeff_Binary_IO, ONLY: CloudCoeff_Binary_ReadFile
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! The shared data
  PUBLIC :: CloudC
  ! ----- Added by Yinghui Lu ------------
  PUBLIC :: CloudC_ice
  PUBLIC :: CloudC_water
  PUBLIC :: CloudC_rain
  PUBLIC :: CloudC_snow
  PUBLIC :: CloudC_graup
  PUBLIC :: CloudC_hail
  ! ----- END Added by Yinghui Lu ------------
  ! Procedures
  PUBLIC :: CRTM_CloudCoeff_Load
  PUBLIC :: CRTM_CloudCoeff_Destroy
  PUBLIC :: CRTM_CloudCoeff_IsLoaded


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: CRTM_CloudCoeff.f90 99117 2017-11-27 18:37:14Z tong.zhu@noaa.gov $'
  ! Message string length
  INTEGER, PARAMETER :: ML = 256


  ! ---------------------------------
  ! The shared cloud coefficient data
  ! ---------------------------------
  TYPE(CloudCoeff_type), TARGET, SAVE :: CloudC
  ! ----- Added by Yinghui Lu ------------
  TYPE(CloudCoeff_type), TARGET, SAVE :: CloudC_ice
  TYPE(CloudCoeff_type), TARGET, SAVE :: CloudC_water
  TYPE(CloudCoeff_type), TARGET, SAVE :: CloudC_rain
  TYPE(CloudCoeff_type), TARGET, SAVE :: CloudC_snow
  TYPE(CloudCoeff_type), TARGET, SAVE :: CloudC_graup
  TYPE(CloudCoeff_type), TARGET, SAVE :: CloudC_hail
  ! ----- END Added by Yinghui Lu ------------


CONTAINS


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_CloudCoeff_Load
!
! PURPOSE:
!       Function to load the CloudCoeff scattering coefficient data into
!       the public data structure CloudC.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_CloudCoeff_Load( &
!                        Filename,                              &
!                        File_Path         = File_Path        , &
!                        Quiet             = Quiet            , &
!                        Process_ID        = Process_ID       , &
!                        Output_Process_ID = Output_Process_ID  )
!
! INPUT ARGUMENTS:
!       Filename:           Name of the Binary format CloudCoeff file.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!
! OPTIONAL INPUT ARGUMENTS:
!       File_Path:          Character string specifying a file path for the
!                           input data file. If not specified, the current
!                           directory is the default.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Quiet:              Set this logical argument to suppress INFORMATION
!                           messages being printed to stdout
!                           If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                              == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                           If not specified, default is .FALSE.
!                           UNITS:      N/A
!                           TYPE:       LOGICAL
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Process_ID:         Set this argument to the MPI process ID that this
!                           function call is running under. This value is used
!                           solely for controlling INFORMATIOn message output.
!                           If MPI is not being used, ignore this argument.
!                           This argument is ignored if the Quiet argument is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Output_Process_ID:  Set this argument to the MPI process ID in which
!                           all INFORMATION messages are to be output. If
!                           the passed Process_ID value agrees with this value
!                           the INFORMATION messages are output. 
!                           This argument is ignored if the Quiet argument
!                           is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
! ----- Added by Yinghui Lu ------------
!       Filename_ice/water/rain/snow/graupel/hail: 
!                           Name of the Binary format CloudCoeff file.
!                           MP-scheme-specific
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
! ----- END Added by Yinghui Lu ------------
!
! FUNCTION RESULT:
!       Error_Status:       The return value is an integer defining the error
!                           status. The error codes are defined in the
!                           Message_Handler module.
!                           If == SUCCESS the CloudCoeff data load was successful
!                              == FAILURE an unrecoverable error occurred.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data structure CloudC.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_CloudCoeff_Load( &
    Filename         , &  ! Input
    File_Path        , &  ! Optional input
    Quiet            , &  ! Optional input
    Process_ID       , &  ! Optional input
    Output_Process_ID, &  ! Optional input
    ! ----- Added by Yinghui Lu ------------
    Filename_ice     , &  ! Optional Input
    Filename_water   , &  ! Optional Input
    Filename_rain    , &  ! Optional Input
    Filename_snow    , &  ! Optional Input
    Filename_graup   , &  ! Optional Input
    Filename_hail    ) &  ! Optional Input
    ! ----- END Added by Yinghui Lu ------------
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN) :: Filename
    CHARACTER(*), OPTIONAL, INTENT(IN) :: File_Path
    LOGICAL     , OPTIONAL, INTENT(IN) :: Quiet             
    INTEGER     , OPTIONAL, INTENT(IN) :: Process_ID
    INTEGER     , OPTIONAL, INTENT(IN) :: Output_Process_ID
    ! ----- Added by Yinghui Lu ------------
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Filename_ice
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Filename_water
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Filename_rain
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Filename_snow
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Filename_graup
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Filename_hail
    ! ----- End Added by Yinghui Lu ------------
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_CloudCoeff_Load'
    ! Local variables
    CHARACTER(ML) :: msg, pid_msg
    CHARACTER(ML) :: CloudCoeff_File
    LOGICAL :: noisy

    ! ----- Added by Yinghui Lu ------------
    CHARACTER(ML) :: CloudCoeff_File_ice
    CHARACTER(ML) :: CloudCoeff_File_water
    CHARACTER(ML) :: CloudCoeff_File_rain
    CHARACTER(ML) :: CloudCoeff_File_snow
    CHARACTER(ML) :: CloudCoeff_File_graup
    CHARACTER(ML) :: CloudCoeff_File_hail
    ! Should distribution specific LUT be used instead of default one?
    LOGICAL :: Use_LUT_ice   
    LOGICAL :: Use_LUT_water
    LOGICAL :: Use_LUT_rain
    LOGICAL :: Use_LUT_snow
    LOGICAL :: Use_LUT_graup
    LOGICAL :: Use_LUT_hail

    CloudCoeff_File_ice   = ''
    CloudCoeff_File_water = ''
    CloudCoeff_File_rain  = ''
    CloudCoeff_File_snow  = ''
    CloudCoeff_File_graup = ''
    CloudCoeff_File_hail  = ''
    ! ----- END Added by Yinghui Lu ------------
    
    ! Setup 
    err_stat = SUCCESS
    ! ...Assign the filename to local variable
    CloudCoeff_File = ADJUSTL(Filename)
    ! ...Add the file path
    IF ( PRESENT(File_Path) ) CloudCoeff_File = TRIM(ADJUSTL(File_Path))//TRIM(CloudCoeff_File)

    ! ----- Added by Yinghui Lu ------------
    ! Let's check if distribution-specific LUT should be used.
    ! Use distribution-specific LUT only if optional input argument
    ! filename_xxx is present and not equal to empty string.

    USE_LUT_ice   = .FALSE.
    IF ( PRESENT(filename_ice  ) ) THEN
      IF (.NOT. (filename_ice   == '')) USE_LUT_ice   = .TRUE.
    END IF

    USE_LUT_water = .FALSE.
    IF ( PRESENT(filename_water) ) THEN
      IF (.NOT. (filename_water == '')) USE_LUT_water = .TRUE.
    END IF

    USE_LUT_rain  = .FALSE.
    IF ( PRESENT(filename_rain ) ) THEN
      IF (.NOT. (filename_rain  == '')) USE_LUT_rain  = .TRUE.
    END IF

    USE_LUT_snow  = .FALSE.
    IF ( PRESENT(filename_snow ) ) THEN
      IF (.NOT. (filename_snow  == '')) USE_LUT_snow  = .TRUE.
    END IF

    USE_LUT_graup = .FALSE.
    IF ( PRESENT(filename_graup) ) THEN
      IF (.NOT. (filename_graup == '')) USE_LUT_graup = .TRUE.
    END IF

    USE_LUT_hail  = .FALSE.
    IF ( PRESENT(filename_hail ) ) THEN
      IF (.NOT. (filename_hail  == '')) USE_LUT_hail  = .TRUE.
    END IF

    IF ( Use_LUT_ice    ) CloudCoeff_File_ice    = TRIM(ADJUSTL(Filename_ice))
    IF ( Use_LUT_water  ) CloudCoeff_File_water  = TRIM(ADJUSTL(Filename_water))
    IF ( Use_LUT_rain   ) CloudCoeff_File_rain   = TRIM(ADJUSTL(Filename_rain))
    IF ( Use_LUT_snow   ) CloudCoeff_File_snow   = TRIM(ADJUSTL(Filename_snow))
    IF ( Use_LUT_graup  ) CloudCoeff_File_graup  = TRIM(ADJUSTL(Filename_graup))
    IF ( Use_LUT_hail   ) CloudCoeff_File_hail   = TRIM(ADJUSTL(Filename_hail))
    
    IF ( PRESENT(File_Path) ) then
      IF ( Use_LUT_ice    ) CloudCoeff_File_ice   = TRIM(ADJUSTL(File_Path))//TRIM(CloudCoeff_File_ice)
      IF ( Use_LUT_water  ) CloudCoeff_File_water = TRIM(ADJUSTL(File_Path))//TRIM(CloudCoeff_File_water)
      IF ( Use_LUT_rain   ) CloudCoeff_File_rain  = TRIM(ADJUSTL(File_Path))//TRIM(CloudCoeff_File_rain)
      IF ( Use_LUT_snow   ) CloudCoeff_File_snow  = TRIM(ADJUSTL(File_Path))//TRIM(CloudCoeff_File_snow)
      IF ( Use_LUT_graup  ) CloudCoeff_File_graup = TRIM(ADJUSTL(File_Path))//TRIM(CloudCoeff_File_graup)
      IF ( Use_LUT_hail   ) CloudCoeff_File_hail  = TRIM(ADJUSTL(File_Path))//TRIM(CloudCoeff_File_hail)
    END IF
    ! ----- End Added by Yinghui Lu ------------
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Check the MPI Process Ids
    IF ( noisy .AND. PRESENT(Process_ID) .AND. PRESENT(Output_Process_ID) ) THEN
      IF ( Process_Id /= Output_Process_Id ) noisy = .FALSE.
    END IF
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF
    
    ! Read the CloudCoeff data file
    err_stat = CloudCoeff_Binary_ReadFile( &
                 CloudCoeff_File, &
                 CloudC, &
                 Quiet = .NOT. noisy )
    IF ( err_stat /= SUCCESS ) THEN
      WRITE( msg,'("Error reading CloudCoeff file ",a)') TRIM(CloudCoeff_File)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
      RETURN
    END IF
    ! ----- Added by Yinghui Lu ------------
    ! Read the MP-scheme-specific CloudCoeff data file, if corresponding input
    ! filename is present
    IF (USE_LUT_ice) THEN
      err_stat = CloudCoeff_Binary_ReadFile( &
                   CloudCoeff_File_ice, &
                   CloudC_ice, &
                   Quiet = .NOT. noisy )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error reading CloudCoeff file ",a)') &
                TRIM(CloudCoeff_File_ice)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
    END IF

    IF (USE_LUT_water)THEN
      err_stat = CloudCoeff_Binary_ReadFile( &
                   CloudCoeff_File_water, &
                   CloudC_water, &
                   Quiet = .NOT. noisy )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error reading CloudCoeff file ",a)') &
                TRIM(CloudCoeff_File_water)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
    END IF

    IF (USE_LUT_rain ) THEN
      err_stat = CloudCoeff_Binary_ReadFile( &
                   CloudCoeff_File_rain, &
                   CloudC_rain, &
                   Quiet = .NOT. noisy )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error reading CloudCoeff file ",a)') &
                TRIM(CloudCoeff_File_rain)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
    END IF

    IF (USE_LUT_snow ) THEN
      err_stat = CloudCoeff_Binary_ReadFile( &
                   CloudCoeff_File_snow, &
                   CloudC_snow, &
                   Quiet = .NOT. noisy )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error reading CloudCoeff file ",a)') &
                TRIM(CloudCoeff_File_snow)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
    END IF

    IF (USE_LUT_graup ) THEN
      err_stat = CloudCoeff_Binary_ReadFile( &
                   CloudCoeff_File_graup, &
                   CloudC_graup, &
                   Quiet = .NOT. noisy )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error reading CloudCoeff file ",a)') &
                TRIM(CloudCoeff_File_graup)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
    END IF

    IF (USE_LUT_hail ) THEN
      err_stat = CloudCoeff_Binary_ReadFile( &
                   CloudCoeff_File_hail, &
                   CloudC_hail, &
                   Quiet = .NOT. noisy )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error reading CloudCoeff file ",a)') &
                TRIM(CloudCoeff_File_hail)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
    END IF
    ! ----- End Added by Yinghui Lu ------------

  END FUNCTION CRTM_CloudCoeff_Load


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_CloudCoeff_Destroy
!
! PURPOSE:
!       Function to deallocate the public data structure CloudC containing
!       the CRTM CloudCoeff scattering coefficient data.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_CloudCoeff_Destroy( Process_ID =Process_ID )
!
! OPTIONAL INPUT ARGUMENTS:
!       Process_ID:       Set this argument to the MPI process ID that this
!                         function call is running under. This value is used
!                         solely for controlling message output. If MPI is not
!                         being used, ignore this argument.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:     The return value is an integer defining the error
!                         status. The error codes are defined in the
!                         Message_Handler module.
!                         If == SUCCESS the deallocation of the public CloudC data
!                                       structure was successful
!                            == FAILURE an unrecoverable error occurred.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data structure CloudC.
!
!------------------------------------------------------------------------------

  FUNCTION CRTM_CloudCoeff_Destroy( Process_ID ) RESULT( err_stat )
    ! Arguments
    INTEGER, OPTIONAL, INTENT(IN) :: Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_CloudCoeff_Destroy'
    ! Local variables
    CHARACTER(ML) :: msg, pid_msg

    ! Setup
    err_stat = SUCCESS
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF

    ! Destroy the structure
    CALL CloudCoeff_Destroy( CloudC )
    IF ( CloudCoeff_Associated( CloudC ) ) THEN
      err_stat = FAILURE
      msg = 'Error deallocating CloudCoeff shared data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,msg,err_stat )
      RETURN
    END IF
   
    ! ----- Added by Yinghui Lu ------------
    ! Destroy the MP-scheme-specific structures

    IF ( CloudCoeff_Associated( CloudC_ice ) ) THEN
      CALL CloudCoeff_Destroy( CloudC_ice )
      IF ( CloudCoeff_Associated( CloudC_ice ) ) THEN
        err_stat = FAILURE
        msg = 'Error deallocating CloudCoeff shared data structure'//TRIM(pid_msg)
        CALL Display_Message( ROUTINE_NAME,msg,err_stat )
        RETURN
      END IF
    END IF

    IF ( CloudCoeff_Associated( CloudC_water ) ) THEN
      CALL CloudCoeff_Destroy( CloudC_water )
      IF ( CloudCoeff_Associated( CloudC_water ) ) THEN
        err_stat = FAILURE
        msg = 'Error deallocating CloudCoeff shared data structure'//TRIM(pid_msg)
        CALL Display_Message( ROUTINE_NAME,msg,err_stat )
        RETURN
      END IF
    END IF

    IF ( CloudCoeff_Associated( CloudC_rain ) ) THEN
      CALL CloudCoeff_Destroy( CloudC_rain )
      IF ( CloudCoeff_Associated( CloudC_rain ) ) THEN
        err_stat = FAILURE
        msg = 'Error deallocating CloudCoeff shared data structure'//TRIM(pid_msg)
        CALL Display_Message( ROUTINE_NAME,msg,err_stat )
        RETURN
      END IF
    END IF

    IF ( CloudCoeff_Associated( CloudC_snow ) ) THEN
      CALL CloudCoeff_Destroy( CloudC_snow )
      IF ( CloudCoeff_Associated( CloudC_snow ) ) THEN
        err_stat = FAILURE
        msg = 'Error deallocating CloudCoeff shared data structure'//TRIM(pid_msg)
        CALL Display_Message( ROUTINE_NAME,msg,err_stat )
        RETURN
      END IF
    END IF

    IF ( CloudCoeff_Associated( CloudC_graup ) ) THEN
      CALL CloudCoeff_Destroy( CloudC_graup )
      IF ( CloudCoeff_Associated( CloudC_graup ) ) THEN
        err_stat = FAILURE
        msg = 'Error deallocating CloudCoeff shared data structure'//TRIM(pid_msg)
        CALL Display_Message( ROUTINE_NAME,msg,err_stat )
        RETURN
      END IF
    END IF

    IF ( CloudCoeff_Associated( CloudC_hail ) ) THEN
      CALL CloudCoeff_Destroy( CloudC_hail )
      IF ( CloudCoeff_Associated( CloudC_hail ) ) THEN
        err_stat = FAILURE
        msg = 'Error deallocating CloudCoeff shared data structure'//TRIM(pid_msg)
        CALL Display_Message( ROUTINE_NAME,msg,err_stat )
        RETURN
      END IF
    END IF
    ! ----- End Added by Yinghui Lu ------------
  END FUNCTION CRTM_CloudCoeff_Destroy


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_CloudCoeff_IsLoaded
!
! PURPOSE:
!       Function to test if the CloudCoeff scattering coefficient data has
!       loaded into the public data structure CloudC.
!
! CALLING SEQUENCE:
!       status = CRTM_CloudCoeff_IsLoaded()
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_CloudCoeff_IsLoaded() RESULT( IsLoaded )
    LOGICAL :: IsLoaded
    IsLoaded = CloudCoeff_Associated( CloudC )
  END FUNCTION CRTM_CloudCoeff_IsLoaded

END MODULE CRTM_CloudCoeff
