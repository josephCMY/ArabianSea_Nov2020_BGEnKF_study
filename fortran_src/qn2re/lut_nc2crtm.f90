module lut_nc2crtm

  use CloudCoeff_Binary_IO, only: CloudCoeff_Binary_ReadFile, &
                                  CloudCoeff_Binary_WriteFile
  use CloudCoeff_Define,    only: CloudCoeff_type, &
                                  CloudCoeff_Associated, &
                                  CloudCoeff_Destroy, &
                                  CloudCoeff_Create, &
                                  CloudCoeff_ValidRelease, &
                                  CloudCoeff_Info
  use netcdf
  implicit none
  
  
  character(len=*), parameter :: SRC_FILENAME = 'lut_nc2crtm.f90' 
  integer, parameter :: LUT_PHASE_SOLID  = 1
  integer, parameter :: LUT_PHASE_LIQUID = 2
contains
  !------------------------------------------!
  subroutine nc2binary( filename_nc, hydro_phase, filename_LUT)
    character(len=*),      intent(in)    :: filename_nc
    integer,               intent(in)    :: hydro_phase
    character(len=*),      intent(in)    :: filename_LUT
   
    ! local variables
    character(len=*), parameter :: SRC_FUNCNAME = 'nc2cloudcoeff'  
    TYPE(CloudCoeff_type) :: cloudcoeff
    integer :: status
  
    call nc2cloudcoeff( filename_nc, hydro_phase, cloudcoeff )

    status = CloudCoeff_Binary_WriteFile( filename_LUT, cloudcoeff ) 

    call CloudCoeff_Destroy( cloudcoeff )

  end subroutine nc2binary
  !------------------------------------------!
  subroutine nc2cloudcoeff(filename_nc, hydro_phase, cloudcoeff)
    character(len=*),      intent(in)    :: filename_nc
    integer,               intent(in)    :: hydro_phase
    TYPE(CloudCoeff_type), intent(inout) :: cloudcoeff

    ! local variables
    character(len=*), parameter :: SRC_FUNCNAME = 'nc2cloudcoeff'  
    
    integer :: ncid
    integer :: n_freq
    integer :: n_temp
    integer :: n_Re
    integer :: n_Legendre_crtm
    integer :: n_Legendre_nc
    
    real, allocatable :: freq_mw(:)
    real, allocatable :: temperature(:)
    real, allocatable :: Re_mw(:)
    real, allocatable :: ke_mw(:,:,:)
    real, allocatable :: w_mw (:,:,:)
    real, allocatable :: g_mw (:,:,:)
    real, allocatable :: pcoeff_mw(:,:,:,:)
   
    integer :: iTemp, iFreq, iReff

    integer, parameter :: LOFFSET_4  = 0
    integer, parameter :: LOFFSET_6  = 5
    integer, parameter :: LOFFSET_8  = 12
    integer, parameter :: LOFFSET_16 = 21

    print *, filename_nc
    if ( CloudCoeff_Associated( cloudcoeff)  ) then
      print *, "Error: Source file: ", SRC_FILENAME, " function name: ", SRC_FUNCNAME
      print *, 'Info: CloudCoeff already associated'
      stop

    end if


    call handle_nc_error( nf90_open( filename_nc, NF90_NOWRITE, ncid ), &
                          SRC_FILENAME, SRC_FUNCNAME )
    call lut_nc_inq_dim( ncid, n_freq, n_temp, n_Re, n_Legendre_nc)
    print *,ncid, n_freq, n_temp, n_Re, n_Legendre_nc

    !TODO: why?
    n_Legendre_crtm = 38
    ! use newly-generated MW LUT data 
    ! For IR data, use some dummy numbers.
    ! These LUT files are not used in IR calculation anyway
    
    ! For ice hydrometers, temperature are also needed -- they are used in the
    ! interpolation routine in CRTM
    if (hydro_phase == LUT_PHASE_SOLID) then
      call CloudCoeff_Create( CloudCoeff, &
                              n_freq, &
                              n_Re, &
                              1, &        ! n_IR_freq
                              1, &        ! n_IR_Radii
                              2, & ! n_temp, dummy
                              1, &        ! n_density
                              n_Legendre_crtm, &
                              1)          ! n_phase_elements
    else if (hydro_phase == LUT_PHASE_LIQUID) then
      call CloudCoeff_Create( CloudCoeff, &
                              n_freq, &
                              n_Re, &
                              1, &        ! n_IR_freq
                              1, &        ! n_IR_Radii
                              n_temp, &
                              1, &        ! n_density
                              n_Legendre_crtm, &
                              1)          ! n_phase_elements
    end if
    ! read dimension 
    CloudCoeff%Frequency_IR = 200. ! dummy
    CloudCoeff%Reff_IR      = 2.   ! dummy
    
    CloudCoeff%Density = 1 ! all LUTs has one density since each ice type has
                           ! its own LUT

    ! read variables
    allocate( freq_mw    ( n_freq), &
              temperature( n_temp), &
              Re_mw      ( n_Re),   &
              ke_mw      ( n_Re, n_temp,                 n_freq ), &
              w_mw       ( n_Re, n_temp,                 n_freq ), &
              g_mw       ( n_Re, n_temp,                 n_freq ), &
              pcoeff_mw  ( n_Re, n_temp, (n_Legendre_nc+1), n_freq ) )

    call lut_nc_read_var( ncid,        &
                          freq_mw,     &
                          temperature, &
                          Re_mw,       &
                          ke_mw,       &
                          w_mw,        &
                          g_mw,        &
                          pcoeff_mw )
    print *, pcoeff_Mw(2,1,:,2)
    ! map dimensions to variables in CloudCoeff 
    select case (hydro_phase)
      case (LUT_PHASE_SOLID)
        cloudcoeff%Frequency_MW = freq_mw
        cloudcoeff%Reff_MW      = Re_mw *1e6
        cloudcoeff%Temperature  = (/273.15, 283.15/) ! dummy
        cloudcoeff%Density      = 1
        do iFreq = 1, n_freq
          cloudcoeff%ke_S_MW(iFreq, :, 1) = ke_mw(:, 1, iFreq)
          cloudcoeff%w_S_MW (iFreq, :, 1) = w_mw (:, 1, iFreq)
          cloudcoeff%g_S_MW (iFreq, :, 1) = g_mw (:, 1, iFreq)
          do iReff = 1, n_Re
            cloudcoeff%pcoeff_S_MW(iFreq, iReff, 1, LOFFSET_4 :(LOFFSET_4+4), 1) = &
                               pcoeff_mw (iReff, 1,         1 :5, iFreq)
            cloudcoeff%pcoeff_S_MW(iFreq, iReff, 1, LOFFSET_6 :(LOFFSET_6+6), 1) = &
                               pcoeff_mw (iReff, 1,         1 :7, iFreq)
            cloudcoeff%pcoeff_S_MW(iFreq, iReff, 1, LOFFSET_8 :(LOFFSET_8+8), 1) = &
                               pcoeff_mw (iReff, 1,          1:9, iFreq)
            cloudcoeff%pcoeff_S_MW(iFreq, iReff, 1, LOFFSET_16:(LOFFSET_16+16), 1) = &
                               pcoeff_mw (iReff, 1,          1:17, iFreq)
          end do
        end do
!print *, cloudcoeff%ke_S_MW(1, :, 1) 
!print *, '-----------'
!print *, ke_mw(:, 1, 1)
print *, '==========='
      case (LUT_PHASE_LIQUID)
        cloudcoeff%Frequency_MW = freq_mw
        cloudcoeff%Reff_MW      = Re_mw *1e6
        cloudcoeff%Temperature  = Temperature
        cloudcoeff%Density      = 1
        do iFreq = 1, n_freq
          do iTemp = 1, n_Temp
            cloudcoeff%ke_L_MW(iFreq, :, iTemp) = ke_mw(:, iTemp, iFreq)
            cloudcoeff%w_L_MW (iFreq, :, iTemp) = w_mw (:, iTemp, iFreq)
            cloudcoeff%g_L_MW (iFreq, :, iTemp) = g_mw (:, iTemp, iFreq)
            do iReff = 1, n_Re
              cloudcoeff%pcoeff_L_MW(iFreq, iReff, iTemp, LOFFSET_4:(LOFFSET_4+4), 1) = &
                                 pcoeff_mw (iReff, iTemp,         1:5, iFreq)
              cloudcoeff%pcoeff_L_MW(iFreq, iReff, iTemp, LOFFSET_6:(LOFFSET_6+6), 1) = &
                                 pcoeff_mw (iReff, iTemp,         1:7, iFreq)
              cloudcoeff%pcoeff_L_MW(iFreq, iReff, iTemp, LOFFSET_8:(LOFFSET_8+8), 1) = &
                                 pcoeff_mw (iReff, iTemp,         1:9, iFreq)
              cloudcoeff%pcoeff_L_MW(iFreq, iReff, iTemp, LOFFSET_16:(LOFFSET_16+16), 1) = &
                                 pcoeff_mw (iReff, iTemp,         1:17, iFreq)
            end do
          end do
        end do
    end select

    call handle_nc_error( nf90_close( ncid ), &
                          SRC_FILENAME, SRC_FUNCNAME )
    deallocate( freq_mw, temperature, Re_mw, &
                ke_mw, w_mw, g_mw, pcoeff_mw )
                 
  end subroutine nc2cloudcoeff


  !------------------------------------------!
  ! inquire NetCDF file dimensions
  subroutine lut_nc_inq_dim( ncid, &
                         n_freq, &
                         n_temp, &
                         n_Re, &
                         n_Legendre)

    integer, intent(in)  :: ncid
    integer, intent(out) :: n_freq
    integer, intent(out) :: n_temp
    integer, intent(out) :: n_Re
    integer, intent(out) :: n_Legendre

    ! local variables
    character(len=*), parameter :: SRC_FUNCNAME = 'lut_nc_inq_dim'  
    integer :: dimid
    integer :: tmp_n_Legendre

    ! frequency
    call handle_nc_error( nf90_inq_dimid(ncid, "f", dimid),&
                           SRC_FILENAME, SRC_FUNCNAME )
    call handle_nc_error( nf90_inquire_dimension(ncid, dimid, len=n_freq),&
                           SRC_FILENAME, SRC_FUNCNAME )
    ! temperature 
    call handle_nc_error( nf90_inq_dimid(ncid, "T", dimid),& 
                           SRC_FILENAME, SRC_FUNCNAME )
    call handle_nc_error( nf90_inquire_dimension(ncid, dimid, len=n_temp),&
                           SRC_FILENAME, SRC_FUNCNAME )
    ! effective radius
    call handle_nc_error( nf90_inq_dimid(ncid, "Re", dimid),&
                           SRC_FILENAME, SRC_FUNCNAME )
    call handle_nc_error( nf90_inquire_dimension(ncid, dimid, len=n_Re),&
                           SRC_FILENAME, SRC_FUNCNAME )
    ! LegendreOrder
    call handle_nc_error( nf90_inq_dimid(ncid, "LegendreOrder", dimid),&
                           SRC_FILENAME, SRC_FUNCNAME )
    call handle_nc_error( nf90_inquire_dimension(ncid, dimid, len=tmp_n_Legendre),&
                           SRC_FILENAME, SRC_FUNCNAME )
    n_Legendre = tmp_n_legendre -1 ! TODO: this is possibly right

  end subroutine lut_nc_inq_dim

  !------------------------------------------!
  ! inquire NetCDF file dimensions   
  subroutine lut_nc_read_var( ncid,        &
                              freq_mw,     &
                              temperature, &
                              Re_mw,       &
                              ke_mw,       &
                              w_mw,        &
                              g_mw,        &
                              pcoeff_mw)

    integer, intent(in)    :: ncid
    real   , intent(inout) :: freq_mw(:)
    real   , intent(inout) :: temperature(:)
    real   , intent(inout) :: Re_mw(:)
    real   , intent(inout) :: ke_mw(:,:,:)
    real   , intent(inout) :: w_mw (:,:,:)
    real   , intent(inout) :: g_mw (:,:,:)
    real   , intent(inout) :: pcoeff_mw(:,:,:,:)

    ! local variables 
    character(len=*), parameter :: SRC_FUNCNAME = 'lut_nc_read_var'  
    integer :: varid

    integer :: dimids(3)
    integer :: dimlen(3)
    integer :: i
    integer :: status

    ! frequency
    call handle_nc_error( nf90_inq_varid( ncid, 'f', varid ) , &
                           SRC_FILENAME, SRC_FUNCNAME // ' inq f' )
    call handle_nc_error( nf90_get_var(   ncid, varid, freq_mw ), &
                           SRC_FILENAME, SRC_FUNCNAME // ' read f')
    ! Temperature
    call handle_nc_error( nf90_inq_varid( ncid, 'T', varid )  , &
                           SRC_FILENAME, SRC_FUNCNAME  // ' inq T' )
    call handle_nc_error( nf90_get_var(   ncid, varid, temperature), &
                           SRC_FILENAME, SRC_FUNCNAME  // ' read T')
    ! effective radius
    call handle_nc_error( nf90_inq_varid( ncid, 'Re', varid )  , &
                           SRC_FILENAME, SRC_FUNCNAME  // ' inq Re' )
    call handle_nc_error( nf90_get_var(   ncid, varid, Re_mw), &
                           SRC_FILENAME, SRC_FUNCNAME  // ' read Re')
    ! ke_mw
    call handle_nc_error( nf90_inq_varid( ncid, 'MassSigExt', varid )  , &
                           SRC_FILENAME, SRC_FUNCNAME  // ' inq ke' )
    call handle_nc_error( nf90_inquire_variable( ncid, varid, dimids = dimids )  , &
                           SRC_FILENAME, SRC_FUNCNAME  // ' inq ke' )
    do i = 1, 3
      status = nf90_inquire_dimension(ncid, dimIDs(i), len = dimlen(i)) 
      print *, dimlen(i) 
    end do


    call handle_nc_error( nf90_get_var(   ncid, varid, ke_mw), &
                           SRC_FILENAME, SRC_FUNCNAME  // ' read ke')
    ! w_mw
    call handle_nc_error( nf90_inq_varid( ncid, 'w0', varid )  , &
                           SRC_FILENAME, SRC_FUNCNAME  // ' inq w0' )
    call handle_nc_error( nf90_get_var(   ncid, varid, w_mw), &
                           SRC_FILENAME, SRC_FUNCNAME  // ' read w0')
    ! g_mw
    call handle_nc_error( nf90_inq_varid( ncid, 'g', varid )  , &
                           SRC_FILENAME, SRC_FUNCNAME  // ' inq g' )
    call handle_nc_error( nf90_get_var(   ncid, varid, g_mw), &
                           SRC_FILENAME, SRC_FUNCNAME  // ' read g')
    ! Legendre coefficients
    call handle_nc_error( nf90_inq_varid( ncid, 'CPn', varid )  , &
                           SRC_FILENAME, SRC_FUNCNAME  // ' inq CPn' )
    call handle_nc_error( nf90_get_var(   ncid, varid, pcoeff_mw), &
                           SRC_FILENAME, SRC_FUNCNAME  // ' read CPn')
    
  end subroutine lut_nc_read_var
  !------------------------------------------!
  subroutine handle_nc_error( status, src_filename, src_funcname )
    
    integer(kind=4),  intent(in) :: status
    character(len=*), intent(in) :: src_filename
    character(len=*), intent(in) :: src_funcname
    
    if ( status /= NF90_NOERR) then
      print *, 'Error in NetCDF operation in ', src_filename, ', ', src_funcname 
      print *, 'error code ', status
      stop
    end if

  end subroutine handle_nc_error  
   


end module

