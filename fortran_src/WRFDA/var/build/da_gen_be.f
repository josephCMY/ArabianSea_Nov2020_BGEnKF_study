












module da_gen_be

   
   
   
   
   
   
   
   
   


   

   use da_control, only : stdout,vertical_ip, t0,es_beta,es_alpha, &
      es_gamma,kappa,rd_over_rv,rd_over_rv1,t_kelvin, gravity, &
      filename_len,vertical_ip_0, trace_use, trace_use_dull, cv_options, use_rf, do_normalize
   use da_reporting, only : da_error, da_warning, da_message, message
   use da_tools_serial, only : da_get_unit, da_free_unit, da_array_print
   use da_lapack, only : dsyev
   use da_wavelet, only: lf,namw,nb,nij,ws,wsd

   implicit none

   real, parameter :: base_pres = 100000.0 
               
contains



subroutine da_trace_entry(name, message, messages, maxnocalls)    
   implicit none
   character (len=*),           intent(in) :: name         
   character (len=*), optional, intent(in) :: message      
   character (len=*), optional, intent(in) :: messages(:)  
   integer, optional,           intent(in) :: maxnocalls   
end subroutine da_trace_entry

subroutine da_trace_exit(name, message, messages, maxnocalls)
   implicit none
   character (len=*), intent(in)           :: name         
   character (len=*), optional, intent(in) :: message      
   character (len=*), optional, intent(in) :: messages(:)  
   integer, optional, intent(in)           :: maxnocalls  
end subroutine da_trace_exit

subroutine da_create_bins(ni, nj, nk, bin_type, num_bins, num_bins2d, bin, bin2d, &
   lat_min, lat_max, binwidth_lat, hgt_min, hgt_max, binwidth_hgt, latitude, height)

   !----------------------------------------------------------------------------
   !
   ! Purpose: To create the bins for calculation of statistics.
   !
   ! Input:
   ! ni, nj, nk   Dimensions
   ! bin_type     0: No binning; 
   !              1: bins for X-direction mean; 
   !              2: bins for each of binwidth_lat/binwidth_hgt.  
   !              3: bins for each of binwidth_lat/nk.  
   !              4: num_hor_bins horizontal bins /nk.  
   !              5: Average over all horizontal points (nk bins for 3D fields)
   !              6: Average over all points (only 1 bin).
   ! Optional for bin_type = 2:
   ! lat_min, lat_max Minimum/maximum latitude ranges for bin_type = 2
   ! binwidth_lat interval between bins for latitude in degree for bin_type = 2
   ! binwidth_hgt interval between bins for height in meter for bin_type = 2
   ! num_bins_hgt the number of height bins for bin_type = 2
   ! latitude     3d field latitude in degree for bin_type = 2
   ! height       3d field height in meter for bin_type = 2
   !
   ! Output:
   ! num_bins,num_bins2d ---- the number of bins for 3d and 2d fields
   ! bin     ,     bin2d ---- Assigned bin to a gridpoints for 3d and 2d fields
   !
   !----------------------------------------------------------------------------

   implicit none

   integer,        intent(in)  :: ni, nj, nk          ! Dimensions read in.
   integer,        intent(in)  :: bin_type            ! Type of bin to average over
   integer,        intent(out) :: num_bins            ! Number of bins.
   integer,        intent(out) :: num_bins2d          ! Number of bins for 2D fields
   integer,        intent(out) :: bin(1:ni,1:nj,1:nk) ! Bin at each 3D point
   integer,        intent(out) :: bin2d(1:ni,1:nj)    ! Bin at each 2D point

   real, optional, intent(in)  :: lat_min, lat_max    ! Used if bin_type = 2 (deg)
   real, optional, intent(in)  :: binwidth_lat        ! Used if bin_type = 2 (deg)
   real, optional, intent(in)  :: hgt_min, hgt_max    ! Used if bin_type = 2 (deg)
   real, optional, intent(in)  :: binwidth_hgt        ! Used if bin_type = 2 (m).
   real, optional, intent(in)  :: latitude(1:ni,1:nj) ! Latitude (degrees).
   real, optional, intent(in)  :: height(1:ni,1:nj,1:nk)     ! Height (m).

   integer           :: b, i, j, k                 ! Loop counters.
   integer           :: count                      ! Counter
   integer           :: num_bins_lat               ! Used if bin_type = 2.
   integer           :: num_bins_hgt               ! Used if bin_type = 2.
   integer           :: bin_lat                    ! Latitude bin.
   integer           :: bin_hgt                    ! Height bin.
   integer           :: num_bins_i, num_bins_j     ! Used if bin_type = 4.
   integer           :: nii, njj                   ! Used if bin_type = 4.
   integer           :: bin_i(1:ni), bin_j(1:nj)   ! Used if bin_type = 4.
   real, allocatable :: binstart_lat(:)            ! Used if bin_type = 2 (deg)
   real, allocatable :: binstart_hgt(:)            ! Used if bin_type = 2 (deg)

   if (trace_use) call da_trace_entry("da_create_bins")

   if (bin_type == 0) then         ! No averaging in space

      num_bins = nk * nj * ni
      num_bins2d = nj * ni    ! Equals number of horizontal points.

      count = 1
      do k = 1, nk
         do j = 1, nj
            do i = 1, ni
               bin(i,j,k) = count
               count = count + 1
            end do
         end do
      end do
      bin2d(:,:) = bin(:,:,1)

   else if (bin_type == 1) then    ! Average over x-direction.

      num_bins = nj * nk
      num_bins2d = nj

      count = 1
      do k = 1, nk
         do j = 1, nj
            bin(1:ni,j,k) = count
            count = count + 1
         end do
      end do
      bin2d(:,:) = bin(:,:,1)

   else if (bin_type == 2) then    ! Global latitude/height bins:

      ! Setup latitude bins:
      write(unit=stdout,fmt='(/a,f12.5)')'   Minimum latitude = ', lat_min
      write(unit=stdout,fmt='(a,f12.5)')'    Maximum latitude = ', lat_max
      write(unit=stdout,fmt='(a,f12.5)') &
         '    Latitude bin width = ', binwidth_lat
      num_bins_lat = nint((lat_max - lat_min) / binwidth_lat)
      write(unit=stdout,fmt='(a,i8)') &
         '    Number of latitude bins = ', num_bins_lat
   
      allocate(binstart_lat(1:num_bins_lat))
      do b = 1, num_bins_lat ! Assume south to north (as in WRF).
         binstart_lat(b) = lat_min + real(b-1) * binwidth_lat
      end do

      ! Setup height bins:
      write(unit=stdout,fmt='(/a,f12.5)')'    Minimum height = ', hgt_min
      write(unit=stdout,fmt='(a,f12.5)')'    Maximum height = ', hgt_max
      write(unit=stdout,fmt='(a,f12.5)')'    Height bin width = ', binwidth_hgt
      num_bins_hgt = nint((hgt_max - hgt_min) / binwidth_hgt)
      write(unit=stdout,fmt='(a,i8)') &
         '    Number of height bins = ', num_bins_hgt

      allocate(binstart_hgt(1:num_bins_hgt))
      do b = 1, num_bins_hgt
         binstart_hgt(b) = hgt_min + real(b-1) * binwidth_hgt
      end do

      num_bins = num_bins_lat * num_bins_hgt
      num_bins2d = num_bins_lat

      ! Select height bins:
      do j = 1, nj
         do i = 1, ni
            do k = 1, nk
               if (height(i,j,k) < binstart_hgt(1)) then 
                  bin_hgt = 1 ! In first bin.
               else if (height(i,j,k) >= binstart_hgt(num_bins_hgt)) then
                  bin_hgt = num_bins_hgt ! In final bin.
               else 
                  do b = 1, num_bins_hgt-1
                     if (height(i,j,k) >= binstart_hgt(b) .and. &
                          height(i,j,k) <  binstart_hgt(b+1)) then
                        bin_hgt = b
                        exit
                     end if
                  end do
               end if

               ! Select latitude bin that point falls in:
               if (k == 1) then
                  do b = 1, num_bins_lat-1
                     if (latitude(i,j) >= binstart_lat(b) .and. &
                        latitude(i,j) < binstart_lat(b+1)) then
                        bin_lat = b
                        exit
                     end if
                  end do
                  if (latitude(i,j) >= binstart_lat(num_bins_lat)) then
                     ! In final bin.
                     bin_lat = num_bins_lat
                  end if
                  bin2d(i,j) = bin_lat
               end if
               bin(i,j,k) = bin_lat + num_bins_lat * (bin_hgt - 1)
            end do
         end do
      end do

      deallocate(binstart_lat)
      deallocate(binstart_hgt)

   else if (bin_type == 3) then    ! Latitude/nk bins:

      ! Setup latitude bins:
      write(unit=stdout,fmt='(/a,f12.5)')'   Minimum latitude = ', lat_min
      write(unit=stdout,fmt='(a,f12.5)')'    Maximum latitude = ', lat_max
      write(unit=stdout,fmt='(a,f12.5)')'    Latitude bin width = ',binwidth_lat
      num_bins_lat = nint((lat_max - lat_min) / binwidth_lat)
      write(unit=stdout,fmt='(a,i8)') &
         '    Number of latitude bins = ', num_bins_lat
   
      allocate(binstart_lat(1:num_bins_lat))
      do b = 1, num_bins_lat ! Assume south to north (as in WRF).
         binstart_lat(b) = lat_min + real(b-1) * binwidth_lat
      end do

      num_bins = num_bins_lat * nk
      num_bins2d = num_bins_lat

      ! Select bins:
      do j = 1, nj
         do i = 1, ni
            do k = 1, nk
               ! Select latitude bin that point falls in:
               if (k == 1) then
                  do b = 1, num_bins_lat-1
                     if (latitude(i,j) >= binstart_lat(b) .and. &
                        latitude(i,j) < binstart_lat(b+1)) then
                        bin_lat = b
                        exit
                     end if
                  end do
                  if (latitude(i,j) >= binstart_lat(num_bins_lat)) then
                     ! In final bin.
                     bin_lat = num_bins_lat
                  end if
                  bin2d(i,j) = bin_lat
               end if
               bin(i,j,k) = bin_lat + num_bins_lat * (k - 1)
            end do
         end do
      end do

      deallocate(binstart_lat)

   else if (bin_type == 4) then    ! binwidth_lat/nk bins:

      ! Setup horizontal bins:
      write(unit=stdout,fmt='(/a,f12.5)') &
         '   Number of grid-cells to average over = ', binwidth_lat
      ! use binwidth_lat, but actually an integer number of points.
 
      num_bins_j = int(real(nj) / real(binwidth_lat))
      njj = int(binwidth_lat) * num_bins_j
      do j = 1, njj
         bin_j(j) = 1 + int(real(j-1) / binwidth_lat)
      end do
      if (nj > njj) bin_j(njj+1:nj) = bin_j(njj)

      num_bins_i = int(real(ni) / real(binwidth_lat))
      nii = int(binwidth_lat) * num_bins_i
      do i = 1, nii
         bin_i(i) = 1 + int(real(i-1) / binwidth_lat)
      end do
      if (ni > nii) bin_i(nii+1:ni) = bin_i(nii)

      num_bins2d = num_bins_i * num_bins_j
      num_bins = num_bins2d * nk

      do j = 1, nj
         do i = 1, ni
            bin2d(i,j) = bin_i(i) + (bin_j(j) - 1) * num_bins_i
            do k = 1, nk
               bin(i,j,k) = bin2d(i,j) + (k - 1) * num_bins2d
            end do
         end do
      end do

   else if (bin_type == 5) then    ! Average over all horizontal points.

      num_bins = nk
      num_bins2d = 1

      do k = 1, nk
         bin(:,:,k) = k
      end do
      bin2d(:,:) = 1

   else if (bin_type == 6) then    ! Average over all points.

      num_bins = 1
      num_bins2d = 1
      bin(:,:,:) = 1
      bin2d(:,:) = 1
   end if

   if (trace_use) call da_trace_exit("da_create_bins")

end subroutine da_create_bins


subroutine da_filter_regcoeffs(ni, nj, nk, num_bins, num_bins2d, num_passes, &
   rf_scale, bin, regcoeff1, regcoeff2, regcoeff3)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: ni, nj, nk                  ! Grid dimensions.
   integer, intent(in)    :: num_bins                    ! Number of bins for 3D coeffs.
   integer, intent(in)    :: num_bins2d                  ! Number of bins for 2D coeffs.
   integer, intent(in)    :: num_passes                  ! Number of passes for RF.
   real,    intent(in)    :: rf_scale                    ! Smoothing scale of RF.
   integer, intent(in)    :: bin(1:ni,1:nj,1:nk)         ! Bin assigned to each point.
   real,    intent(inout) :: regcoeff1(1:num_bins)       ! psi/chi regression cooefficient.
   real,    intent(inout) :: regcoeff2(1:nk,1:num_bins2d)! psi/ps regression cooefficient.
   real,    intent(inout) :: regcoeff3(1:nk,1:nk,1:num_bins2d) ! psi/T regression cooefficient.

   integer :: i, j, k              ! Loop counters.
   integer :: b                    ! Bin index.
   real*8  :: field(1:ni,1:nj)     ! Field for recursive filter.

   if (trace_use) call da_trace_entry("da_filter_regcoeffs")

   !----------------------------------------------------------------------------
   ! [1] Filter psi/chi coefficient:
   !----------------------------------------------------------------------------

   do k = 1, nk
      ! Prepare field for filtering:
      do j = 1, nj
         do i = 1, ni
            b = bin(i,j,k)
            field(i,j) = regcoeff1(b)
         end do
      end do

      call da_perform_2drf(ni, nj, num_passes, rf_scale, field)

      do j = 1, nj
         do i = 1, ni
            b = bin(i,j,k)
            regcoeff1(b) = field(i,j)
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_filter_regcoeffs")

end subroutine da_filter_regcoeffs


subroutine da_get_field( input_file, var, field_dims, dim1, dim2, dim3,k,field)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

!     NetCDF-3.
!
! netcdf version 3 fortran interface:
!

!
! external netcdf data types:
!
      integer nf_byte
      integer nf_int1
      integer nf_char
      integer nf_short
      integer nf_int2
      integer nf_int
      integer nf_float
      integer nf_real
      integer nf_double
      integer nf_ubyte
      integer nf_ushort
      integer nf_uint
      integer nf_int64
      integer nf_uint64

      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)

!
! default fill values:
!
      integer           nf_fill_byte
      integer           nf_fill_int1
      integer           nf_fill_char
      integer           nf_fill_short
      integer           nf_fill_int2
      integer           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double

      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690e+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690d+36)

!
! mode flags for opening and creating a netcdf dataset:
!
      integer nf_nowrite
      integer nf_write
      integer nf_clobber
      integer nf_noclobber
      integer nf_fill
      integer nf_nofill
      integer nf_lock
      integer nf_share
      integer nf_64bit_offset
      integer nf_64bit_data
      integer nf_cdf5
      integer nf_sizehint_default
      integer nf_align_chunk
      integer nf_format_classic
      integer nf_format_64bit
      integer nf_format_64bit_offset
      integer nf_format_64bit_data
      integer nf_format_cdf5
      integer nf_diskless
      integer nf_mmap

      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)

!
! size argument for defining an unlimited dimension:
!
      integer nf_unlimited
      parameter (nf_unlimited = 0)

!
! global attribute id:
!
      integer nf_global
      parameter (nf_global = 0)

!
! implementation limits:
!
      integer nf_max_dims
      integer nf_max_attrs
      integer nf_max_vars
      integer nf_max_name
      integer nf_max_var_dims

      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)

!
! error codes:
!
      integer nf_noerr
      integer nf_ebadid
      integer nf_eexist
      integer nf_einval
      integer nf_eperm
      integer nf_enotindefine
      integer nf_eindefine
      integer nf_einvalcoords
      integer nf_emaxdims
      integer nf_enameinuse
      integer nf_enotatt
      integer nf_emaxatts
      integer nf_ebadtype
      integer nf_ebaddim
      integer nf_eunlimpos
      integer nf_emaxvars
      integer nf_enotvar
      integer nf_eglobal
      integer nf_enotnc
      integer nf_ests
      integer nf_emaxname
      integer nf_eunlimit
      integer nf_enorecvars
      integer nf_echar
      integer nf_eedge
      integer nf_estride
      integer nf_ebadname
      integer nf_erange
      integer nf_enomem
      integer nf_evarsize
      integer nf_edimsize
      integer nf_etrunc

      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
!
! error handling modes:
!
      integer  nf_fatal
      integer nf_verbose

      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)

!
! miscellaneous routines:
!
      character*80   nf_inq_libvers
      external       nf_inq_libvers

      character*80   nf_strerror
!                         (integer             ncerr)
      external       nf_strerror

      logical        nf_issyserr
!                         (integer             ncerr)
      external       nf_issyserr

!
! control routines:
!
      integer         nf_inq_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_inq_base_pe

      integer         nf_set_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_set_base_pe

      integer         nf_create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             ncid)
      external        nf_create

      integer         nf__create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create

      integer         nf__create_mp
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create_mp

      integer         nf_open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             ncid)
      external        nf_open

      integer         nf__open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open

      integer         nf__open_mp
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open_mp

      integer         nf_set_fill
!                         (integer             ncid,
!                          integer             fillmode,
!                          integer             old_mode)
      external        nf_set_fill

      integer         nf_set_default_format
!                          (integer             format,
!                          integer             old_format)
      external        nf_set_default_format

      integer         nf_redef
!                         (integer             ncid)
      external        nf_redef

      integer         nf_enddef
!                         (integer             ncid)
      external        nf_enddef

      integer         nf__enddef
!                         (integer             ncid,
!                          integer             h_minfree,
!                          integer             v_align,
!                          integer             v_minfree,
!                          integer             r_align)
      external        nf__enddef

      integer         nf_sync
!                         (integer             ncid)
      external        nf_sync

      integer         nf_abort
!                         (integer             ncid)
      external        nf_abort

      integer         nf_close
!                         (integer             ncid)
      external        nf_close

      integer         nf_delete
!                         (character*(*)       ncid)
      external        nf_delete

!
! general inquiry routines:
!

      integer         nf_inq
!                         (integer             ncid,
!                          integer             ndims,
!                          integer             nvars,
!                          integer             ngatts,
!                          integer             unlimdimid)
      external        nf_inq

! new inquire path

      integer nf_inq_path
      external nf_inq_path

      integer         nf_inq_ndims
!                         (integer             ncid,
!                          integer             ndims)
      external        nf_inq_ndims

      integer         nf_inq_nvars
!                         (integer             ncid,
!                          integer             nvars)
      external        nf_inq_nvars

      integer         nf_inq_natts
!                         (integer             ncid,
!                          integer             ngatts)
      external        nf_inq_natts

      integer         nf_inq_unlimdim
!                         (integer             ncid,
!                          integer             unlimdimid)
      external        nf_inq_unlimdim

      integer         nf_inq_format
!                         (integer             ncid,
!                          integer             format)
      external        nf_inq_format

!
! dimension routines:
!

      integer         nf_def_dim
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             len,
!                          integer             dimid)
      external        nf_def_dim

      integer         nf_inq_dimid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             dimid)
      external        nf_inq_dimid

      integer         nf_inq_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_dim

      integer         nf_inq_dimname
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_inq_dimname

      integer         nf_inq_dimlen
!                         (integer             ncid,
!                          integer             dimid,
!                          integer             len)
      external        nf_inq_dimlen

      integer         nf_rename_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_rename_dim

!
! general attribute routines:
!

      integer         nf_inq_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len)
      external        nf_inq_att

      integer         nf_inq_attid
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             attnum)
      external        nf_inq_attid

      integer         nf_inq_atttype
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype)
      external        nf_inq_atttype

      integer         nf_inq_attlen
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_attlen

      integer         nf_inq_attname
!                         (integer             ncid,
!                          integer             varid,
!                          integer             attnum,
!                          character(*)        name)
      external        nf_inq_attname

      integer         nf_copy_att
!                         (integer             ncid_in,
!                          integer             varid_in,
!                          character(*)        name,
!                          integer             ncid_out,
!                          integer             varid_out)
      external        nf_copy_att

      integer         nf_rename_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        curname,
!                          character(*)        newname)
      external        nf_rename_att

      integer         nf_del_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_del_att

!
! attribute put/get routines:
!

      integer         nf_put_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len,
!                          character(*)        text)
      external        nf_put_att_text

      integer         nf_get_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          character(*)        text)
      external        nf_get_att_text

      integer         nf_put_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int1_t           i1vals(1))
      external        nf_put_att_int1

      integer         nf_get_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int1_t           i1vals(1))
      external        nf_get_att_int1

      integer         nf_put_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int2_t           i2vals(1))
      external        nf_put_att_int2

      integer         nf_get_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int2_t           i2vals(1))
      external        nf_get_att_int2

      integer         nf_put_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          integer             ivals(1))
      external        nf_put_att_int

      integer         nf_get_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             ivals(1))
      external        nf_get_att_int

      integer         nf_put_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int8_t           i8vals(1))
      external        nf_put_att_int64

      integer         nf_get_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int8_t           i8vals(1))
      external        nf_get_att_int64

      integer         nf_put_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          real                rvals(1))
      external        nf_put_att_real

      integer         nf_get_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          real                rvals(1))
      external        nf_get_att_real

      integer         nf_put_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          double              dvals(1))
      external        nf_put_att_double

      integer         nf_get_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          double              dvals(1))
      external        nf_get_att_double

!
! general variable routines:
!

      integer         nf_def_var
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             varid)
      external        nf_def_var

      integer         nf_inq_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             natts)
      external        nf_inq_var

      integer         nf_inq_varid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             varid)
      external        nf_inq_varid

      integer         nf_inq_varname
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_inq_varname

      integer         nf_inq_vartype
!                         (integer             ncid,
!                          integer             varid,
!                          integer             xtype)
      external        nf_inq_vartype

      integer         nf_inq_varndims
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ndims)
      external        nf_inq_varndims

      integer         nf_inq_vardimid
!                         (integer             ncid,
!                          integer             varid,
!                          integer             dimids(1))
      external        nf_inq_vardimid

      integer         nf_inq_varnatts
!                         (integer             ncid,
!                          integer             varid,
!                          integer             natts)
      external        nf_inq_varnatts

      integer         nf_rename_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_rename_var

      integer         nf_copy_var
!                         (integer             ncid_in,
!                          integer             varid,
!                          integer             ncid_out)
      external        nf_copy_var

!
! entire variable put/get routines:
!

      integer         nf_put_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_put_var_text

      integer         nf_get_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_get_var_text

      integer         nf_put_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_put_var_int1

      integer         nf_get_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_get_var_int1

      integer         nf_put_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_put_var_int2

      integer         nf_get_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_get_var_int2

      integer         nf_put_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_put_var_int

      integer         nf_get_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_get_var_int

      integer         nf_put_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_put_var_real

      integer         nf_get_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_get_var_real

      integer         nf_put_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_put_var_double

      integer         nf_get_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_get_var_double

!
! single variable put/get routines:
!

      integer         nf_put_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_put_var1_text

      integer         nf_get_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_get_var1_text

      integer         nf_put_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_put_var1_int1

      integer         nf_get_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_get_var1_int1

      integer         nf_put_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_put_var1_int2

      integer         nf_get_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_get_var1_int2

      integer         nf_put_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_put_var1_int

      integer         nf_get_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_get_var1_int

      integer         nf_put_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_put_var1_real

      integer         nf_get_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_get_var1_real

      integer         nf_put_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_put_var1_double

      integer         nf_get_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_get_var1_double

!
! variable array put/get routines:
!

      integer         nf_put_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_put_vara_text

      integer         nf_get_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_get_vara_text

      integer         nf_put_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vara_int1

      integer         nf_get_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vara_int1

      integer         nf_put_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vara_int2

      integer         nf_get_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vara_int2

      integer         nf_put_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_put_vara_int

      integer         nf_get_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_get_vara_int

      integer         nf_put_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_put_vara_real

      integer         nf_get_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_get_vara_real

      integer         nf_put_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vara_double

      integer         nf_get_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vara_double

!
! strided variable put/get routines:
!

      integer         nf_put_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_put_vars_text

      integer         nf_get_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_get_vars_text

      integer         nf_put_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vars_int1

      integer         nf_get_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vars_int1

      integer         nf_put_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vars_int2

      integer         nf_get_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vars_int2

      integer         nf_put_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_put_vars_int

      integer         nf_get_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_get_vars_int

      integer         nf_put_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_put_vars_real

      integer         nf_get_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_get_vars_real

      integer         nf_put_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vars_double

      integer         nf_get_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vars_double

!
! mapped variable put/get routines:
!

      integer         nf_put_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_put_varm_text

      integer         nf_get_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_get_varm_text

      integer         nf_put_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_varm_int1

      integer         nf_get_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_varm_int1

      integer         nf_put_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_varm_int2

      integer         nf_get_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_varm_int2

      integer         nf_put_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_put_varm_int

      integer         nf_get_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_get_varm_int

      integer         nf_put_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_put_varm_real

      integer         nf_get_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_get_varm_real

      integer         nf_put_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_put_varm_double

      integer         nf_get_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_get_varm_double

!     64-bit int functions.
      integer nf_put_var1_int64
      external nf_put_var1_int64
      integer nf_put_vara_int64
      external nf_put_vara_int64
      integer nf_put_vars_int64
      external nf_put_vars_int64
      integer nf_put_varm_int64
      external nf_put_varm_int64
      integer nf_put_var_int64
      external nf_put_var_int64
      integer nf_get_var1_int64
      external nf_get_var1_int64
      integer nf_get_vara_int64
      external nf_get_vara_int64
      integer nf_get_vars_int64
      external nf_get_vars_int64
      integer nf_get_varm_int64
      external nf_get_varm_int64
      integer nf_get_var_int64
      external nf_get_var_int64


!     NetCDF-4.
!     This is part of netCDF-4. Copyright 2006, UCAR, See COPYRIGHT
!     file for distribution information.

!     Netcdf version 4 fortran interface.

!     $Id: netcdf4.inc,v 1.28 2010/05/25 13:53:02 ed Exp $

!     New netCDF-4 types.
      integer nf_string
      integer nf_vlen
      integer nf_opaque
      integer nf_enum
      integer nf_compound

      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)

!     New netCDF-4 fill values.
      integer           nf_fill_ubyte
      integer           nf_fill_ushort
!      real              nf_fill_uint
!      real              nf_fill_int64
!      real              nf_fill_uint64
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)

!     New constants.
      integer nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)

      integer nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)

      integer nf_netcdf4
      parameter (nf_netcdf4 = 4096)

      integer nf_classic_model
      parameter (nf_classic_model = 256)

      integer nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)

      integer nf_endian_native
      parameter (nf_endian_native = 0)
      integer nf_endian_little
      parameter (nf_endian_little = 1)
      integer nf_endian_big
      parameter (nf_endian_big = 2)

!     For NF_DEF_VAR_CHUNKING
      integer nf_chunked
      parameter (nf_chunked = 0)
      integer nf_contiguous
      parameter (nf_contiguous = 1)
      integer nf_compact
      parameter (nf_compact = 2)

!     For NF_DEF_VAR_FLETCHER32
      integer nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer nf_fletcher32
      parameter (nf_fletcher32 = 1)

!     For NF_DEF_VAR_DEFLATE
      integer nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer nf_shuffle
      parameter (nf_shuffle = 1)

!     For NF_DEF_VAR_SZIP
      integer nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)

!     For parallel I/O.
      integer nf_mpiio      
      parameter (nf_mpiio = 8192)
      integer nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer nf_pnetcdf
      parameter (nf_pnetcdf = 32768)

!     For NF_VAR_PAR_ACCESS.
      integer nf_independent
      parameter (nf_independent = 0)
      integer nf_collective
      parameter (nf_collective = 1)

!     New error codes.
      integer nf_ehdferr        ! Error at 1 layer. 
      parameter (nf_ehdferr = -101)
      integer nf_ecantread      ! Can't read. 
      parameter (nf_ecantread = -102)
      integer nf_ecantwrite     ! Can't write. 
      parameter (nf_ecantwrite = -103)
      integer nf_ecantcreate    ! Can't create. 
      parameter (nf_ecantcreate = -104)
      integer nf_efilemeta      ! Problem with file metadata. 
      parameter (nf_efilemeta = -105)
      integer nf_edimmeta       ! Problem with dimension metadata. 
      parameter (nf_edimmeta = -106)
      integer nf_eattmeta       ! Problem with attribute metadata. 
      parameter (nf_eattmeta = -107)
      integer nf_evarmeta       ! Problem with variable metadata. 
      parameter (nf_evarmeta = -108)
      integer nf_enocompound    ! Not a compound type. 
      parameter (nf_enocompound = -109)
      integer nf_eattexists     ! Attribute already exists. 
      parameter (nf_eattexists = -110)
      integer nf_enotnc4        ! Attempting netcdf-4 operation on netcdf-3 file.   
      parameter (nf_enotnc4 = -111)
      integer nf_estrictnc3     ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
      parameter (nf_estrictnc3 = -112)
      integer nf_enotnc3        ! Attempting netcdf-3 operation on netcdf-4 file.   
      parameter (nf_enotnc3 = -113)
      integer nf_enopar         ! Parallel operation on file opened for non-parallel access.   
      parameter (nf_enopar = -114)
      integer nf_eparinit       ! Error initializing for parallel access.   
      parameter (nf_eparinit = -115)
      integer nf_ebadgrpid      ! Bad group ID.   
      parameter (nf_ebadgrpid = -116)
      integer nf_ebadtypid      ! Bad type ID.   
      parameter (nf_ebadtypid = -117)
      integer nf_etypdefined    ! Type has already been defined and may not be edited. 
      parameter (nf_etypdefined = -118)
      integer nf_ebadfield      ! Bad field ID.   
      parameter (nf_ebadfield = -119)
      integer nf_ebadclass      ! Bad class.   
      parameter (nf_ebadclass = -120)
      integer nf_emaptype       ! Mapped access for atomic types only.   
      parameter (nf_emaptype = -121)
      integer nf_elatefill      ! Attempt to define fill value when data already exists. 
      parameter (nf_elatefill = -122)
      integer nf_elatedef       ! Attempt to define var properties, like deflate, after enddef. 
      parameter (nf_elatedef = -123)
      integer nf_edimscale      ! Probem with 1 dimscales. 
      parameter (nf_edimscale = -124)
      integer nf_enogrp       ! No group found.
      parameter (nf_enogrp = -125)


!     New functions.

!     Parallel I/O.
      integer nf_create_par
      external nf_create_par

      integer nf_open_par
      external nf_open_par

      integer nf_var_par_access
      external nf_var_par_access

!     Functions to handle groups.
      integer nf_inq_ncid
      external nf_inq_ncid

      integer nf_inq_grps
      external nf_inq_grps

      integer nf_inq_grpname
      external nf_inq_grpname

      integer nf_inq_grpname_full
      external nf_inq_grpname_full

      integer nf_inq_grpname_len
      external nf_inq_grpname_len

      integer nf_inq_grp_parent
      external nf_inq_grp_parent

      integer nf_inq_grp_ncid
      external nf_inq_grp_ncid

      integer nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid

      integer nf_inq_varids
      external nf_inq_varids

      integer nf_inq_dimids
      external nf_inq_dimids

      integer nf_def_grp
      external nf_def_grp

!     New rename grp function

      integer nf_rename_grp
      external nf_rename_grp

!     New options for netCDF variables.
      integer nf_def_var_deflate
      external nf_def_var_deflate

      integer nf_inq_var_deflate
      external nf_inq_var_deflate

      integer nf_def_var_szip
      external nf_def_var_szip

      integer nf_inq_var_szip
      external nf_inq_var_szip

      integer nf_def_var_fletcher32
      external nf_def_var_fletcher32

      integer nf_inq_var_fletcher32
      external nf_inq_var_fletcher32

      integer nf_def_var_chunking
      external nf_def_var_chunking

      integer nf_inq_var_chunking
      external nf_inq_var_chunking

      integer nf_def_var_fill
      external nf_def_var_fill

      integer nf_inq_var_fill
      external nf_inq_var_fill

      integer nf_def_var_endian
      external nf_def_var_endian

      integer nf_inq_var_endian
      external nf_inq_var_endian

      integer nf_def_var_filter
      external nf_def_var_filter

      integer nf_inq_var_filter
      external nf_inq_var_filter

!     User defined types.
      integer nf_inq_typeids
      external nf_inq_typeids

      integer nf_inq_typeid
      external nf_inq_typeid

      integer nf_inq_type
      external nf_inq_type

      integer nf_inq_user_type
      external nf_inq_user_type

!     User defined types - compound types.
      integer nf_def_compound
      external nf_def_compound

      integer nf_insert_compound
      external nf_insert_compound

      integer nf_insert_array_compound
      external nf_insert_array_compound

      integer nf_inq_compound
      external nf_inq_compound

      integer nf_inq_compound_name
      external nf_inq_compound_name

      integer nf_inq_compound_size
      external nf_inq_compound_size

      integer nf_inq_compound_nfields
      external nf_inq_compound_nfields

      integer nf_inq_compound_field
      external nf_inq_compound_field

      integer nf_inq_compound_fieldname
      external nf_inq_compound_fieldname

      integer nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex

      integer nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset

      integer nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype

      integer nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims

      integer nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes

!     User defined types - variable length arrays.
      integer nf_def_vlen
      external nf_def_vlen

      integer nf_inq_vlen
      external nf_inq_vlen

      integer nf_free_vlen
      external nf_free_vlen

!     User defined types - enums.
      integer nf_def_enum
      external nf_def_enum

      integer nf_insert_enum
      external nf_insert_enum

      integer nf_inq_enum
      external nf_inq_enum

      integer nf_inq_enum_member
      external nf_inq_enum_member

      integer nf_inq_enum_ident
      external nf_inq_enum_ident

!     User defined types - opaque.
      integer nf_def_opaque
      external nf_def_opaque

      integer nf_inq_opaque
      external nf_inq_opaque

!     Write and read attributes of any type, including user defined
!     types.
      integer nf_put_att
      external nf_put_att
      integer nf_get_att
      external nf_get_att

!     Write and read variables of any type, including user defined
!     types.
      integer nf_put_var
      external nf_put_var
      integer nf_put_var1
      external nf_put_var1
      integer nf_put_vara
      external nf_put_vara
      integer nf_put_vars
      external nf_put_vars
      integer nf_get_var
      external nf_get_var
      integer nf_get_var1
      external nf_get_var1
      integer nf_get_vara
      external nf_get_vara
      integer nf_get_vars
      external nf_get_vars

!     For helping F77 users with VLENs.
      integer nf_get_vlen_element
      external nf_get_vlen_element
      integer nf_put_vlen_element
      external nf_put_vlen_element

!     For dealing with file level chunk cache.
      integer nf_set_chunk_cache
      external nf_set_chunk_cache
      integer nf_get_chunk_cache
      external nf_get_chunk_cache

!     For dealing with per variable chunk cache.
      integer nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer nf_get_var_chunk_cache
      external nf_get_var_chunk_cache

!     NetCDF-2.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!      
! functions in the fortran interface
!
      integer nccre
      integer ncopn
      integer ncddef
      integer ncdid
      integer ncvdef
      integer ncvid
      integer nctlen
      integer ncsfil

      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil


      integer ncrdwr
      integer nccreat
      integer ncexcl
      integer ncindef
      integer ncnsync
      integer nchsync
      integer ncndirty
      integer nchdirty
      integer nclink
      integer ncnowrit
      integer ncwrite
      integer ncclob
      integer ncnoclob
      integer ncglobal
      integer ncfill
      integer ncnofill
      integer maxncop
      integer maxncdim
      integer maxncatt
      integer maxncvar
      integer maxncnam
      integer maxvdims
      integer ncnoerr
      integer ncebadid
      integer ncenfile
      integer nceexist
      integer nceinval
      integer nceperm
      integer ncenotin
      integer nceindef
      integer ncecoord
      integer ncemaxds
      integer ncename
      integer ncenoatt
      integer ncemaxat
      integer ncebadty
      integer ncebadd
      integer ncests
      integer nceunlim
      integer ncemaxvs
      integer ncenotvr
      integer nceglob
      integer ncenotnc
      integer ncfoobar
      integer ncsyserr
      integer ncfatal
      integer ncverbos
      integer ncentool


!
! netcdf data types:
!
      integer ncbyte
      integer ncchar
      integer ncshort
      integer nclong
      integer ncfloat
      integer ncdouble

      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)

!     
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!     

!     read/write, 0 => readonly 
      parameter(ncrdwr = 1)
!     in create phase, cleared by ncendef 
      parameter(nccreat = 2)
!     on create destroy existing file 
      parameter(ncexcl = 4)
!     in define mode, cleared by ncendef 
      parameter(ncindef = 8)
!     synchronise numrecs on change (x'10')
      parameter(ncnsync = 16)
!     synchronise whole header on change (x'20')
      parameter(nchsync = 32)
!     numrecs has changed (x'40')
      parameter(ncndirty = 64)  
!     header info has changed (x'80')
      parameter(nchdirty = 128)
!     prefill vars on endef and increase of record, the default behavior
      parameter(ncfill = 0)
!     do not fill vars on endef and increase of record (x'100')
      parameter(ncnofill = 256)
!     isa link (x'8000')
      parameter(nclink = 32768)

!     
!     'mode' arguments for nccreate and ncopen
!     
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)

!     
!     'size' argument to ncdimdef for an unlimited dimension
!     
      integer ncunlim
      parameter(ncunlim = 0)

!     
!     attribute id to put/get a global attribute
!     
      parameter(ncglobal  = 0)

!     
!     advisory maximums:
!     
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
!     not enforced 
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)

!     
!     global netcdf error status variable
!     initialized in error.c
!     

!     no error 
      parameter(ncnoerr = nf_noerr)
!     not a netcdf id 
      parameter(ncebadid = nf_ebadid)
!     too many netcdfs open 
      parameter(ncenfile = -31)   ! nc_syserr
!     netcdf file exists && ncnoclob
      parameter(nceexist = nf_eexist)
!     invalid argument 
      parameter(nceinval = nf_einval)
!     write to read only 
      parameter(nceperm = nf_eperm)
!     operation not allowed in data mode 
      parameter(ncenotin = nf_enotindefine )   
!     operation not allowed in define mode 
      parameter(nceindef = nf_eindefine)   
!     coordinates out of domain 
      parameter(ncecoord = nf_einvalcoords)
!     maxncdims exceeded 
      parameter(ncemaxds = nf_emaxdims)
!     string match to name in use 
      parameter(ncename = nf_enameinuse)   
!     attribute not found 
      parameter(ncenoatt = nf_enotatt)
!     maxncattrs exceeded 
      parameter(ncemaxat = nf_emaxatts)
!     not a netcdf data type 
      parameter(ncebadty = nf_ebadtype)
!     invalid dimension id 
      parameter(ncebadd = nf_ebaddim)
!     ncunlimited in the wrong index 
      parameter(nceunlim = nf_eunlimpos)
!     maxncvars exceeded 
      parameter(ncemaxvs = nf_emaxvars)
!     variable not found 
      parameter(ncenotvr = nf_enotvar)
!     action prohibited on ncglobal varid 
      parameter(nceglob = nf_eglobal)
!     not a netcdf file 
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname) 
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)

!     
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!     
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)

!
!     default fill values.  these must be the same as in the c interface.
!
      integer filbyte
      integer filchar
      integer filshort
      integer fillong
      real filfloat
      doubleprecision fildoub

      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690e+36)
      parameter (fildoub = 9.9692099683868690e+36)
   
   character(len=200), intent(in)  :: input_file       ! 1 file nane.
   character(len=10),  intent(in)  :: var              ! Variable to search for.
   integer,            intent(in)  :: field_dims       ! # Dimensions of field. 
   integer,            intent(in)  :: dim1             ! Dimension 1 of field. 
   integer,            intent(in)  :: dim2             ! Dimension 2 of field. 
   integer,            intent(in)  :: dim3             ! Dimension 3 of field. 
   integer,            intent(in)  :: k                ! Vertical index.
   real,               intent(out) :: field(1:dim1,1:dim2) ! Output field

   integer                   :: cdfid            ! 1 file id.
   integer                   :: rcode            ! Return code(0=ok).
   integer                   :: length           ! Length of filename.
   integer                   :: id_var           ! 1 variable ID. 

   integer                   :: istart(4)        ! Start value of arrays.
   integer                   :: iend(4)          ! End value of arrays.
   real(kind=4), allocatable :: field1d(:)       ! Used if 1D field read. 
   real(kind=4), allocatable :: field2d(:,:)     ! Used if 2D field read. 
   real(kind=4), allocatable :: field3d(:,:,:)   ! Used if 3D field read. 

   if (trace_use_dull) call da_trace_entry("da_get_field")

   length = len_trim(input_file)
   rcode = nf_open( input_file(1:length), NF_NOwrite, cdfid)
   if (rcode /= 0) then
      write(message(1),'(3a,i0)')' nf_open(',input_file(1:length),') returned ',rcode
      call da_error("da_get_field.inc",37,message(1:1))
   end if

   !  Check variable is in file:
   rcode = nf_inq_varid( cdfid, var, id_var)
   if (rcode /= 0) then
      write(message(1),'(3a)')var,' variable is not in input file ',input_file(1:length)
      call da_error("da_get_field.inc",44,message(1:1))
   end if

   istart = 1
   iend(1) = dim1
   iend(2) = dim2
   iend(4) = 1          ! Single time assumed.

   if (field_dims == 1) then
      iend(2) = 1
      iend(3) = 1
      allocate( field1d(1:dim1))
      call ncvgt( cdfid, id_var, istart, iend, field1d, rcode)
      field(:,1) = field1d(:)
      rcode = nf_close( cdfid)
      deallocate( field1d)
   else if (field_dims == 2) then
      iend(3) = 1
      allocate( field2d(1:dim1,1:dim2))
      call ncvgt( cdfid, id_var, istart, iend, field2d, rcode)
      field(:,:) = field2d(:,:)
      rcode = nf_close( cdfid)
      deallocate( field2d)
   else if (field_dims == 3) then
      iend(3) = dim3
      allocate( field3d(1:dim1,1:dim2,1:dim3))
      call ncvgt( cdfid, id_var, istart, iend, field3d, rcode)
      field(:,:) = field3d(:,:,k)
      deallocate( field3d)
   end if

   rcode = nf_close( cdfid)

   if (trace_use_dull) call da_trace_exit("da_get_field")

end subroutine da_get_field


subroutine da_get_height( input_file, dim1, dim2, dim3, height)

   !---------------------------------------------------------------------------
   ! Purpose: Calculates T, RH from input WRF file.
   !---------------------------------------------------------------------------
   
   implicit none
   
   character(len=200), intent(in)  :: input_file       ! 1 file nane.
   integer,            intent(in)  :: dim1, dim2, dim3          ! Dimensions.
   real,               intent(out) :: height(1:dim1,1:dim2,1:dim3) ! Height.

   character(len=10) :: var                            ! Variable to search for.
   integer           :: k                              ! Loop counter.
   real              :: gravity_inv                    ! 1/gravity.
   real              :: phb(1:dim1,1:dim2)             ! Base state geopotential.
   real              :: ph(1:dim1,1:dim2)              ! Perturbation geopotential.
   real              :: phfull(1:dim1,1:dim2,1:dim3+1) ! Geopotential.

   if (trace_use) call da_trace_entry("da_get_height")

   gravity_inv = 1.0 / gravity

   do k = 1, dim3+1
      var = "PHB"
      call da_get_field( input_file, var, 3, dim1, dim2, dim3+1, k, phb)
      var = "PH"
      call da_get_field( input_file, var, 3, dim1, dim2, dim3+1, k, ph)

      phfull(:,:,k) = phb + ph ! Calculate full geopotential on full(w) model levels:
   end do

   do k = 1, dim3
      height(:,:,k) = 0.5 *( phfull(:,:,k+1) + phfull(:,:,k)) * gravity_inv
   end do

   if (trace_use) call da_trace_exit("da_get_height")

end subroutine da_get_height


subroutine da_get_trh( input_file, dim1, dim2, dim3, k, temp, rh )

   !---------------------------------------------------------------------------
   ! Purpose: Calculates T, RH from input WRF file.
   !---------------------------------------------------------------------------

   implicit none

   character(len=200), intent(in)  :: input_file       ! 1 file name.
   integer,            intent(in)  :: dim1, dim2, dim3          ! Dimensions.
   integer,            intent(in)  :: k                         ! Model level.
   real,               intent(out) :: temp(1:dim1,1:dim2)       ! Temperature.
   real,               intent(out) :: rh(1:dim1,1:dim2)         ! Relative humidity.

   character(len=10) :: var                       ! Variable to search for. var = "T"
   integer           :: i, j                      ! Loop counters.

   real              :: thetap(1:dim1,1:dim2)     ! Perturbation potential temperature.
   real              :: pb(1:dim1,1:dim2)         ! Base state pressure.
   real              :: pp(1:dim1,1:dim2)         ! Pressure perturbation.
   real              :: x(1:dim1,1:dim2)          ! Vapor mixing ratio.

   real              :: theta                     ! Potential temperature.
   real              :: p                         ! Pressure.
   real              :: q                         ! Specific humidity.
   real              :: t_c                       ! Temp(Celsius).
   real              :: es                        ! Saturation vapor pressure.
   real              :: qs                        ! Saturation specific humidity.

   if (trace_use) call da_trace_entry("da_get_trh")

   var = "T" ! Perturbation potential temperature in WRF.
   call da_get_field( input_file, var, 3, dim1, dim2, dim3, k, thetap)

   var = "PB"  ! Base state pressure in WRF.
   call da_get_field( input_file, var, 3, dim1, dim2, dim3, k, pb)

   var = "P" ! Perturbation pressure in WRF.
   call da_get_field( input_file, var, 3, dim1, dim2, dim3, k, pp)

   var = "QVAPOR"  ! Water vapor mixing ratio.
   call da_get_field( input_file, var, 3, dim1, dim2, dim3, k, x)

   do j = 1, dim2
      do i = 1, dim1
         ! Convert p', theta' to T:
         theta = t0 + thetap(i,j)                 ! Theta = Theta0 + Thetap
         p = pb(i,j) + pp(i,j)                     ! p = base p + pert p.
         temp(i,j) = theta *( p/base_pres)**kappa ! Theta to T.

         ! Convert to specific humidity.
         q = x(i,j) /( 1.0 + x(i,j))

         ! Calculate relative humidity:
         t_c = temp(i,j) - t_kelvin
         es = es_alpha * exp( es_beta * t_c /( t_c + es_gamma))
         qs = rd_over_rv * es /( p - rd_over_rv1 * es)
         rh(i,j) = 100.0 * q / qs
      end do
   end do

   if (trace_use) call da_trace_exit("da_get_trh")

end subroutine da_get_trh


subroutine da_print_be_stats_h_global(outunit, variable, k, max_wavenumber, &
  total_power)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none 

   integer,      intent(inout) :: outunit        ! Output file unit.
   character*10, intent(in)    :: variable       ! Variable name.
   integer,      intent(in)    :: k              ! Vertical index.
   integer,      intent(in)    :: max_wavenumber ! Smallest scale required (ni/2-1)

   ! Total Power spectrum (averaged over time/members)
   real,         intent(in) :: total_power(0:max_wavenumber)

   integer :: n                     ! Loop counter.
   real    :: accum_power, variance ! Accumulated power, variance.

   if (trace_use) call da_trace_entry("da_print_be_stats_h_global")

   accum_power = 0.0
   variance = sum(total_power(0:max_wavenumber))

   write(unit=stdout,fmt='(3a,i5,a,i5)')' Power spectra for variable ', trim(variable), &
                          ' and level ', k, ' in unit ', outunit

   do n = 0, max_wavenumber
      accum_power = accum_power + total_power(n)
      write(unit=outunit,fmt='(2i4,2f18.5,f8.3)')k, n, total_power(n), accum_power, &
                                        accum_power / variance
   end do

   outunit = outunit + 1
   write(unit=stdout,fmt=*) ''

   if (trace_use) call da_trace_exit("da_print_be_stats_h_global")

end subroutine da_print_be_stats_h_global


subroutine da_print_be_stats_h_regional(outunit, variable, nk, scale_length)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none 

   integer,      intent(inout) :: outunit            ! Output file unit.
   character*10, intent(in)    :: variable           ! Variable name.
   integer,      intent(in)    :: nk                 ! Dimension of vertical index
   real,         intent(in)    :: scale_length(1:nk) ! Correlation scale lengths

   integer :: k                  ! Loop counter.

   if (trace_use) call da_trace_entry("da_print_be_stats_h_regional")

   write(unit=stdout,fmt='(3a,i5)') &
      ' Scale length for variable ', trim(variable), ' in unit ', outunit

   do k = 1, nk
     write(unit=outunit,fmt='(i4,1pe15.5)')k, scale_length(k)
   end do

   outunit = outunit + 1
   write(unit=stdout,fmt=*) ' '

   if (trace_use) call da_trace_exit("da_print_be_stats_h_regional")

end subroutine da_print_be_stats_h_regional


subroutine da_print_be_stats_p(outunit, ni, nj, nk, num_bins, num_bins2d, bin, &
   bin2d, regcoeff1, regcoeff2, regcoeff3)

   !----------------------------------------------------------------------------
   ! Purpose: To print out the regression coefficients for the physical
   !           transform.
   !
   ! Input   : outunit              --- Fortran unit for writing out
   !           ni,nj,nk             --- Dimensions
   !           num_bins, num_bins2d --- Number of the 3d and 2d bins
   !           bin, bin2d           --- bin index for the gridpoints
   !           regcoeff1   --- Reg coefs for chi = regcoeff1 * psi
   !           regcoeff2   --- Reg coefs for Ps  = sum(regcoeff2(1:k)*psi(1:k))
   !           regcoeff3   --- Reg coefs for T(k)= sum(regcoeff3(k,1:k)*psi(1:k))
   !
   ! Output  : fort.171    --- regcoef1
   !           fort.172    --- regcoeff2
   !           fort.173    --- regcoeff3
   !----------------------------------------------------------------------------

   implicit none

   integer, intent(inout) :: outunit                      ! Output file unit.
   integer, intent(in) :: ni, nj, nk                      ! Grid dimensions. 
   integer, intent(in) :: num_bins                        ! Number of 3D field bins.
   integer, intent(in) :: num_bins2d                      ! Number of 2D field bins.
   integer, intent(in) :: bin(1:ni,1:nj,1:nk)             ! Bin assigned to each 3D point.
   integer, intent(in) :: bin2d(1:ni,1:nj)                ! Bin assigned to each 3D point.
   real, intent(in)    :: regcoeff1(1:num_bins)           ! psi/chi regression cooefficient.
   real, intent(in)    :: regcoeff2(1:nk,1:num_bins2d)    ! psi/ps regression cooefficient.
   real, intent(in)    :: regcoeff3(1:nk,1:nk,1:num_bins2d) ! psi/T regression cooefficient.

   integer             :: k1, k2, i, k, j, b, number      ! Loop counters.

   if (trace_use) call da_trace_entry("da_print_be_stats_p")

   write(unit=stdout,fmt='(a,i5)')' psi/chi regression coefficient in unit ', outunit

   open(unit=outunit)
   write(unit=outunit,fmt='(a,i5)')' psi/chi regression coefficient in unit ', outunit
   do b = 1, num_bins
      number = 0
      do i = 1,ni
         do j = 1,nj
            do k = 1,nk
               if (bin(i,j,k) == b) then
                  number = number + 1
               end if
            end do
         end do
      end do
      write(unit=outunit,fmt='("bin=",i6," number_in_bin=",i6,5X,1pE12.5)') &
         b, number, regcoeff1(b)
   end do
   close(unit=outunit)

   outunit = outunit + 1

   write(unit=stdout,fmt='(a,i5)') &
      ' psi/ps  regression coefficient in unit ', outunit

   open(unit=outunit)
   write(unit=outunit,fmt='(a,i5)') &
      ' psi/ps  regression coefficient in unit ', outunit

   do b = 1, num_bins2d
      number = 0
      do i = 1,ni
         do j = 1,nj
            if (bin2d(i,j) == b) then
               number = number + 1
            end if
         end do
      end do
      write(unit=outunit,fmt='(/"bin=",i6," number_in_bin=",i6)') b, number
      write(unit=outunit,fmt='((2X,9(i3,2x,1pE12.5)))') (k,regcoeff2(k,b), k=1,nk)
   end do
   close(unit=outunit)

   outunit = outunit + 1

   write(unit=stdout,fmt='(a,i5)') &
      ' psi/T   regression coefficient in unit ', outunit

   open(unit=outunit)
   write(unit=outunit,fmt='(a,i5)') &
      ' psi/T   regression coefficient in unit ', outunit
   do b = 1, num_bins2d
      number = 0
      do i = 1,ni
         do j = 1,nj
            if (bin2d(i,j) == b) then
               number = number + 1
            end if
         end do
      end do

      write(unit=outunit,fmt='(/"bin=",i6," number_in_bin=",i6)') b, number
      do k1 = 1,nk
         write(unit=outunit,fmt='(2X,"Temperature at k1=",i3)') k1
         write(unit=outunit,fmt='((2X,9(i3,2x,1pE12.5)))') (k2,regcoeff3(k1,k2,b), k2=1,nk)
      end do
   end do
   close(unit=outunit)

   write(unit=stdout,fmt=*) ' '

   outunit = outunit + 1

   if (trace_use) call da_trace_exit("da_print_be_stats_p")

end subroutine da_print_be_stats_p


subroutine da_print_be_stats_v(outunit, variable, nk, num_bins2d, &
   e_vec, e_val, e_vec_loc, e_val_loc)

   !----------------------------------------------------------------------------
   ! Purpose: To print out the first global and local eigenvector and 
   !           eigenvalues for 'variable'.
   !
   ! Input   : outunit    --- output unit
   !           variable   --- variable name: psi, chi_u, t_u, rh
   !           nk         --- the number of vertical modes or levels
   !           num_bins2d --- the number of the 2d bins
   !           e_vec, e_val         --- global eigenvector and eigenvalues
   !           e_vec_loc, e_val_loc --- local eigenvectors and eigenvalues
   !
   ! Output  : fort.174,178,182,186,(190) --- The first 5 global eigenvectors
   !           fort.175,179,183,187,(191) --- The global eigenvalues
   !           fort.176,180,184,188,(192) --- The first 5 local eigenvectors
   !           fort.177,181,185,189,(193) --- The first 5 local eigenvalues
   !      
   !          * in parenthisis, the units for 2d fields ps_u (ps).
   !----------------------------------------------------------------------------

   implicit none

   integer, intent(inout)   :: outunit                    ! Output file unit.
   character*10, intent(in) :: variable                   ! Variable name
   integer, intent(in)      :: nk                         ! Vertical dimension.
   integer, intent(in)      :: num_bins2d                 ! # bins for 2D fields.
   real, intent(in)         :: e_vec(1:nk,1:nk)           ! Domain-averaged eigenvectors.
   real, intent(in)         :: e_val(1:nk)                ! Domain-averaged eigenvalues.
   real, intent(in) :: e_vec_loc(1:nk,1:nk,1:num_bins2d)  ! Latitudinally varying eigenvectors.
   real, intent(in) :: e_val_loc(1:nk,1:num_bins2d)       ! Latitudinally varying eigenvalues.

   integer                  :: k, m, b, mn                ! Loop counters.

   if (trace_use) call da_trace_entry("da_print_be_stats_v")

   if (nk > 5) then
      mn = 5
   else
      mn = nk
   end if

   ! 1, Global vectors:
   write(unit=stdout,fmt='(3a,i5)')' First 5 Global eigenvectors for variable ', trim(variable), &
                     ' in unit ', outunit

   open(unit=outunit)
   do k = 1, nk
      write(unit=outunit,fmt='(i4,5f15.8)') k, (e_vec(k,m), m = 1, mn)
   end do
   close(unit=outunit)

   ! 2, Global values:
   outunit = outunit + 1

   write(unit=stdout,fmt='(3a,i5)')' Global eigenvalues for variable ', trim(variable), &
                     ' in unit ', outunit

   open(unit=outunit)
   do k = 1, nk
      write(unit=outunit,fmt='(i4,1pe18.5)') k, e_val(k)
   end do
   close(unit=outunit)

   ! 3, Local vectors:

   outunit = outunit + 1

   write(unit=stdout,fmt='(3a,i5)')' First 5 local eigenvectors for variable ', trim(variable), &
                     ' in unit ', outunit

   open(unit=outunit)
   do b = 1, num_bins2d
     write(unit=outunit,fmt='(/"bin =",i6)') b
     do k = 1, nk
       write(unit=outunit,fmt='(i4,5f15.8)') k, (e_vec_loc(k,m,b), m = 1, mn)
     end do
   end do
   close(unit=outunit)

   ! 4. Local values:

   outunit = outunit + 1

   write(unit=stdout,fmt='(3a,i5)')' First 5 local eigenvalues for variable ', trim(variable), &
                     ' in unit ', outunit

   open(unit=outunit)
   do b = 1, num_bins2d
      write(unit=outunit,fmt='(i4,5(1pe18.5))') b, (e_val_loc(m,b), m = 1, mn)
   end do
   close(unit=outunit)

   outunit = outunit + 1 
   write(unit=stdout,fmt=*) ' '

   if (trace_use) call da_trace_exit("da_print_be_stats_v")

end subroutine da_print_be_stats_v


subroutine da_readwrite_be_stage1(outunit, nk)

   ! ----------------------------------------------------------------------
   ! Purpose: Read and write the dimensions and bin information from stage 1
   !
   !  Note: Please acknowledge author/institute in work that uses this code.
   !
   ! ----------------------------------------------------------------------

   implicit none

   integer, intent(in)      :: outunit                    ! Output unit number.
   integer, intent(out)     :: nk                         ! Number of vertical levels/modes.
   character*10        :: start_date, end_date       ! Starting and ending dates.
   character*10        :: date                       ! Current date (ccyymmddhh).
   character*10        :: be_method                  ! BE Method (NMC or ENS)
   integer             :: interval                   ! Period between dates (hours).
   integer             :: ne                         ! Number of ensemble members.
 
   character(len=filename_len) :: filename                   ! Input filename.
   character(len=filename_len) :: dat_dir                    ! Input data directory.

   integer                  :: ni, nj                     ! Number of points in x/y direction.
   integer                  :: bin_type                   ! Type of bin to average over. !!!DALE ADD.
   integer                  :: num_bins                   ! Number of 3D bins.
   integer                  :: num_bins2d                 ! Number of 2D bins.

   real                :: hgt_min, hgt_max           ! Used if bin_type = 2 (m).
   real                :: lat_min, lat_max           ! Used if bin_type = 2 (degrees).
   real                :: binwidth_lat               ! Used if bin_type = 2 (degrees). !!!DALE ADD..
   real                :: binwidth_hgt               ! Used if bin_type = 2 (m). !!!DALE ADD..
   integer, allocatable:: bin(:,:,:)                 ! Bin assigned to each 3D point.
   integer, allocatable:: bin2d(:,:)                 ! Bin assigned to each 2D point.


   integer :: iunit, namelist_unit

   namelist / gen_be_stage1_nl / start_date, end_date, interval, &
                                 be_method, ne, bin_type, cv_options, &
                                 lat_min, lat_max, binwidth_lat, &
                                 hgt_min, hgt_max, binwidth_hgt, dat_dir

   start_date = '2004030312'
   end_date = '2004033112'
   interval = 24
   be_method = 'NMC'
   ne = 1
   bin_type = 5         ! 0 = Every pt, 1 = x direction, 2 = latitude, ....
   cv_options = 5
   lat_min = -90.0
   lat_max = 90.0
   binwidth_lat = 10.0
   hgt_min = 0.0
   hgt_max = 20000.0
   binwidth_hgt = 1000.0
   dat_dir = 'NO_DIRECTORY_SPECIFIED'

   call da_get_unit(namelist_unit)

   ! Reading Namelist:            
   open(unit=namelist_unit, file='gen_be_stage1_nl.nl', &
        form='formatted', status='old', action='read')
   read(namelist_unit, gen_be_stage1_nl)
   close(namelist_unit)
   call da_free_unit(namelist_unit)
  
   ! Read domain info from stage 1 since we skipped stage 2
   call da_get_unit(iunit)
   filename = trim(dat_dir)//'/pert.'//start_date(1:10)//'.e001'
   open (iunit, file = trim(filename), form='unformatted')
   read(iunit)date, ni, nj, nk
   close(iunit)

   ! Read in the bin information:

   allocate(bin(1:ni,1:nj,1:nk))
   allocate(bin2d(1:ni,1:nj))
 
   filename = 'bin.data'
   open (iunit, file = filename, form='unformatted')

   read(iunit)bin_type
   read(iunit)lat_min, lat_max, binwidth_lat
   read(iunit)hgt_min, hgt_max, binwidth_hgt
   read(iunit)num_bins, num_bins2d
   read(iunit)bin(1:ni,1:nj,1:nk)
   read(iunit)bin2d(1:ni,1:nj)
   close(iunit)
   call da_free_unit(iunit)

   ! Write out the dimensions and bin information:

   write(outunit)ni, nj, nk

   write(outunit)bin_type
   write(outunit)lat_min, lat_max, binwidth_lat
   write(outunit)hgt_min, hgt_max, binwidth_hgt
   write(outunit)num_bins, num_bins2d
   write(outunit)bin(1:ni,1:nj,1:nk)
   write(outunit)bin2d(1:ni,1:nj)

   deallocate(bin)
   deallocate(bin2d)

end subroutine da_readwrite_be_stage1


subroutine da_readwrite_be_stage2(outunit, nk)

   ! ----------------------------------------------------------------------
   ! Purpose: Read and write the dimensions, bin information, and 
   !          regression coefficients.
   !  Update: Multivariate BE option (cv_options=6)
   !          Syed RH Rizvi (MMM/NESL/NCAR)   Date: 02/01/2010
   !
   !  Note: Please acknowledge author/institute in work that uses this code.
   !
   ! ----------------------------------------------------------------------

   implicit none

   integer, intent(in)      :: outunit                    ! Output unit number.
   integer, intent(out)     :: nk                         ! Number of vertical levels/modes.
   character*10        :: start_date, end_date       ! Starting and ending dates.
   character*10        :: date, new_date             ! Current date (ccyymmddhh).
   integer             :: interval                   ! Period between dates (hours).
   integer             :: ne                         ! Number of ensemble members.
   logical             :: testing_eofs               ! True if testing EOF decomposition.
   real                :: rf_scale                   ! Recursive filter scale.
   integer             :: num_passes                 ! Recursive filter passes.
 
   character(len=filename_len) :: filename                   ! Input filename.
   character*80             :: variable_psi_chi
   character*80             :: variable_psi_t
   character*80             :: variable_psi_ps
   character*80             :: variable_psi_rh
   character*80             :: variable_chi_u_t
   character*80             :: variable_chi_u_ps
   character*80             :: variable_chi_u_rh
   character*80             :: variable_t_u_rh
   character*80             :: variable_ps_u_rh

   integer                  :: ni, nj                     ! Number of points in x/y direction.
   integer                  :: bin_type                   ! Type of bin to average over. !!!DALE ADD.
   integer                  :: num_bins                   ! Number of 3D bins.
   integer                  :: num_bins2d                 ! Number of 2D bins.

   real                :: hgt_min, hgt_max           ! Used if bin_type = 2 (m).
   real                :: lat_min, lat_max           ! Used if bin_type = 2 (degrees).
   real                :: binwidth_lat               ! Used if bin_type = 2 (degrees). !!!DALE ADD..
   real                :: binwidth_hgt               ! Used if bin_type = 2 (m). !!!DALE ADD..
   integer, allocatable:: bin(:,:,:)                 ! Bin assigned to each 3D point.
   integer, allocatable:: bin2d(:,:)                 ! Bin assigned to each 2D point.


   real, allocatable   :: regcoeff_psi_chi(:)        ! chi/psi regression cooefficient.
   real, allocatable   :: regcoeff_psi_ps(:,:)       ! ps/psi regression cooefficient.
   real, allocatable   :: regcoeff_psi_t(:,:,:)      ! t/psi regression cooefficient.
   real, allocatable   :: regcoeff_psi_rh(:,:,:)     ! rh/psi regression cooefficient.

   real, allocatable   :: regcoeff_chi_u_ps(:,:)     ! ps/chi_u regression coefficient
   real, allocatable   :: regcoeff_chi_u_t(:,:,:)    ! t/chi_u regression coefficient
   real, allocatable   :: regcoeff_ps_u_rh(:,:)      ! rh/ps_u regression coefficient
   real, allocatable   :: regcoeff_chi_u_rh(:,:,:)   ! rh/chi_u regression coefficient
   real, allocatable   :: regcoeff_t_u_rh(:,:,:)     ! rh/t_u regression coefficient

   integer :: iunit, namelist_unit

   namelist / gen_be_stage2_nl / start_date, end_date, interval, &
                                 ne, testing_eofs, num_passes, rf_scale, cv_options

   start_date = '2004030312'
   end_date = '2004033112'
   interval = 24
   ne = 1
   testing_eofs = .true.
   num_passes = 0
   rf_scale = 1.0
   cv_options = 5

   call da_get_unit(namelist_unit)

   ! Reading Namelist:            
   open(unit=namelist_unit, file='gen_be_stage2_nl.nl', &
        form='formatted', status='old', action='read')
   read(namelist_unit, gen_be_stage2_nl)
   close(namelist_unit)
   call da_free_unit(namelist_unit)
  
   ! Read in the coefficients:

   call da_get_unit(iunit)
   filename = 'gen_be_stage2.dat'
   open (iunit, file = filename, form='unformatted')
   read(iunit)ni, nj, nk
   read(iunit)num_bins, num_bins2d

   allocate( regcoeff_psi_chi(1:num_bins) )
   allocate( regcoeff_psi_ps(1:nk,1:num_bins2d) )
   allocate( regcoeff_psi_t(1:nk,1:nk,1:num_bins2d) )
   allocate( regcoeff_psi_rh(1:nk,1:nk,1:num_bins2d) )

   allocate( regcoeff_chi_u_ps(1:nk,1:num_bins2d) )
   allocate( regcoeff_chi_u_t(1:nk,1:nk,1:num_bins2d) )
   allocate( regcoeff_chi_u_rh(1:nk,1:nk,1:num_bins2d) )

   allocate( regcoeff_ps_u_rh(1:nk,1:num_bins2d) )
   allocate( regcoeff_t_u_rh(1:nk,1:nk,1:num_bins2d) )

   if( cv_options == 5) then
   read(iunit)regcoeff_psi_chi
   read(iunit)regcoeff_psi_ps
   read(iunit)regcoeff_psi_t
   end if

   if( cv_options == 7) then
   read(iunit)regcoeff_psi_chi
   read(iunit)regcoeff_psi_ps
   read(iunit)regcoeff_psi_t
   end if

   if( cv_options == 6) then
   read(iunit)variable_psi_chi
   read(iunit)regcoeff_psi_chi

   read(iunit)variable_psi_t
   read(iunit)regcoeff_psi_t

   read(iunit)variable_psi_ps
   read(iunit)regcoeff_psi_ps

   read(iunit)variable_psi_rh
   read(iunit)regcoeff_psi_rh

   read(iunit)variable_chi_u_t
   read(iunit)regcoeff_chi_u_t

   read(iunit)variable_chi_u_ps
   read(iunit)regcoeff_chi_u_ps

   read(iunit)variable_chi_u_rh
   read(iunit)regcoeff_chi_u_rh

   read(iunit)variable_t_u_rh
   read(iunit)regcoeff_t_u_rh

   read(iunit)variable_ps_u_rh
   read(iunit)regcoeff_ps_u_rh

   end if
   close(iunit)

   ! Read in the bin information:

   allocate(bin(1:ni,1:nj,1:nk))
   allocate(bin2d(1:ni,1:nj))
 
   filename = 'bin.data'
   open (iunit, file = filename, form='unformatted')

   read(iunit)bin_type
   read(iunit)lat_min, lat_max, binwidth_lat
   read(iunit)hgt_min, hgt_max, binwidth_hgt
   read(iunit)num_bins, num_bins2d
   read(iunit)bin(1:ni,1:nj,1:nk)
   read(iunit)bin2d(1:ni,1:nj)
   close(iunit)
   call da_free_unit(iunit)

   ! Write out the dimensions and bin information:

   write(outunit)ni, nj, nk

   write(outunit)bin_type
   write(outunit)lat_min, lat_max, binwidth_lat
   write(outunit)hgt_min, hgt_max, binwidth_hgt
   write(outunit)num_bins, num_bins2d
   write(outunit)bin(1:ni,1:nj,1:nk)
   write(outunit)bin2d(1:ni,1:nj)

   deallocate(bin)
   deallocate(bin2d)

   ! Write out the coefficients:

   if( cv_options == 5) then
   write(outunit)regcoeff_psi_chi
   write(outunit)regcoeff_psi_ps
   write(outunit)regcoeff_psi_t
   end if

   if( cv_options == 7) then
   write(outunit)regcoeff_psi_chi
   write(outunit)regcoeff_psi_ps
   write(outunit)regcoeff_psi_t
   end if

   if( cv_options == 6) then
   write(outunit)variable_psi_chi
   write(outunit)regcoeff_psi_chi

   write(outunit)variable_psi_t
   write(outunit)regcoeff_psi_t

   write(outunit)variable_psi_ps
   write(outunit)regcoeff_psi_ps

   write(outunit)variable_psi_rh
   write(outunit)regcoeff_psi_rh

   write(outunit)variable_chi_u_t
   write(outunit)regcoeff_chi_u_t

   write(outunit)variable_chi_u_ps
   write(outunit)regcoeff_chi_u_ps

   write(outunit)variable_chi_u_rh
   write(outunit)regcoeff_chi_u_rh

   write(outunit)variable_t_u_rh
   write(outunit)regcoeff_t_u_rh

   write(outunit)variable_ps_u_rh
   write(outunit)regcoeff_ps_u_rh
   end if

   deallocate(regcoeff_psi_chi)
   deallocate(regcoeff_psi_ps)
   deallocate(regcoeff_psi_t)
   deallocate(regcoeff_psi_rh)

   deallocate( regcoeff_chi_u_ps)
   deallocate( regcoeff_chi_u_t)
   deallocate( regcoeff_chi_u_rh)

   deallocate( regcoeff_ps_u_rh)
   deallocate( regcoeff_t_u_rh)


end subroutine da_readwrite_be_stage2


subroutine da_readwrite_be_stage3(outunit, nk, variable)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)      :: outunit                    ! Output unit number.
   integer, intent(in)      :: nk                         ! Number of vertical levels/modes.
   character*10, intent(in) :: variable                   ! Variable name

   character(len=filename_len) :: filename                   ! Input filename.
   character*10             :: vardum                     ! Dummy variable name.

   integer                  :: num_bins2d                 ! Number of Eigenvector bins.
   integer                  :: nkdum                      ! Dummy nk variable.
   real                     :: e_vec(1:nk,1:nk)           ! Domain-averaged eigenvectors.
   real                     :: e_val(1:nk)                ! Domain-averaged eigenvalues.  
   real, allocatable        :: e_vec_loc(:,:,:)           ! Latitudinally varying eigenvectors.
   real, allocatable        :: e_val_loc(:,:)             ! Latitudinally varying eigenvalues.

   integer :: iunit
 
   call da_get_unit(iunit)
   filename = 'gen_be_stage3.'//trim(variable)//'.dat'
   open (iunit, file = filename, form='unformatted')
   read(iunit)vardum
   if (trim(vardum) /= trim(variable)) then
     call da_error("da_readwrite_be_stage3.inc",30, &
       (/"Inconsistent variable name"/))
   end if

   read(iunit)nkdum, num_bins2d
   if (nkdum /= nk) then
     call da_error("da_readwrite_be_stage3.inc",36, &
       (/"Inconsistent nk between regression coefficients and vertical modes"/))
   end if

   allocate(e_vec_loc(1:nk,1:nk,1:num_bins2d))
   allocate(e_val_loc(1:nk,1:num_bins2d))

   read(iunit)e_vec
   read(iunit)e_val
   read(iunit)e_vec_loc
   read(iunit)e_val_loc
   close(iunit)
   call da_free_unit(iunit)

   write(outunit)variable
   write(outunit)nk, num_bins2d
   write(outunit)e_vec
   write(outunit)e_val
   write(outunit)e_vec_loc
   write(outunit)e_val_loc

   deallocate(e_vec_loc)
   deallocate(e_val_loc)

end subroutine da_readwrite_be_stage3


subroutine da_readwrite_be_stage4(outunit, nk, uh_method, n_smth_sl, variable)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, parameter          :: spike_tolerance = 1.5      ! Threshold for detecting spikes in data. 

   integer, intent(in)      :: outunit                    ! Output unit number.
   integer, intent(in)      :: nk                         ! Number of vertical levels/modes.
   integer    , intent(in)  :: n_smth_sl                  ! Number of smoothing. 0: no smoothig
   character*7, intent(in)  :: uh_method                  ! Uh method (power, scale, wavelet)
   character*10, intent(in) :: variable                   ! Variable name.
   character, save          :: namw1                      ! Test namw.
   character(len=filename_len) :: filename                ! Input filename.
   character*10, save       :: cvar                       ! Dummy variable name.
   character*2              :: ck                         ! Loop index -> character.
   integer                  :: iunit
   integer                  :: k                          ! Loop counter.
   integer                  :: kdum                       ! Dummy vertical index.
   integer, save            :: lf1,nb1,nij1(0:1)          ! Test lf, nb & nij.
   integer                  :: max_wavenumber             ! Smallest scale required (ni/2 - 1).
   integer                  :: stride                     ! Only for filename.
   logical                  :: ex                         ! Test if input file exists.
   logical, save            :: do_normalize1              ! Test do_normalize.
   real                     :: ml                         ! Gradient (read from file but not used).
   real                     :: mean_scale                 ! Average scale between points.
   logical, save            :: first_time=.true.          ! True if first time called.
   real                     :: scale_length(nk)
   real                     :: sl_smth(nk)
   real, allocatable        :: sd(:,:)                    ! Horizontal standard deviations.
   real, allocatable        :: total_power(:)             ! Total Power spectrum.

   call da_get_unit(iunit)

   if (trim(uh_method) .eq. 'power') then
      do k = 1, nk
         write(ck,'(i2)')k
         if (k < 10) ck = '0'//ck(2:2)

         filename = trim(variable)//'/'//trim(variable) &
                                  //'.'//ck//'.spectrum'

         open (iunit, file = filename, form='unformatted')
         read(iunit)cvar
         if (trim(cvar) /=  trim(variable)) then
           call da_error("da_readwrite_be_stage4.inc",49, &
             (/"Variable name inconsistency"/))
         end if

         read(iunit)max_wavenumber, kdum
         if (kdum /= k) then
           call da_error("da_readwrite_be_stage4.inc",55, &
             (/"Inconsistent vertical label"/))
         end if

         allocate(total_power(0:max_wavenumber))

         read(iunit)total_power
         close(iunit)

         write(outunit)variable
         write(outunit)max_wavenumber, k
         write(outunit) .false. ! preserve file format
         write(outunit)total_power

         deallocate(total_power)

      end do

   else if (trim(uh_method) == 'scale') then

      if( .not.use_rf ) then
         call da_error("da_readwrite_be_stage4.inc",76,(/"{uh_method,use_rf}={''scale'',.false.} error"/))
      end if

      filename = trim(variable)//'/'//'sl_print.'//trim(variable)
      open (iunit, file = filename)

      do k=1, nk
         read(iunit,'(a,2e20.8)') ck, ml, scale_length(k)
         ! If missing value is encountered, use the value from the last
         ! mode or level (YRG, 06/30/2005):
         if (ml == 0.0 .and. scale_length(k) == 0.0) &
             scale_length(k) = scale_length(k-1)
      end do

      ! Remove spikes in lengthscales (extrapolate if spike detected):
      do k = 2, nk-1
         mean_scale = 0.5 * ( scale_length(k-1) + scale_length(k+1) )
         if ( scale_length(k) > spike_tolerance * mean_scale ) then
            scale_length(k) = mean_scale
         end if
      end do

      ! Smoothing the scale_length
      sl_smth =  scale_length
      do kdum = 1, n_smth_sl
         do k = 2, nk-1
            sl_smth(k) = scale_length(k) &
               + 0.25*(scale_length(k-1)+scale_length(k+1)-2.0*scale_length(k))
         end do
         scale_length = sl_smth 
      end do
     
      write(outunit) variable
      write(outunit) scale_length
      close(iunit)

      if( do_normalize )then
         allocate(nij(0:0,0:1,0:0))
         do k=1, nk
            write(ck,'(i0)')k
            ex=.false.
            stride=1
            do while( .not.ex )		! loop thus only because gen_be_stage4_regional uses stride:
               write(filename,'("gen_be_stage4_regional.",i0)')stride
               inquire(exist=ex,file=trim(filename))
               if( .not.ex ) then
                  write(unit=message(1),fmt='(a,a)') trim(filename),' does not exist.'
                  call da_warning("da_readwrite_be_stage4.inc",123,message(1:1))
               end if
               stride=stride+1
            enddo
            filename = trim(filename)//"/dir."//trim(variable)//trim(ck)//"/mom"
            open(iunit,action="read",file=trim(filename),form="unformatted",status="old")
            read(iunit)do_normalize1
            if( do_normalize.neqv.do_normalize1 ) &
               call da_error("da_readwrite_be_stage4.inc",131,(/"do_normalize.neqv.do_normalize1 for uh_method==''scale''"/))
            read(iunit)nij
            if( k==1 )then
               write(outunit)do_normalize
               write(outunit)nij
               allocate(sd(nij(0,1,0),nij(0,0,0)))
            endif
            read(iunit)sd
            write(outunit)sd
            close(iunit)
         enddo
         deallocate(nij,sd)
      endif

   elseif( trim(uh_method) == 'wavelet' )then
      if( use_rf )call da_error("da_readwrite_be_stage4.inc",146,(/"{uh_method,use_rf}={''wavelet'',.true.} error"/))
      if( first_time )then
         write(outunit) len_trim('wavelet')
         write(outunit) 'wavelet'
         cvar=trim(variable)
      endif
      write(outunit) len_trim(variable),nk
      write(outunit) trim(variable)
      do k=1, nk
         write(ck,'(i0)')k
         filename = "gen_be_stage4_regional/dir."//trim(variable)//trim(ck)//"/momw"
         open(iunit,action="read",file=trim(filename),form="unformatted",status="old")
         read(iunit)do_normalize1,namw,lf,nb
         if( first_time )then		! Use same basis forall {variable,k}:
            if( do_normalize.neqv.do_normalize1 ) &
               call da_error("da_readwrite_be_stage4.inc",161,(/"do_normalize.neqv.do_normalize1 in "//filename/))
            namw1=namw
            lf1=lf
            nb1=nb
            write(outunit)do_normalize,namw,lf,nb
            allocate(nij(0:nb,0:1,0:2))
         elseif( do_normalize.neqv.do_normalize1 .or. namw/=namw1 .or. lf/=lf1 .or. nb/=nb1 )then
            write(message(1), &
               '("{nb namw lf do_normalize}[",a,",k=",i0,"]/{nb namw lf do_normalize}[",a,",1]={",i0,a,i0,l2,"}/{",i0,a,i0,l2,"}.")') &
               trim(variable),k,trim(cvar),nb,namw,lf,do_normalize,nb1,namw1,lf1,do_normalize1
            call da_error("da_readwrite_be_stage4.inc",171,message(1:1))
         endif
         read(iunit)nij
         if( first_time )then
            nij1=nij(0,:,0)
            write(outunit)nij
            allocate(wsd(nij(0,1,2),nij(0,0,2)))
            if( do_normalize )allocate(sd(nij(0,1,0),nij(0,0,0)))
            first_time=.false.
         elseif( any(nij(0,:,0)/=nij1) )then
            write(message(1), &
               '("{ni nj}[",a,",k=",i0,"]/{ni nj}[",a,",1]={",i0,1x,i0,"}/{",i0,1x,i0,"}.")') &
               trim(variable),k,trim(cvar),nij(0,:,0),nij1
            call da_error("da_readwrite_be_stage4.inc",184,message(1:1))
         endif
         read(iunit)wsd
         write(outunit)wsd
         if( do_normalize )then
            read(iunit)sd
            write(outunit)sd
         endif
         close(iunit)
      end do				! do k=1,nk
   end if

   call da_free_unit(iunit)

end subroutine da_readwrite_be_stage4
subroutine da_stage0_initialize(input_file, var, dim1, dim2, dim3, ds, mp_physics)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

!     NetCDF-3.
!
! netcdf version 3 fortran interface:
!

!
! external netcdf data types:
!
      integer nf_byte
      integer nf_int1
      integer nf_char
      integer nf_short
      integer nf_int2
      integer nf_int
      integer nf_float
      integer nf_real
      integer nf_double
      integer nf_ubyte
      integer nf_ushort
      integer nf_uint
      integer nf_int64
      integer nf_uint64

      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)

!
! default fill values:
!
      integer           nf_fill_byte
      integer           nf_fill_int1
      integer           nf_fill_char
      integer           nf_fill_short
      integer           nf_fill_int2
      integer           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double

      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690e+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690d+36)

!
! mode flags for opening and creating a netcdf dataset:
!
      integer nf_nowrite
      integer nf_write
      integer nf_clobber
      integer nf_noclobber
      integer nf_fill
      integer nf_nofill
      integer nf_lock
      integer nf_share
      integer nf_64bit_offset
      integer nf_64bit_data
      integer nf_cdf5
      integer nf_sizehint_default
      integer nf_align_chunk
      integer nf_format_classic
      integer nf_format_64bit
      integer nf_format_64bit_offset
      integer nf_format_64bit_data
      integer nf_format_cdf5
      integer nf_diskless
      integer nf_mmap

      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)

!
! size argument for defining an unlimited dimension:
!
      integer nf_unlimited
      parameter (nf_unlimited = 0)

!
! global attribute id:
!
      integer nf_global
      parameter (nf_global = 0)

!
! implementation limits:
!
      integer nf_max_dims
      integer nf_max_attrs
      integer nf_max_vars
      integer nf_max_name
      integer nf_max_var_dims

      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)

!
! error codes:
!
      integer nf_noerr
      integer nf_ebadid
      integer nf_eexist
      integer nf_einval
      integer nf_eperm
      integer nf_enotindefine
      integer nf_eindefine
      integer nf_einvalcoords
      integer nf_emaxdims
      integer nf_enameinuse
      integer nf_enotatt
      integer nf_emaxatts
      integer nf_ebadtype
      integer nf_ebaddim
      integer nf_eunlimpos
      integer nf_emaxvars
      integer nf_enotvar
      integer nf_eglobal
      integer nf_enotnc
      integer nf_ests
      integer nf_emaxname
      integer nf_eunlimit
      integer nf_enorecvars
      integer nf_echar
      integer nf_eedge
      integer nf_estride
      integer nf_ebadname
      integer nf_erange
      integer nf_enomem
      integer nf_evarsize
      integer nf_edimsize
      integer nf_etrunc

      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
!
! error handling modes:
!
      integer  nf_fatal
      integer nf_verbose

      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)

!
! miscellaneous routines:
!
      character*80   nf_inq_libvers
      external       nf_inq_libvers

      character*80   nf_strerror
!                         (integer             ncerr)
      external       nf_strerror

      logical        nf_issyserr
!                         (integer             ncerr)
      external       nf_issyserr

!
! control routines:
!
      integer         nf_inq_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_inq_base_pe

      integer         nf_set_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_set_base_pe

      integer         nf_create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             ncid)
      external        nf_create

      integer         nf__create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create

      integer         nf__create_mp
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create_mp

      integer         nf_open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             ncid)
      external        nf_open

      integer         nf__open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open

      integer         nf__open_mp
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open_mp

      integer         nf_set_fill
!                         (integer             ncid,
!                          integer             fillmode,
!                          integer             old_mode)
      external        nf_set_fill

      integer         nf_set_default_format
!                          (integer             format,
!                          integer             old_format)
      external        nf_set_default_format

      integer         nf_redef
!                         (integer             ncid)
      external        nf_redef

      integer         nf_enddef
!                         (integer             ncid)
      external        nf_enddef

      integer         nf__enddef
!                         (integer             ncid,
!                          integer             h_minfree,
!                          integer             v_align,
!                          integer             v_minfree,
!                          integer             r_align)
      external        nf__enddef

      integer         nf_sync
!                         (integer             ncid)
      external        nf_sync

      integer         nf_abort
!                         (integer             ncid)
      external        nf_abort

      integer         nf_close
!                         (integer             ncid)
      external        nf_close

      integer         nf_delete
!                         (character*(*)       ncid)
      external        nf_delete

!
! general inquiry routines:
!

      integer         nf_inq
!                         (integer             ncid,
!                          integer             ndims,
!                          integer             nvars,
!                          integer             ngatts,
!                          integer             unlimdimid)
      external        nf_inq

! new inquire path

      integer nf_inq_path
      external nf_inq_path

      integer         nf_inq_ndims
!                         (integer             ncid,
!                          integer             ndims)
      external        nf_inq_ndims

      integer         nf_inq_nvars
!                         (integer             ncid,
!                          integer             nvars)
      external        nf_inq_nvars

      integer         nf_inq_natts
!                         (integer             ncid,
!                          integer             ngatts)
      external        nf_inq_natts

      integer         nf_inq_unlimdim
!                         (integer             ncid,
!                          integer             unlimdimid)
      external        nf_inq_unlimdim

      integer         nf_inq_format
!                         (integer             ncid,
!                          integer             format)
      external        nf_inq_format

!
! dimension routines:
!

      integer         nf_def_dim
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             len,
!                          integer             dimid)
      external        nf_def_dim

      integer         nf_inq_dimid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             dimid)
      external        nf_inq_dimid

      integer         nf_inq_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_dim

      integer         nf_inq_dimname
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_inq_dimname

      integer         nf_inq_dimlen
!                         (integer             ncid,
!                          integer             dimid,
!                          integer             len)
      external        nf_inq_dimlen

      integer         nf_rename_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_rename_dim

!
! general attribute routines:
!

      integer         nf_inq_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len)
      external        nf_inq_att

      integer         nf_inq_attid
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             attnum)
      external        nf_inq_attid

      integer         nf_inq_atttype
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype)
      external        nf_inq_atttype

      integer         nf_inq_attlen
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_attlen

      integer         nf_inq_attname
!                         (integer             ncid,
!                          integer             varid,
!                          integer             attnum,
!                          character(*)        name)
      external        nf_inq_attname

      integer         nf_copy_att
!                         (integer             ncid_in,
!                          integer             varid_in,
!                          character(*)        name,
!                          integer             ncid_out,
!                          integer             varid_out)
      external        nf_copy_att

      integer         nf_rename_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        curname,
!                          character(*)        newname)
      external        nf_rename_att

      integer         nf_del_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_del_att

!
! attribute put/get routines:
!

      integer         nf_put_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len,
!                          character(*)        text)
      external        nf_put_att_text

      integer         nf_get_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          character(*)        text)
      external        nf_get_att_text

      integer         nf_put_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int1_t           i1vals(1))
      external        nf_put_att_int1

      integer         nf_get_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int1_t           i1vals(1))
      external        nf_get_att_int1

      integer         nf_put_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int2_t           i2vals(1))
      external        nf_put_att_int2

      integer         nf_get_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int2_t           i2vals(1))
      external        nf_get_att_int2

      integer         nf_put_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          integer             ivals(1))
      external        nf_put_att_int

      integer         nf_get_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             ivals(1))
      external        nf_get_att_int

      integer         nf_put_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int8_t           i8vals(1))
      external        nf_put_att_int64

      integer         nf_get_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int8_t           i8vals(1))
      external        nf_get_att_int64

      integer         nf_put_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          real                rvals(1))
      external        nf_put_att_real

      integer         nf_get_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          real                rvals(1))
      external        nf_get_att_real

      integer         nf_put_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          double              dvals(1))
      external        nf_put_att_double

      integer         nf_get_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          double              dvals(1))
      external        nf_get_att_double

!
! general variable routines:
!

      integer         nf_def_var
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             varid)
      external        nf_def_var

      integer         nf_inq_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             natts)
      external        nf_inq_var

      integer         nf_inq_varid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             varid)
      external        nf_inq_varid

      integer         nf_inq_varname
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_inq_varname

      integer         nf_inq_vartype
!                         (integer             ncid,
!                          integer             varid,
!                          integer             xtype)
      external        nf_inq_vartype

      integer         nf_inq_varndims
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ndims)
      external        nf_inq_varndims

      integer         nf_inq_vardimid
!                         (integer             ncid,
!                          integer             varid,
!                          integer             dimids(1))
      external        nf_inq_vardimid

      integer         nf_inq_varnatts
!                         (integer             ncid,
!                          integer             varid,
!                          integer             natts)
      external        nf_inq_varnatts

      integer         nf_rename_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_rename_var

      integer         nf_copy_var
!                         (integer             ncid_in,
!                          integer             varid,
!                          integer             ncid_out)
      external        nf_copy_var

!
! entire variable put/get routines:
!

      integer         nf_put_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_put_var_text

      integer         nf_get_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_get_var_text

      integer         nf_put_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_put_var_int1

      integer         nf_get_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_get_var_int1

      integer         nf_put_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_put_var_int2

      integer         nf_get_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_get_var_int2

      integer         nf_put_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_put_var_int

      integer         nf_get_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_get_var_int

      integer         nf_put_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_put_var_real

      integer         nf_get_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_get_var_real

      integer         nf_put_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_put_var_double

      integer         nf_get_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_get_var_double

!
! single variable put/get routines:
!

      integer         nf_put_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_put_var1_text

      integer         nf_get_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_get_var1_text

      integer         nf_put_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_put_var1_int1

      integer         nf_get_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_get_var1_int1

      integer         nf_put_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_put_var1_int2

      integer         nf_get_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_get_var1_int2

      integer         nf_put_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_put_var1_int

      integer         nf_get_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_get_var1_int

      integer         nf_put_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_put_var1_real

      integer         nf_get_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_get_var1_real

      integer         nf_put_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_put_var1_double

      integer         nf_get_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_get_var1_double

!
! variable array put/get routines:
!

      integer         nf_put_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_put_vara_text

      integer         nf_get_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_get_vara_text

      integer         nf_put_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vara_int1

      integer         nf_get_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vara_int1

      integer         nf_put_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vara_int2

      integer         nf_get_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vara_int2

      integer         nf_put_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_put_vara_int

      integer         nf_get_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_get_vara_int

      integer         nf_put_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_put_vara_real

      integer         nf_get_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_get_vara_real

      integer         nf_put_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vara_double

      integer         nf_get_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vara_double

!
! strided variable put/get routines:
!

      integer         nf_put_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_put_vars_text

      integer         nf_get_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_get_vars_text

      integer         nf_put_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vars_int1

      integer         nf_get_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vars_int1

      integer         nf_put_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vars_int2

      integer         nf_get_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vars_int2

      integer         nf_put_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_put_vars_int

      integer         nf_get_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_get_vars_int

      integer         nf_put_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_put_vars_real

      integer         nf_get_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_get_vars_real

      integer         nf_put_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vars_double

      integer         nf_get_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vars_double

!
! mapped variable put/get routines:
!

      integer         nf_put_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_put_varm_text

      integer         nf_get_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_get_varm_text

      integer         nf_put_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_varm_int1

      integer         nf_get_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_varm_int1

      integer         nf_put_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_varm_int2

      integer         nf_get_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_varm_int2

      integer         nf_put_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_put_varm_int

      integer         nf_get_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_get_varm_int

      integer         nf_put_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_put_varm_real

      integer         nf_get_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_get_varm_real

      integer         nf_put_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_put_varm_double

      integer         nf_get_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_get_varm_double

!     64-bit int functions.
      integer nf_put_var1_int64
      external nf_put_var1_int64
      integer nf_put_vara_int64
      external nf_put_vara_int64
      integer nf_put_vars_int64
      external nf_put_vars_int64
      integer nf_put_varm_int64
      external nf_put_varm_int64
      integer nf_put_var_int64
      external nf_put_var_int64
      integer nf_get_var1_int64
      external nf_get_var1_int64
      integer nf_get_vara_int64
      external nf_get_vara_int64
      integer nf_get_vars_int64
      external nf_get_vars_int64
      integer nf_get_varm_int64
      external nf_get_varm_int64
      integer nf_get_var_int64
      external nf_get_var_int64


!     NetCDF-4.
!     This is part of netCDF-4. Copyright 2006, UCAR, See COPYRIGHT
!     file for distribution information.

!     Netcdf version 4 fortran interface.

!     $Id: netcdf4.inc,v 1.28 2010/05/25 13:53:02 ed Exp $

!     New netCDF-4 types.
      integer nf_string
      integer nf_vlen
      integer nf_opaque
      integer nf_enum
      integer nf_compound

      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)

!     New netCDF-4 fill values.
      integer           nf_fill_ubyte
      integer           nf_fill_ushort
!      real              nf_fill_uint
!      real              nf_fill_int64
!      real              nf_fill_uint64
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)

!     New constants.
      integer nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)

      integer nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)

      integer nf_netcdf4
      parameter (nf_netcdf4 = 4096)

      integer nf_classic_model
      parameter (nf_classic_model = 256)

      integer nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)

      integer nf_endian_native
      parameter (nf_endian_native = 0)
      integer nf_endian_little
      parameter (nf_endian_little = 1)
      integer nf_endian_big
      parameter (nf_endian_big = 2)

!     For NF_DEF_VAR_CHUNKING
      integer nf_chunked
      parameter (nf_chunked = 0)
      integer nf_contiguous
      parameter (nf_contiguous = 1)
      integer nf_compact
      parameter (nf_compact = 2)

!     For NF_DEF_VAR_FLETCHER32
      integer nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer nf_fletcher32
      parameter (nf_fletcher32 = 1)

!     For NF_DEF_VAR_DEFLATE
      integer nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer nf_shuffle
      parameter (nf_shuffle = 1)

!     For NF_DEF_VAR_SZIP
      integer nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)

!     For parallel I/O.
      integer nf_mpiio      
      parameter (nf_mpiio = 8192)
      integer nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer nf_pnetcdf
      parameter (nf_pnetcdf = 32768)

!     For NF_VAR_PAR_ACCESS.
      integer nf_independent
      parameter (nf_independent = 0)
      integer nf_collective
      parameter (nf_collective = 1)

!     New error codes.
      integer nf_ehdferr        ! Error at 1 layer. 
      parameter (nf_ehdferr = -101)
      integer nf_ecantread      ! Can't read. 
      parameter (nf_ecantread = -102)
      integer nf_ecantwrite     ! Can't write. 
      parameter (nf_ecantwrite = -103)
      integer nf_ecantcreate    ! Can't create. 
      parameter (nf_ecantcreate = -104)
      integer nf_efilemeta      ! Problem with file metadata. 
      parameter (nf_efilemeta = -105)
      integer nf_edimmeta       ! Problem with dimension metadata. 
      parameter (nf_edimmeta = -106)
      integer nf_eattmeta       ! Problem with attribute metadata. 
      parameter (nf_eattmeta = -107)
      integer nf_evarmeta       ! Problem with variable metadata. 
      parameter (nf_evarmeta = -108)
      integer nf_enocompound    ! Not a compound type. 
      parameter (nf_enocompound = -109)
      integer nf_eattexists     ! Attribute already exists. 
      parameter (nf_eattexists = -110)
      integer nf_enotnc4        ! Attempting netcdf-4 operation on netcdf-3 file.   
      parameter (nf_enotnc4 = -111)
      integer nf_estrictnc3     ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
      parameter (nf_estrictnc3 = -112)
      integer nf_enotnc3        ! Attempting netcdf-3 operation on netcdf-4 file.   
      parameter (nf_enotnc3 = -113)
      integer nf_enopar         ! Parallel operation on file opened for non-parallel access.   
      parameter (nf_enopar = -114)
      integer nf_eparinit       ! Error initializing for parallel access.   
      parameter (nf_eparinit = -115)
      integer nf_ebadgrpid      ! Bad group ID.   
      parameter (nf_ebadgrpid = -116)
      integer nf_ebadtypid      ! Bad type ID.   
      parameter (nf_ebadtypid = -117)
      integer nf_etypdefined    ! Type has already been defined and may not be edited. 
      parameter (nf_etypdefined = -118)
      integer nf_ebadfield      ! Bad field ID.   
      parameter (nf_ebadfield = -119)
      integer nf_ebadclass      ! Bad class.   
      parameter (nf_ebadclass = -120)
      integer nf_emaptype       ! Mapped access for atomic types only.   
      parameter (nf_emaptype = -121)
      integer nf_elatefill      ! Attempt to define fill value when data already exists. 
      parameter (nf_elatefill = -122)
      integer nf_elatedef       ! Attempt to define var properties, like deflate, after enddef. 
      parameter (nf_elatedef = -123)
      integer nf_edimscale      ! Probem with 1 dimscales. 
      parameter (nf_edimscale = -124)
      integer nf_enogrp       ! No group found.
      parameter (nf_enogrp = -125)


!     New functions.

!     Parallel I/O.
      integer nf_create_par
      external nf_create_par

      integer nf_open_par
      external nf_open_par

      integer nf_var_par_access
      external nf_var_par_access

!     Functions to handle groups.
      integer nf_inq_ncid
      external nf_inq_ncid

      integer nf_inq_grps
      external nf_inq_grps

      integer nf_inq_grpname
      external nf_inq_grpname

      integer nf_inq_grpname_full
      external nf_inq_grpname_full

      integer nf_inq_grpname_len
      external nf_inq_grpname_len

      integer nf_inq_grp_parent
      external nf_inq_grp_parent

      integer nf_inq_grp_ncid
      external nf_inq_grp_ncid

      integer nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid

      integer nf_inq_varids
      external nf_inq_varids

      integer nf_inq_dimids
      external nf_inq_dimids

      integer nf_def_grp
      external nf_def_grp

!     New rename grp function

      integer nf_rename_grp
      external nf_rename_grp

!     New options for netCDF variables.
      integer nf_def_var_deflate
      external nf_def_var_deflate

      integer nf_inq_var_deflate
      external nf_inq_var_deflate

      integer nf_def_var_szip
      external nf_def_var_szip

      integer nf_inq_var_szip
      external nf_inq_var_szip

      integer nf_def_var_fletcher32
      external nf_def_var_fletcher32

      integer nf_inq_var_fletcher32
      external nf_inq_var_fletcher32

      integer nf_def_var_chunking
      external nf_def_var_chunking

      integer nf_inq_var_chunking
      external nf_inq_var_chunking

      integer nf_def_var_fill
      external nf_def_var_fill

      integer nf_inq_var_fill
      external nf_inq_var_fill

      integer nf_def_var_endian
      external nf_def_var_endian

      integer nf_inq_var_endian
      external nf_inq_var_endian

      integer nf_def_var_filter
      external nf_def_var_filter

      integer nf_inq_var_filter
      external nf_inq_var_filter

!     User defined types.
      integer nf_inq_typeids
      external nf_inq_typeids

      integer nf_inq_typeid
      external nf_inq_typeid

      integer nf_inq_type
      external nf_inq_type

      integer nf_inq_user_type
      external nf_inq_user_type

!     User defined types - compound types.
      integer nf_def_compound
      external nf_def_compound

      integer nf_insert_compound
      external nf_insert_compound

      integer nf_insert_array_compound
      external nf_insert_array_compound

      integer nf_inq_compound
      external nf_inq_compound

      integer nf_inq_compound_name
      external nf_inq_compound_name

      integer nf_inq_compound_size
      external nf_inq_compound_size

      integer nf_inq_compound_nfields
      external nf_inq_compound_nfields

      integer nf_inq_compound_field
      external nf_inq_compound_field

      integer nf_inq_compound_fieldname
      external nf_inq_compound_fieldname

      integer nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex

      integer nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset

      integer nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype

      integer nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims

      integer nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes

!     User defined types - variable length arrays.
      integer nf_def_vlen
      external nf_def_vlen

      integer nf_inq_vlen
      external nf_inq_vlen

      integer nf_free_vlen
      external nf_free_vlen

!     User defined types - enums.
      integer nf_def_enum
      external nf_def_enum

      integer nf_insert_enum
      external nf_insert_enum

      integer nf_inq_enum
      external nf_inq_enum

      integer nf_inq_enum_member
      external nf_inq_enum_member

      integer nf_inq_enum_ident
      external nf_inq_enum_ident

!     User defined types - opaque.
      integer nf_def_opaque
      external nf_def_opaque

      integer nf_inq_opaque
      external nf_inq_opaque

!     Write and read attributes of any type, including user defined
!     types.
      integer nf_put_att
      external nf_put_att
      integer nf_get_att
      external nf_get_att

!     Write and read variables of any type, including user defined
!     types.
      integer nf_put_var
      external nf_put_var
      integer nf_put_var1
      external nf_put_var1
      integer nf_put_vara
      external nf_put_vara
      integer nf_put_vars
      external nf_put_vars
      integer nf_get_var
      external nf_get_var
      integer nf_get_var1
      external nf_get_var1
      integer nf_get_vara
      external nf_get_vara
      integer nf_get_vars
      external nf_get_vars

!     For helping F77 users with VLENs.
      integer nf_get_vlen_element
      external nf_get_vlen_element
      integer nf_put_vlen_element
      external nf_put_vlen_element

!     For dealing with file level chunk cache.
      integer nf_set_chunk_cache
      external nf_set_chunk_cache
      integer nf_get_chunk_cache
      external nf_get_chunk_cache

!     For dealing with per variable chunk cache.
      integer nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer nf_get_var_chunk_cache
      external nf_get_var_chunk_cache

!     NetCDF-2.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!      
! functions in the fortran interface
!
      integer nccre
      integer ncopn
      integer ncddef
      integer ncdid
      integer ncvdef
      integer ncvid
      integer nctlen
      integer ncsfil

      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil


      integer ncrdwr
      integer nccreat
      integer ncexcl
      integer ncindef
      integer ncnsync
      integer nchsync
      integer ncndirty
      integer nchdirty
      integer nclink
      integer ncnowrit
      integer ncwrite
      integer ncclob
      integer ncnoclob
      integer ncglobal
      integer ncfill
      integer ncnofill
      integer maxncop
      integer maxncdim
      integer maxncatt
      integer maxncvar
      integer maxncnam
      integer maxvdims
      integer ncnoerr
      integer ncebadid
      integer ncenfile
      integer nceexist
      integer nceinval
      integer nceperm
      integer ncenotin
      integer nceindef
      integer ncecoord
      integer ncemaxds
      integer ncename
      integer ncenoatt
      integer ncemaxat
      integer ncebadty
      integer ncebadd
      integer ncests
      integer nceunlim
      integer ncemaxvs
      integer ncenotvr
      integer nceglob
      integer ncenotnc
      integer ncfoobar
      integer ncsyserr
      integer ncfatal
      integer ncverbos
      integer ncentool


!
! netcdf data types:
!
      integer ncbyte
      integer ncchar
      integer ncshort
      integer nclong
      integer ncfloat
      integer ncdouble

      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)

!     
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!     

!     read/write, 0 => readonly 
      parameter(ncrdwr = 1)
!     in create phase, cleared by ncendef 
      parameter(nccreat = 2)
!     on create destroy existing file 
      parameter(ncexcl = 4)
!     in define mode, cleared by ncendef 
      parameter(ncindef = 8)
!     synchronise numrecs on change (x'10')
      parameter(ncnsync = 16)
!     synchronise whole header on change (x'20')
      parameter(nchsync = 32)
!     numrecs has changed (x'40')
      parameter(ncndirty = 64)  
!     header info has changed (x'80')
      parameter(nchdirty = 128)
!     prefill vars on endef and increase of record, the default behavior
      parameter(ncfill = 0)
!     do not fill vars on endef and increase of record (x'100')
      parameter(ncnofill = 256)
!     isa link (x'8000')
      parameter(nclink = 32768)

!     
!     'mode' arguments for nccreate and ncopen
!     
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)

!     
!     'size' argument to ncdimdef for an unlimited dimension
!     
      integer ncunlim
      parameter(ncunlim = 0)

!     
!     attribute id to put/get a global attribute
!     
      parameter(ncglobal  = 0)

!     
!     advisory maximums:
!     
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
!     not enforced 
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)

!     
!     global netcdf error status variable
!     initialized in error.c
!     

!     no error 
      parameter(ncnoerr = nf_noerr)
!     not a netcdf id 
      parameter(ncebadid = nf_ebadid)
!     too many netcdfs open 
      parameter(ncenfile = -31)   ! nc_syserr
!     netcdf file exists && ncnoclob
      parameter(nceexist = nf_eexist)
!     invalid argument 
      parameter(nceinval = nf_einval)
!     write to read only 
      parameter(nceperm = nf_eperm)
!     operation not allowed in data mode 
      parameter(ncenotin = nf_enotindefine )   
!     operation not allowed in define mode 
      parameter(nceindef = nf_eindefine)   
!     coordinates out of domain 
      parameter(ncecoord = nf_einvalcoords)
!     maxncdims exceeded 
      parameter(ncemaxds = nf_emaxdims)
!     string match to name in use 
      parameter(ncename = nf_enameinuse)   
!     attribute not found 
      parameter(ncenoatt = nf_enotatt)
!     maxncattrs exceeded 
      parameter(ncemaxat = nf_emaxatts)
!     not a netcdf data type 
      parameter(ncebadty = nf_ebadtype)
!     invalid dimension id 
      parameter(ncebadd = nf_ebaddim)
!     ncunlimited in the wrong index 
      parameter(nceunlim = nf_eunlimpos)
!     maxncvars exceeded 
      parameter(ncemaxvs = nf_emaxvars)
!     variable not found 
      parameter(ncenotvr = nf_enotvar)
!     action prohibited on ncglobal varid 
      parameter(nceglob = nf_eglobal)
!     not a netcdf file 
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname) 
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)

!     
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!     
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)

!
!     default fill values.  these must be the same as in the c interface.
!
      integer filbyte
      integer filchar
      integer filshort
      integer fillong
      real filfloat
      doubleprecision fildoub

      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690e+36)
      parameter (fildoub = 9.9692099683868690e+36)

   character (len=200), intent(in):: input_file       ! 1 file name.
   character (len=10), intent(in) :: var              ! Variable to search for.
   integer, intent(out)           :: dim1             ! Dimensions of field. 
   integer, intent(out)           :: dim2             ! Dimensions of field. 
   integer, intent(out)           :: dim3             ! Dimensions of field. 
   real, intent(out)              :: ds               ! Grid resolution.
   integer, intent(out), optional :: mp_physics       ! microphysics option

   character (len=80)             :: att_name         ! Attribute name.
   integer                        :: i                ! Loop counters.
   integer                        :: attlen           ! Length of attribute.
   integer                        :: cdfid            ! 1 file id.
   integer                        :: rcode            ! Return code (0=ok).
   integer                        :: length           ! Length of filename.
   integer                        :: id_var           ! 1 variable ID. 
   integer                        :: ivtype           ! 4=integer, 5=real, 6=d.p.
   integer                        :: ndims            !
   integer                        :: natts            ! Number of attributes.

   integer                        :: dimids(10)       !
   integer                        :: dims(1:4)        ! Dimensions of field. 
   real (kind=4), allocatable     :: value_real(:)    ! real attribute value. 

   ! if (trace_use) call da_trace_entry("da_stage0_initialize")

   ! Check file exists:
   length = len_trim(input_file)
   rcode = nf_open(input_file(1:length), NF_NOwrite, cdfid)
   if (rcode /= 0) then
      write(unit=message(1),fmt='(A,A)') &
         ' Error opening netcdf file ', input_file(1:length)
      call da_error("da_stage0_initialize.inc",42,message(1:1))
   end if

   ! Check variable is in file:
   rcode = nf_inq_varid (cdfid, var, id_var)
   if (rcode /= 0) then 
      write(unit=message(1),fmt='(A,A)') &
         var, ' variable is not in input file'
      call da_error("da_stage0_initialize.inc",50,message(1:1))
   end if

   ! Get metadata (time, dimensions, global attributes):
   rcode = nf_inq_var(cdfid, id_var, var, ivtype, ndims, dimids, natts)
   do i = 1, ndims
     rcode = nf_inq_dimlen(cdfid, dimids(i), dims(i))
   end do
   dim1 = dims(1)
   dim2 = dims(2)
   dim3 = dims(3)

   ! Get resolution:
   att_name = 'DX'  ! Assume dx= dy.
   rcode = nf_inq_att(cdfid, nf_global, att_name, ivtype, attlen)
   allocate(value_real(1:attlen))
   rcode = NF_GET_ATT_real(cdfid, nf_global, att_name, value_real)
   ds = value_real(1)
   deallocate(value_real)

   ! get MP_PHYSICS option
   if ( present(mp_physics) ) then
      mp_physics = 0 !initialize
      ! get mp_physics from the global attribute
      rcode = NF_GET_ATT_int(cdfid, nf_global, 'MP_PHYSICS', mp_physics)
   end if

   rcode = nf_close(cdfid)

   ! if (trace_use) call da_trace_exit("da_stage0_initialize")

end subroutine da_stage0_initialize



   
subroutine da_transform_vptovv (evec, eval, vertical_wgt, vp, vv, mz, &
   kds,kde, ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte)

   !---------------------------------------------------------------------------
   ! Purpose: Transform from fields on vertical levels to fields on vertical 
   ! EOFS.
   !
   ! Method: Perform vv(i,j,m) = L^{-1/2} E^{T} vp(i,j,k) transform.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: mz                         ! # vertical modes.
   integer, intent(in)    :: kds,kde  ! domain dims.
   integer, intent(in)    :: ims,ime, jms,jme, kms,kme  ! memory dims.
   integer, intent(in)    :: its,ite, jts,jte, kts,kte  ! tile   dims
   real*8,  intent(in)    :: evec(jts:jte,kds:kde,1:mz) ! Eigenvectors.
   real*8,  intent(in)    :: eval(jts:jte,1:mz)         ! Eigenvalues.
   real,    intent(in)    :: vertical_wgt(ims:ime,jms:jme,kms:kme) ! Weighting.
   real,    intent(inout) :: vp(ims:ime,jms:jme,kms:kme)! CV in level space.
   real,    intent(out)   :: vv(ims:ime,jms:jme,1:mz)   ! CV in EOF space.
   
   integer :: i, j, m                    ! Loop counters.
   real    :: ETVp                       ! E(k,m)^{T}*vp(i,j,k)

   ! Do not add trace, as routine used by gen_be
   ! if (trace_use) call da_trace_entry("da_transform_vptovv")
   
   !-------------------------------------------------------------------
   ! [1.0] Apply inner-product weighting if vertical_ip /= vertical_ip_0
   !------------------------------------------------------------------- 

   if (vertical_ip /= vertical_ip_0) then
      vp(its:ite,jts:jte,kts:kte) = vp(its:ite,jts:jte,kts:kte) * &
                                    vertical_wgt(its:ite,jts:jte,kts:kte)
   end if

   !-------------------------------------------------------------------
   ! [2.0] Perform vv(i,j,m) = L^{-1/2} E^T vp(i,j,k) transform:
   !-------------------------------------------------------------------

   do m = 1, mz
      do j = jts, jte
         do i = its, ite
            ETVp = sum(evec(j,kts:kte,m) * vp(i,j,kts:kte))
            vv(i,j,m) = ETVp / eval(j,m)
         end do
      end do
   end do

end subroutine da_transform_vptovv


subroutine da_eof_decomposition (kz, bx, e, l)
   
   !---------------------------------------------------------------------------
   ! Purpose: Compute eigenvectors E and eigenvalues L of vertical covariance 
   !          matrix
   !          B_{x} defined by equation:  E^{T} B_{x} E = L, given input kz x kz 
   !          BE field.
   !---------------------------------------------------------------------------
   
   implicit none

   integer, intent(in)  :: kz               ! Dimension of error matrix. 
   real,    intent(in)  :: bx(1:kz,1:kz)    ! Vert. background error.
   real*8,  intent(out) :: e(1:kz,1:kz)     ! Eigenvectors of Bx.
   real*8,  intent(out) :: l(1:kz)          ! Eigenvalues of Bx.

   integer :: work             ! Size of work array.
   integer :: m                ! Loop counters
   integer :: info             ! Info code.

   real*8  :: work_array(1:3*kz-1)
   real*8  :: ecopy(1:kz,1:kz)
   real*8  :: lcopy(1:kz)   

   if (trace_use) call da_trace_entry("da_eof_decomposition")    

   !-------------------------------------------------------------------------
   ! [5.0]: Perform global eigenvalue decomposition using LAPACK software:
   !-------------------------------------------------------------------------
   
   work = 3 * kz - 1   
   ecopy(1:kz,1:kz) = bx(1:kz,1:kz)
   lcopy(1:kz) = 0.0

   call dsyev( 'V', 'U', kz, ecopy, kz, lcopy, work_array, work, info )
   
   if ( info /= 0 ) then
      write(unit=message(1),fmt='(A,I4)') &
         "Error in decomposition, info = ", info
      call da_error("da_eof_decomposition.inc",40,message(1:1))
   end if
   
   ! Swap order of eigenvalues, vectors so 1st is one with most variance:
   
   do m = 1, kz
      l(m) = lcopy(kz+1-m)
      e(1:kz,m) = ecopy(1:kz,kz+1-m)
   end do  

   if (trace_use) call da_trace_exit("da_eof_decomposition")    
   
end subroutine da_eof_decomposition


subroutine da_eof_decomposition_test (kz, bx, e, l)
   
   !------------------------------------------------------------------------------
   ! Purpose: 
   ! [1] Print eigenvalues:
   ! [2] Test orthogonality of eigenvectors - sum_k (e_m(k) e_n(k)) = delta_mn:
   ! [3] Test eigenvectors completeness - sum_m (e_m(k1) e_m(k2)) = delta_k1k2:
   ! [4] Check B correctness: B = E*L*E^T
   !------------------------------------------------------------------------------
   
   implicit none

   integer, intent(in) :: kz               ! Dimension of BE matrix   
   real,    intent(in) :: bx(1:kz,1:kz)    ! Global vert. background error.
   real*8,  intent(in) :: e(1:kz,1:kz)     ! Eigenvectors of Bx.
   real*8,  intent(in) :: l(1:kz)          ! Eigenvalues of Bx.
   
   integer                  :: k, k1, k2, m     ! Loop counters
   real                     :: tot_variance     ! Total variance.
   real                     :: cum_variance     ! Cumulative variance.
   real                     :: max_off_diag     ! Maximum off-diagonal.

   real                     :: work(1:kz,1:kz)  ! 2D Work matrix.
   real                     :: bc(1:kz,1:kz)    ! 2D Work matrix.
   logical                  :: array_mask(1:kz) ! Array mask for MAXVAL.

   if (trace_use) call da_trace_entry("da_eof_decomposition_test")

   !------------------------------------------------------------------------- 
   ! [1] Print eigenvalues:
   !-------------------------------------------------------------------------

   tot_variance = sum(l(1:kz))
   cum_variance = 0.0
   
   write(unit=stdout,fmt='(A)')'  Mode    Eigenvalue     Cumulative Variance      e(k,k)'

   do k = 1, kz
      cum_variance = cum_variance + l(k)
      write(unit=stdout,fmt='(I4,4x,e12.4,10x,f8.4,4x,e12.4)') &
            k, l(k), cum_variance / tot_variance, e(k,k)
   end do

   write(unit=stdout,fmt=*)
   
   call da_array_print( 1, e, 'Global Eigenvectors' )

   !-------------------------------------------------------------------------
   ! [2] Test orthogonality of eigenvectors - sum_k (e_m(k) e_n(k)) = delta_mn:
   !-------------------------------------------------------------------------
   
   write(unit=stdout,fmt='(A)')' Eigenvector orthogonality check:'
   write(unit=stdout,fmt='(A)')' Mode     Diagonal         Maximum off-diagonal'

   do k1 = 1, kz
      do k2 = 1, kz
         work(k1,k2) = sum(e(1:kz,k1) * e(1:kz,k2))
      end do
   
      array_mask(1:kz) =.true.
      array_mask(k1) = .false.
      max_off_diag = maxval(abs(work(k1,:)),mask=array_mask(:))
      write(unit=stdout,fmt='(I4,4x,1pe12.4,10x,1pe12.4)')k1, work(k1,k1), max_off_diag
   end do
   write(unit=stdout,fmt=*)

   !-------------------------------------------------------------------------   
   ! [3] Test eigenvectors completeness - sum_m (e_m(k1) e_m(k2)) = delta_k1k2:
   !-------------------------------------------------------------------------   
   
   write(unit=stdout,fmt='(A)')' Eigenvector completeness check:'
   write(unit=stdout,fmt='(A)')' Level    Diagonal         Maximum off-diagonal'

   do k1 = 1, kz
      do k2 = 1, kz
         work(k1,k2) = sum(e(k1,1:kz) * e(k2,1:kz))
      end do
   
      array_mask(1:kz) =.true.
      array_mask(k1) = .false.
      max_off_diag = maxval(abs(work(k1,:)),mask=array_mask(:))
      write(unit=stdout,fmt='(I4,4x,1pe12.4,10x,1pe12.4)')k1, work(k1,k1), max_off_diag
   end do
   write(unit=stdout,fmt=*)

   !-------------------------------------------------------------------------
   ! [4]  check B correctness: B = E*L*E^T
   !-------------------------------------------------------------------------

   write(unit=stdout,fmt='(a/a)') &
        'real and Calculated B (diagonal)', &
        'lvl                 real-B                    Calculated-B'

   do k=1,kz
      do m=1,kz
         work(k,m)=l(k)*e(m,k)
         bc(k,m)=0.0
      end do
   end do
   
   do k1=1,kz
      do k2=1,kz
         do m=1,kz
            bc(k1,k2)=bc(k1,k2)+e(k1,m)*work(m,k2)
         end do
      end do

      write(unit=stdout,fmt='(I5,2F20.5)') k1, bx(k1,k1), bc(k1,k1)
   end do

   do k2=1,kz
      write(unit=stdout, fmt='(a,i4/a)') &
           'real and Calculated B (off diagonal):', k2, &
           'lvl                 real-B                    Calculated-B'

      do k1=1,kz
        write(unit=stdout,fmt='(I5,2F20.5)') k1, bx(k1,k2), bc(k1,k2)
      end do
   end do

   if (trace_use) call da_trace_exit("da_eof_decomposition_test")
   
end subroutine da_eof_decomposition_test


subroutine da_perform_2drf(ni, nj, num_passes, rf_scale, field)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: ni               ! Array dimension 1.
   integer, intent(in)    :: nj               ! Array dimension 2.
   integer, intent(in)    :: num_passes       ! Number of passes of RF.
   real,    intent(in)    :: rf_scale         ! Recursive filter scaling parameter.
   real*8,  intent(inout) :: field(1:ni,1:nj) ! Field to be filtered.

   integer               :: i, j, pass       ! Loop counters.
   real*8                :: e, alpha         ! Recursive filter parameters.
   real                  :: mean_field       ! Mean field.

   if (trace_use) call da_trace_entry("da_perform_2drf")

   e = 0.25 * num_passes / (rf_scale * rf_scale)
   alpha = 1 + e - sqrt(e * (e + 2.0))

   mean_field = sum(field(1:ni,1:nj)) / real(ni*nj)

   do pass = 1, num_passes
      ! Perform filter in I-direction:
      do j = 1, nj
         call da_recursive_filter_1d(pass, alpha, field(1:ni,j), ni)
      end do

      ! Perform filter in J-direction:
      do i = 1, ni
         call da_recursive_filter_1d(pass, alpha, field(i,1:nj), nj)
      end do
   end do

   if (trace_use) call da_trace_exit("da_perform_2drf")

end subroutine da_perform_2drf


subroutine da_recursive_filter_1d(pass, alpha, field, n)

   !---------------------------------------------------------------------------
   ! Purpose: Perform one pass of recursive filter on 1D array.
   !
   ! Method:  Perform right-moving filter followed by left-moving filter.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: pass           ! Current pass of filter.
   real*8,  intent(in)    :: alpha          ! Alpha coefficient for RF.
   real*8,  intent(inout) :: field(:)       ! Array to be filtered.
   integer, intent(in)    :: n              ! Size of field array.

   integer :: j              ! Loop counter.
   real    :: one_alpha      ! 1 - alpha.
   real    :: a(1:n)         ! Input field.
   real    :: b(1:n)         ! Field after left-right pass.
   real    :: c(1:n)         ! Field after right-left pass.

   if (trace_use_dull) call da_trace_entry("da_recursive_filter_1d")
   
   !-------------------------------------------------------------------------
   ! [1.0] Initialise:
   !-------------------------------------------------------------------------

   one_alpha = 1.0 - alpha
   
   a(1:n) = field(1:n)

   !-------------------------------------------------------------------------
   ! [2.0] Perform right-moving filter:
   !-------------------------------------------------------------------------

   ! use turning conditions as in the appendix of Hayden & Purser (1995):

   if (pass == 1) then
      b(1) = one_alpha * a(1)
   else if (pass == 2) then
      b(1) = a(1) / (1.0 + alpha)
   else
      b(1) = one_alpha * (a(1) - alpha**3 * a(2)) / (1.0 - alpha**2)**2
   end if

   ! [2.2] Perform pass left to right:

   do j = 2, n
      b(j) = alpha * b(j-1) + one_alpha * a(j)
   end do

   !-------------------------------------------------------------------------
   ! [3.0] Perform left-moving filter:
   !-------------------------------------------------------------------------

   ! use turning conditions as in the appendix of Hayden & Purser (1995):

   if (pass == 1) then
      c(n) = b(n) / (1.0 + alpha)
   else
      c(n) = one_alpha * (b(n) - alpha**3 * b(n-1)) / (1.0 - alpha**2)**2
   end if

   ! [3.2] Perform pass left to right:

   do j = n-1, 1, -1
      c(j) = alpha * c(j+1) + one_alpha * b(j)
   end do
        
   field(1:n) = c(1:n)

   if (trace_use_dull) call da_trace_exit("da_recursive_filter_1d")
   
end subroutine da_recursive_filter_1d



end module da_gen_be

subroutine wrf_abort
   stop
end subroutine wrf_abort

   LOGICAL FUNCTION wrf_dm_on_monitor()
      IMPLICIT NONE
      wrf_dm_on_monitor = .TRUE.
      RETURN
   END FUNCTION wrf_dm_on_monitor

