












module da_par_util

   
   
   
   
   

   use da_control, only: use_rf
   use module_domain, only : domain, xpose_type

   use module_dm, only : local_communicator_x, &
      local_communicator_y, ntasks_x, ntasks_y, data_order_xyz
   use da_par_util1, only : true_mpi_real, true_mpi_complex




   use da_define_structures, only : be_subtype, &
      x_type, vp_type, residual_synop_type, residual_sound_type, iv_type, &
      y_type, count_obs_number_type, maxmin_field_type

   use da_control, only : trace_use,num_ob_indexes, myproc, root, comm, ierr, &
      rootproc, num_procs, stdout, print_detail_parallel, its,ite, jts, jte, &
      kts,kte,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,ips,ipe,jps,jpe, &
      kps, kpe, grid_stagger, grid_ordering, trace_use_dull, &
      sound, synop, pilot, satem, geoamv, polaramv, airep, gpspw, gpsref, &
      metar, ships, ssmi_rv, ssmi_tb, ssmt1, ssmt2, qscat, profiler, buoy, bogus, &
      pseudo, radar, radiance, airsr, sonde_sfc, trace_use_frequent, &
      its_int,ite_int,jts_int,jte_int,kts_int,kte_int, &
      ims_int,ime_int,jms_int,jme_int,kms_int,kme_int, &
      ids_int,ide_int,jds_int,jde_int,kds_int,kde_int, &
      ips_int,ipe_int,jps_int,jpe_int,kps_int,kpe_int, &
      anal_type_hybrid_dual_res
   use da_reporting, only : da_error
   use da_tracing, only : da_trace_entry, da_trace_exit
   use da_wrf_interfaces, only : &
      wrf_dm_xpose_z2x, wrf_dm_xpose_x2y, wrf_dm_xpose_y2x, wrf_dm_xpose_x2z, &
      wrf_dm_xpose_z2y, wrf_dm_xpose_y2z, wrf_patch_to_global_real, wrf_debug

   implicit none

   include 'mpif.h'

   ! generic vector type
   type generic_vector_type
      real, pointer :: ptr(:)
   end type generic_vector_type

   ! Generic residual type contains lists of vectors and scalars.  
   ! Implementation notes:   
   !  - Vector values are always assumed to be stored by reference in
   !    self%values(1:UBOUND(self%values,1))%ptr.  The size of each vector
   !    is the same.  These pointers are used to reference arrays that 
   !    will be deallocated elsewhere, they are not deallocated by 
   !    destructor residual_generic_free().  When 
   !    UBOUND(self%values,1)==0, there are no vector values.
   !  - Scalars are copied into and out of self%values(0)%ptr.  This pointer
   !    is be allocated and deallocated by this class.  When
   !    LBOUND(self%values,1)==1, there are no scalar values.  The
   !    duplication of scalars is unavoidable without changing
   !    implementations of the specific residual_*_type classes.
   !  - It would be better to store references to model_loc_type and
   !    info_type objects here rather than copying proc_domain and
   !    obs_global_index.  Make this change later when there is
   !    time to modify other objects to also store references(there should 
   !    be only one object responsible for each model_loc_type and info_type
   !    object, but in Fortran90 this must be implemented as a pointer).
   type residual_generic_type
      logical                             :: proc_domain
      integer                             :: obs_global_index
      type(generic_vector_type), pointer :: values(:) ! vectors & scalars
   end type residual_generic_type

   ! template for allocating memory
   type residual_template_type
      integer                             :: lbnd  ! lower bound
      integer                             :: ubnd  ! upper bound
   end type residual_template_type

   ! single-obs generic y_type
   type y_facade_type
      integer                               :: num_obs
      integer                               :: num_obs_glo
      type(residual_generic_type), pointer :: obs(:)
   end type y_facade_type



   interface da_patch_to_global
      module procedure da_patch_to_global_2d
      module procedure da_patch_to_global_3d
   end interface

   contains

subroutine da_cv_to_vv (cv_size, rcv, mzs, vv)

   !---------------------------------------------------------------------------
   ! Purpose: Fill (local) vv arrays from 1D (local) cv array.
   !---------------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: cv_size     ! Length of CV array.
   real,           intent(in)    :: rcv(1:cv_size) ! Control variables v.
   integer,        intent(in)    :: mzs(:) ! Background error structure levels.
   type (vp_type), intent(inout) :: vv          ! Grdipt/EOF cv_array.

   integer   :: is,ie       ! Local grid range in y coordinate.
   integer   :: js,je       ! Local grid range in x coordinate.
   integer   :: ix,jy       ! Local grid horizontal dimensions.
   integer   :: mz          ! Max vertical coordinate for v1 through v5 arrays.
   integer   :: ne          ! Ensemble size.
   integer   :: cv_s,cv_e   ! Starting and ending indices into CV array.
   integer   :: ijm         ! Size of interior of v1 through v5 arrays.
   integer   :: ijmn        ! Size of interior of alpha cv arrays.

   if (trace_use) call da_trace_entry("da_cv_to_vv")

   if( use_rf )then
      is = its
      ie = ite
      js = jts
      je = jte
   else
      call da_error("da_cv_to_vv.inc",31,(/"This subroutine should not be called for use_rf = .false."/))
   endif
   ix = ie-is+1
   jy = je-js+1
   cv_e = 0

   !--------------------------------------------------------------------------
   ! [1] Transfer components of Jb control variable:
   !--------------------------------------------------------------------------

   ! Fill v1
   mz = mzs(1)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      vv % v1(is:ie,js:je,1:mz) = RESHAPE(rcv(cv_s:cv_e),(/ix, jy, mz/))
   end if

   ! Fill v2
   mz = mzs(2)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      vv % v2(is:ie,js:je,1:mz) = RESHAPE(rcv(cv_s:cv_e),(/ix, jy, mz/))
   end if

   ! Fill v3
   mz = mzs(3)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      vv % v3(is:ie,js:je,1:mz) = RESHAPE(rcv(cv_s:cv_e),(/ix, jy, mz/))
   end if

   ! Fill v4
   mz = mzs(4)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      vv % v4(is:ie,js:je,1:mz) = RESHAPE(rcv(cv_s:cv_e),(/ix, jy, mz/))
   end if

   ! Fill v5
   mz = mzs(5)
   if (mz == 1) then ! Can only be 0 or 1 (2D ps_u field)
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      vv % v5(is:ie,js:je,1:mz) = RESHAPE(rcv(cv_s:cv_e),(/ix, jy,mz/))
   end if
   !--------------------------------------------------------------------------
   ! [2] Transfer components of Je control variable:
   !--------------------------------------------------------------------------
   mz = mzs(6)
   ne = mzs(7)
   if ( ne > 0 ) then
      ix = ite_int - its_int + 1
      jy = jte_int - jts_int + 1
      ijmn = ix * jy * mz * ne
      cv_s = cv_e + 1
      cv_e = cv_s + ijmn - 1
!     vv % alpha(is:ie,js:je,1:mz,1:ne) = RESHAPE(rcv(cv_s:cv_e),(/ix, jy, mz, ne/))
      vv % alpha(its_int:ite_int,jts_int:jte_int,1:mz,1:ne) = RESHAPE(rcv(cv_s:cv_e),(/ix, jy, mz, ne/))
   end if

   if (trace_use) call da_trace_exit("da_cv_to_vv")

end subroutine da_cv_to_vv


subroutine da_vv_to_cv(vv, xp, mzs, cv_size, rcv)

   !---------------------------------------------------------------------------
   ! Purpose: Fill (local) 1D cv array from 3D (local) vv arrays.
   !---------------------------------------------------------------------------

   implicit none

   type (vp_type), intent(in)    :: vv          ! Grdipt/EOF rcv.
   type (xpose_type), intent(in) :: xp          ! Dimensions and xpose buffers.
   integer,        intent(in)    :: mzs(:)      ! Background error structure levels.
   integer,        intent(in)    :: cv_size     ! Length of CV array.
   real,           intent(inout) :: rcv(1:cv_size) ! Control variables v.

   integer   :: is,ie       ! Local grid range in y coordinate.
   integer   :: js,je       ! Local grid range in x coordinate.
   integer   :: ix,jy       ! Local grid horizontal dimensions.
   integer   :: mz          ! Max vertical coordinate for v1 through v5 arrays.
   integer   :: ne          ! Ensemble size.
   integer   :: cv_s,cv_e   ! Starting and ending indices into CV array.
   integer   :: ijm         ! Size of interior of v1 through v5 arrays.
   integer   :: ijmn        ! Size of interior of alpha cv arrays.

   integer   :: i,j,k,ijk,m,n

   if (trace_use) call da_trace_entry("da_vv_to_cv")

   if( use_rf )then
      is = xp % its
      ie = xp % ite
      js = xp % jts
      je = xp % jte
   else
      call da_error("da_vv_to_cv.inc",36,(/"This subroutine should not be called for use_rf = .false."/))
   endif
   ix = ie-is+1
   jy = je-js+1
   cv_e = 0
   ijk=0

   ! Store v1
   mz = mzs(1)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v1(i,j,k)
            end do
        end do
     end do
   end if
   ! Store v2
   mz = mzs(2)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v2(i,j,k)
            end do
         end do
      end do
   end if
   ! Store v3
   mz = mzs(3)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1

      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v3(i,j,k)
            end do
         end do
      end do
   end if
   ! Store v4
   mz = mzs(4)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v4(i,j,k)
            end do
         end do
      end do
   end if
   ! Store v5
   mz = mzs(5)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v5(i,j,k)
            end do
         end do
      end do
   end if
   ! Store alpha:
   mz = mzs(6)
   ne = mzs(7)
   if ( ne > 0 ) then
      ix = ite_int - its_int + 1
      jy = jte_int - jts_int + 1
      ijmn = ix * jy * mz * ne
      cv_s = cv_e + 1
      cv_e = cv_s + ijmn - 1
      do n = 1, ne
         do m = 1, mz
            do j = jts_int, jte_int !js, je
               do i = its_int, ite_int ! is, ie
                  ijk = ijk + 1
                  rcv(ijk) = vv%alpha(i,j,m,n)
               end do
            end do
         end do
      end do
   end if

   if (trace_use) call da_trace_exit("da_vv_to_cv")

end subroutine da_vv_to_cv


subroutine da_alloc_and_copy_be_arrays (vg , v, jts,jte, kts,kte)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (be_subtype), intent(inout) :: vg       ! Global backgrd error struct.
   type (be_subtype), intent(inout) :: v        ! Local backgrd error struct.
   integer, intent(in)              :: jts,jte  ! Tile dimension.
   integer, intent(in)              :: kts,kte  ! Tile dimension.

   ! Allocate local-grid structure.
   v % mz = vg % mz
   if (v % mz > 0) then
      allocate  (v % val(jts:jte,1:v % mz))
      allocate  (v % evec(jts:jte,kts:kte,1:v % mz))
      allocate  (v % rf_alpha(1:v % mz))
   end if

   ! Make local copies from global arrays.
   if (v % mz > 0) then
       v % val(jts:jte,1:v % mz) = vg % val(jts:jte,1:v % mz)
       v % evec(jts:jte,kts:kte,1:v % mz) = vg % evec(jts:jte,kts:kte,1:v % mz)
       v % rf_alpha(1:v % mz) = vg % rf_alpha(1:v % mz)
   end if

   ! Deallocate global arrays.
   if (v % mz > 0) then 
      deallocate  (vg % val )
      deallocate  (vg % evec)
      deallocate  (vg % rf_alpha)
   end if

end subroutine da_alloc_and_copy_be_arrays 


subroutine da_copy_dims(grid)

   !---------------------------------------------------------------------------
   ! Purpose: Copy dimensioning information from grid structure.
   !---------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)         :: grid

   if (trace_use_dull) call da_trace_entry("da_copy_dims")

   ! De-reference dimension information stored in the grid data structure.

   ids = grid%sd31 
   ide = grid%ed31 - 1
   jds = grid%sd32 
   jde = grid%ed32 - 1
   kds = grid%sd33 
   kde = grid%ed33 - 1

   ims = grid%sm31 
   ime = grid%em31 
   jms = grid%sm32 
   jme = grid%em32 
   kms = grid%sm33 
   kme = grid%em33 

   ips = grid%sp31 
   ipe = grid%ep31 
   jps = grid%sp32 
   jpe = grid%ep32 
   kps = grid%sp33 
   kpe = grid%ep33 

   ! Indices for yz decomposition

   grid%xp%idsx = grid%sd31
   grid%xp%idex = grid%ed31 - 1
   grid%xp%jdsx = grid%sd32
   grid%xp%jdex = grid%ed32 - 1
   grid%xp%kdsx = grid%sd33
   grid%xp%kdex = grid%ed33 - 1

   grid%xp%imsx = grid%sm31x
   grid%xp%imex = grid%em31x
   grid%xp%jmsx = grid%sm32x
   grid%xp%jmex = grid%em32x
   grid%xp%kmsx = grid%sm33x
   grid%xp%kmex = grid%em33x

   grid%xp%itsx = grid%sp31x
   grid%xp%itex = grid%ep31x
   grid%xp%jtsx = grid%sp32x
   grid%xp%jtex = grid%ep32x
   grid%xp%ktsx = grid%sp33x
   grid%xp%ktex = grid%ep33x

   grid%xp%ipsx = grid%sp31x
   grid%xp%ipex = grid%ep31x
   grid%xp%jpsx = grid%sp32x
   grid%xp%jpex = grid%ep32x
   grid%xp%kpsx = grid%sp33x
   grid%xp%kpex = grid%ep33x

   ! Indices for xz decomposition

   grid%xp%idsy = grid%sd31
   grid%xp%idey = grid%ed31 - 1
   grid%xp%jdsy = grid%sd32
   grid%xp%jdey = grid%ed32 - 1
   grid%xp%kdsy = grid%sd33
   grid%xp%kdey = grid%ed33 - 1

   grid%xp%imsy = grid%sm31y
   grid%xp%imey = grid%em31y
   grid%xp%jmsy = grid%sm32y
   grid%xp%jmey = grid%em32y
   grid%xp%kmsy = grid%sm33y
   grid%xp%kmey = grid%em33y

   grid%xp%itsy = grid%sp31y
   grid%xp%itey = grid%ep31y
   grid%xp%jtsy = grid%sp32y
   grid%xp%jtey = grid%ep32y
   grid%xp%ktsy = grid%sp33y
   grid%xp%ktey = grid%ep33y

   grid%xp%ipsy = grid%sp31y
   grid%xp%ipey = grid%ep31y
   grid%xp%jpsy = grid%sp32y
   grid%xp%jpey = grid%ep32y
   grid%xp%kpsy = grid%sp33y
   grid%xp%kpey = grid%ep33y

   if (ipe > ide) ipe = ide
   if (jpe > jde) jpe = jde
   if (kpe > kde) kpe = kde

   ! Indices for yz decomposition

   if (grid%xp%itex > ide) grid%xp%itex = ide
   if (grid%xp%jtex > jde) grid%xp%jtex = jde
   if (grid%xp%ktex > kde) grid%xp%ktex = kde

   if (grid%xp%ipex > ide) grid%xp%ipex = ide
   if (grid%xp%jpex > jde) grid%xp%jpex = jde
   if (grid%xp%kpex > kde) grid%xp%kpex = kde

   ! Indices for xz decomposition

   if (grid%xp%itey > ide) grid%xp%itey = ide
   if (grid%xp%jtey > jde) grid%xp%jtey = jde
   if (grid%xp%ktey > kde) grid%xp%ktey = kde

   if (grid%xp%ipey > ide) grid%xp%ipey = ide
   if (grid%xp%jpey > jde) grid%xp%jpey = jde
   if (grid%xp%kpey > kde) grid%xp%kpey = kde

   ! Copy grid%xpose dimensions from grid structure to grid%xp structure.

   ! Indices for xy decomposition

   grid%xp%ids = ids
   grid%xp%ide = ide
   grid%xp%jds = jds
   grid%xp%jde = jde
   grid%xp%kds = kds
   grid%xp%kde = kde

   grid%xp%ims = ims
   grid%xp%ime = ime
   grid%xp%jms = jms
   grid%xp%jme = jme
   grid%xp%kms = kms
   grid%xp%kme = kme

   grid%xp%ips = ips
   grid%xp%ipe = ipe
   grid%xp%jps = jps
   grid%xp%jpe = jpe
   grid%xp%kps = kps
   grid%xp%kpe = kpe

   if (print_detail_parallel) then
      write(unit=stdout, fmt='(2(a, i4, 5x))') &
           'grid%xp%ids =', grid%xp%ids , 'grid%xp%ide =', grid%xp%ide , &
           'grid%xp%jds =', grid%xp%jds , 'grid%xp%jde =', grid%xp%jde , &
           'grid%xp%kds =', grid%xp%kds , 'grid%xp%kde =', grid%xp%kde
      write(unit=stdout, fmt='(//)')

      write(unit=stdout, fmt='(2(a, i4, 5x))') &
           'grid%xp%ims =', grid%xp%ims , 'grid%xp%ime =', grid%xp%ime , &
           'grid%xp%jms =', grid%xp%jms , 'grid%xp%jme =', grid%xp%jme , &
           'grid%xp%kms =', grid%xp%kms , 'grid%xp%kme =', grid%xp%kme
      write(unit=stdout, fmt='(//)')

      write(unit=stdout, fmt='(2(a, i4, 5x))') &
           'grid%xp%ips =', grid%xp%ips , 'grid%xp%ipe =', grid%xp%ipe , &
           'grid%xp%jps =', grid%xp%jps , 'grid%xp%jpe =', grid%xp%jpe , &
           'grid%xp%kps =', grid%xp%kps , 'grid%xp%kpe =', grid%xp%kpe
      write(unit=stdout, fmt='(//)')

      write(unit=stdout, fmt='(2(a, i4, 5x))') &
           'grid%xp%imsx =', grid%xp%imsx, 'grid%xp%imex=', grid%xp%imex, &
           'grid%xp%jmsx =', grid%xp%jmsx, 'grid%xp%jmex=', grid%xp%jmex, &
           'grid%xp%kmsx =', grid%xp%kmsx, 'grid%xp%kmex=', grid%xp%kmex
      write(unit=stdout, fmt='(//)')

      write(unit=stdout, fmt='(2(a, i4, 5x))') &
           'grid%xp%ipsx =', grid%xp%ipsx, 'grid%xp%ipex=', grid%xp%ipex, &
           'grid%xp%jpsx =', grid%xp%jpsx, 'grid%xp%jpex=', grid%xp%jpex, &
           'grid%xp%kpsx =', grid%xp%kpsx, 'grid%xp%kpex=', grid%xp%kpex
      write(unit=stdout, fmt='(//)')

      write(unit=stdout, fmt='(2(a, i4, 5x))') &
           'grid%xp%imsy =', grid%xp%imsy, 'grid%xp%imey=', grid%xp%imey, &
           'grid%xp%jmsy =', grid%xp%jmsy, 'grid%xp%jmey=', grid%xp%jmey, &
           'grid%xp%kmsy =', grid%xp%kmsy, 'grid%xp%kmey=', grid%xp%kmey
      write(unit=stdout, fmt='(//)')

      write(unit=stdout, fmt='(2(a, i4, 5x))') &
           'grid%xp%ipsy =', grid%xp%ipsy, 'grid%xp%ipey=', grid%xp%ipey, &
           'grid%xp%jpsy =', grid%xp%jpsy, 'grid%xp%jpey=', grid%xp%jpey, &
           'grid%xp%kpsy =', grid%xp%kpsy, 'grid%xp%kpey=', grid%xp%kpey
   end if

   if (trace_use_dull) call da_trace_exit("da_copy_dims")

end subroutine da_copy_dims


subroutine da_copy_tile_dims(grid)

   !---------------------------------------------------------------------------
   ! Purpose: Copy tile dimensions from grid structure.
   !---------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)         :: grid

   integer                  :: ij   ! Loop counter

   if (trace_use_dull) call da_trace_entry("da_copy_tile_dims")

   ! De-reference tile indices stored in the grid data structure.

   its = MINVAL( grid%i_start(1:grid%num_tiles) )
   ite = MAXVAL( grid%i_end(1:grid%num_tiles) )
   jts = MINVAL( grid%j_start(1:grid%num_tiles) )
   jte = MAXVAL( grid%j_end(1:grid%num_tiles) )
   kts = grid%xp%kds
   kte = grid%xp%kde

   grid%xp%its = its
   grid%xp%ite = ite
   grid%xp%jts = jts
   grid%xp%jte = jte
   grid%xp%kts = kts
   grid%xp%kte = kte

   if (grid%xp%ite > grid%xp%ide) grid%xp%ite = grid%xp%ide
   if (grid%xp%jte > grid%xp%jde) grid%xp%jte = grid%xp%jde
   if (grid%xp%kte > grid%xp%kde) grid%xp%kte = grid%xp%kde

   if (ite > grid%xp%ide) ite = grid%xp%ide
   if (jte > grid%xp%jde) jte = grid%xp%jde
   if (kte > grid%xp%kde) kte = grid%xp%kde

   do ij = 1 , grid%num_tiles
      if (print_detail_parallel) then
         write(unit=stdout, fmt='(/)')
         write(unit=stdout, fmt='(A,i3,A,5x,3(i3,A,i3,5x))') 'Tile ',ij, &
                 ' size:', its,':',ite, jts,':',jte, kts,':',kte
      end if
   end do

   if (trace_use_dull) call da_trace_exit("da_copy_tile_dims")

end subroutine da_copy_tile_dims


subroutine da_pack_count_obs (num_obs, offset, value)

   !---------------------------------------------------------------------------
   ! Purpose: Pack the 4 integer num_obs values into value(offset) to 
   !          value(offset+3).
   !---------------------------------------------------------------------------

   implicit none

   type(count_obs_number_type), intent(in)     :: num_obs
   integer,                     intent(inout)  :: offset
   integer,                     intent(inout)  :: value(*)

   if (trace_use_dull) call da_trace_entry("da_pack_count_obs")

   value(offset)   = num_obs % num_used
   value(offset+1) = num_obs % num_outside_iyjx
   value(offset+2) = num_obs % num_max_err_chk
   value(offset+3) = num_obs % num_missing

   offset = offset + 4

   if (trace_use_dull) call da_trace_exit("da_pack_count_obs")

end subroutine da_pack_count_obs


subroutine da_unpack_count_obs (num_obs, offset, value)

   !---------------------------------------------------------------------------
   ! Purpose: Unpack the 4 integer values starting at value(offset) into the
   !          num_obs structure.
   !---------------------------------------------------------------------------

   implicit none

   type (count_obs_number_type), intent(out)    :: num_obs
   integer                     , intent(inout)  :: offset
   integer                     , intent(in)     :: value(*)

   if (trace_use) call da_trace_entry("da_unpack_count_obs")

   num_obs % num_used         = value(offset)
   num_obs % num_outside_iyjx = value(offset+1)
   num_obs % num_max_err_chk  = value(offset+2)
   num_obs % num_missing      = value(offset+3)
   offset = offset + 4

   if (trace_use) call da_trace_exit("da_unpack_count_obs")

end subroutine da_unpack_count_obs


subroutine da_transpose_x2y (grid)

   implicit none

   type(domain), intent(inout)               :: grid
   integer                                   :: ij, i, j, k

   if (trace_use_dull) call da_trace_entry("da_transpose_x2y")

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/XPOSE_V1_x2y.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
  call trans_x2y ( ntasks_y, local_communicator_y, 1, 8, 4, DATA_ORDER_XYZ , &
                   grid%xp%v1x, &  ! variable in X decomp
                   grid%sd31, grid%ed31, grid%sd32, grid%ed32, grid%sd33, grid%ed33, &
                   grid%sp31x, grid%ep31x, grid%sp32x, grid%ep32x, grid%sp33x, grid%ep33x, &
                   grid%sm31x, grid%em31x, grid%sm32x, grid%em32x, grid%sm33x, grid%em33x, &
                   grid%xp%v1y, &  ! variable in Y decomp
                   grid%sp31y, grid%ep31y, grid%sp32y, grid%ep32y, grid%sp33y, grid%ep33y, &
                   grid%sm31y, grid%em31y, grid%sm32y, grid%em32y, grid%sm33y, grid%em33y ) 
!ENDOFREGISTRYGENERATEDINCLUDE

   if (trace_use_dull) call da_trace_exit("da_transpose_x2y")

end subroutine da_transpose_x2y


subroutine da_transpose_y2x (grid)

   implicit none

   type(domain), intent(inout)               :: grid
   integer                                   :: ij, i, j, k

   if (trace_use_dull) call da_trace_entry("da_transpose_y2x")

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/XPOSE_V1_y2x.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
  call trans_x2y ( ntasks_y, local_communicator_y, 0, 8, 4, DATA_ORDER_XYZ , &
                   grid%xp%v1x, &  ! variable in X decomp
                   grid%sd31, grid%ed31, grid%sd32, grid%ed32, grid%sd33, grid%ed33, &
                   grid%sp31x, grid%ep31x, grid%sp32x, grid%ep32x, grid%sp33x, grid%ep33x, &
                   grid%sm31x, grid%em31x, grid%sm32x, grid%em32x, grid%sm33x, grid%em33x, &
                   grid%xp%v1y, &  ! variable in Y decomp
                   grid%sp31y, grid%ep31y, grid%sp32y, grid%ep32y, grid%sp33y, grid%ep33y, &
                   grid%sm31y, grid%em31y, grid%sm32y, grid%em32y, grid%sm33y, grid%em33y ) 
!ENDOFREGISTRYGENERATEDINCLUDE

   if (trace_use_dull) call da_trace_exit("da_transpose_y2x")

end subroutine da_transpose_y2x


subroutine da_transpose_z2x (grid)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)               :: grid
   integer                                   :: ij, i, j, k

   if (trace_use_dull) call da_trace_entry("da_transpose_z2x")
   
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/XPOSE_V1_z2x.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
  call trans_z2x ( ntasks_x, local_communicator_x, 1, 8, 4, DATA_ORDER_XYZ , &
                   grid%xp%v1z, &  ! variable in Z decomp
                   grid%sd31, grid%ed31, grid%sd32, grid%ed32, grid%sd33, grid%ed33, &
                   grid%sp31, grid%ep31, grid%sp32, grid%ep32, grid%sp33, grid%ep33, &
                   grid%sm31, grid%em31, grid%sm32, grid%em32, grid%sm33, grid%em33, &
                   grid%xp%v1x, &  ! variable in X decomp
                   grid%sp31x, grid%ep31x, grid%sp32x, grid%ep32x, grid%sp33x, grid%ep33x, &
                   grid%sm31x, grid%em31x, grid%sm32x, grid%em32x, grid%sm33x, grid%em33x ) 
!ENDOFREGISTRYGENERATEDINCLUDE

   if (trace_use_dull) call da_trace_exit("da_transpose_z2x")

end subroutine da_transpose_z2x


subroutine da_transpose_x2z (grid)

   implicit none

   type(domain), intent(inout)               :: grid
   integer                                   :: ij, i, j, k

   if (trace_use_dull) call da_trace_entry("da_transpose_x2z")

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/XPOSE_V1_x2z.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
  call trans_z2x ( ntasks_x, local_communicator_x, 0, 8, 4, DATA_ORDER_XYZ , &
                   grid%xp%v1z, &  ! variable in Z decomp
                   grid%sd31, grid%ed31, grid%sd32, grid%ed32, grid%sd33, grid%ed33, &
                   grid%sp31, grid%ep31, grid%sp32, grid%ep32, grid%sp33, grid%ep33, &
                   grid%sm31, grid%em31, grid%sm32, grid%em32, grid%sm33, grid%em33, &
                   grid%xp%v1x, &  ! variable in X decomp
                   grid%sp31x, grid%ep31x, grid%sp32x, grid%ep32x, grid%sp33x, grid%ep33x, &
                   grid%sm31x, grid%em31x, grid%sm32x, grid%em32x, grid%sm33x, grid%em33x ) 
!ENDOFREGISTRYGENERATEDINCLUDE

   if (trace_use_dull) call da_trace_exit("da_transpose_x2z")

end subroutine da_transpose_x2z


subroutine da_transpose_y2z (grid)

   implicit none

   type(domain), intent(inout)               :: grid
   integer                                   :: ij, i, j, k

   if (trace_use_dull) call da_trace_entry("da_transpose_y2z")
   
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/XPOSE_V1_y2z.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
  call trans_x2y ( ntasks_y, local_communicator_y, 0, 8, 4, DATA_ORDER_XYZ , &
                   grid%xp%v1x, &  ! variable in X decomp
                   grid%sd31, grid%ed31, grid%sd32, grid%ed32, grid%sd33, grid%ed33, &
                   grid%sp31x, grid%ep31x, grid%sp32x, grid%ep32x, grid%sp33x, grid%ep33x, &
                   grid%sm31x, grid%em31x, grid%sm32x, grid%em32x, grid%sm33x, grid%em33x, &
                   grid%xp%v1y, &  ! variable in Y decomp
                   grid%sp31y, grid%ep31y, grid%sp32y, grid%ep32y, grid%sp33y, grid%ep33y, &
                   grid%sm31y, grid%em31y, grid%sm32y, grid%em32y, grid%sm33y, grid%em33y ) 
  call trans_z2x ( ntasks_x, local_communicator_x, 0, 8, 4, DATA_ORDER_XYZ , &
                   grid%xp%v1z, &  ! variable in Z decomp
                   grid%sd31, grid%ed31, grid%sd32, grid%ed32, grid%sd33, grid%ed33, &
                   grid%sp31, grid%ep31, grid%sp32, grid%ep32, grid%sp33, grid%ep33, &
                   grid%sm31, grid%em31, grid%sm32, grid%em32, grid%sm33, grid%em33, &
                   grid%xp%v1x, &  ! variable in X decomp
                   grid%sp31x, grid%ep31x, grid%sp32x, grid%ep32x, grid%sp33x, grid%ep33x, &
                   grid%sm31x, grid%em31x, grid%sm32x, grid%em32x, grid%sm33x, grid%em33x)
!ENDOFREGISTRYGENERATEDINCLUDE

   if (trace_use_dull) call da_trace_exit("da_transpose_y2z")

end subroutine da_transpose_y2z


subroutine da_transpose_z2y (grid)

   implicit none

   type(domain), intent(inout)               :: grid
   integer                                   :: ij, i, j, k

   if (trace_use_dull) call da_trace_entry("da_transpose_z2y")

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/XPOSE_V1_z2y.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
  call trans_z2x ( ntasks_x, local_communicator_x, 1, 8, 4, DATA_ORDER_XYZ , &
                   grid%xp%v1z, &  ! variable in Z decomp
                   grid%sd31, grid%ed31, grid%sd32, grid%ed32, grid%sd33, grid%ed33, &
                   grid%sp31, grid%ep31, grid%sp32, grid%ep32, grid%sp33, grid%ep33, &
                   grid%sm31, grid%em31, grid%sm32, grid%em32, grid%sm33, grid%em33, &
                   grid%xp%v1x, &  ! variable in X decomp
                   grid%sp31x, grid%ep31x, grid%sp32x, grid%ep32x, grid%sp33x, grid%ep33x, &
                   grid%sm31x, grid%em31x, grid%sm32x, grid%em32x, grid%sm33x, grid%em33x )
  call trans_x2y ( ntasks_y, local_communicator_y, 1, 8, 4, DATA_ORDER_XYZ , &
                   grid%xp%v1x, &  ! variable in X decomp
                   grid%sd31, grid%ed31, grid%sd32, grid%ed32, grid%sd33, grid%ed33, &
                   grid%sp31x, grid%ep31x, grid%sp32x, grid%ep32x, grid%sp33x, grid%ep33x, &
                   grid%sm31x, grid%em31x, grid%sm32x, grid%em32x, grid%sm33x, grid%em33x, &
                   grid%xp%v1y, &  ! variable in Y decomp
                   grid%sp31y, grid%ep31y, grid%sp32y, grid%ep32y, grid%sp33y, grid%ep33y, &
                   grid%sm31y, grid%em31y, grid%sm32y, grid%em32y, grid%sm33y, grid%em33y ) 
!ENDOFREGISTRYGENERATEDINCLUDE

   if (trace_use_dull) call da_trace_exit("da_transpose_z2y")

end subroutine da_transpose_z2y


subroutine da_transpose_x2y_v2 (grid)

   implicit none

   type(domain), intent(inout)               :: grid
   integer                                   :: ij, i, j, k

   if (trace_use) call da_trace_entry("da_transpose_x2y_v2")
   
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/XPOSE_V2_x2y.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
  call trans_x2y ( ntasks_y, local_communicator_y, 1, 8, 4, DATA_ORDER_XYZ , &
                   grid%xp%v2x, &  ! variable in X decomp
                   grid%sd31, grid%ed31, grid%sd32, grid%ed32, grid%sd33, grid%ed33, &
                   grid%sp31x, grid%ep31x, grid%sp32x, grid%ep32x, grid%sp33x, grid%ep33x, &
                   grid%sm31x, grid%em31x, grid%sm32x, grid%em32x, grid%sm33x, grid%em33x, &
                   grid%xp%v2y, &  ! variable in Y decomp
                   grid%sp31y, grid%ep31y, grid%sp32y, grid%ep32y, grid%sp33y, grid%ep33y, &
                   grid%sm31y, grid%em31y, grid%sm32y, grid%em32y, grid%sm33y, grid%em33y ) 
!ENDOFREGISTRYGENERATEDINCLUDE

   if (trace_use) call da_trace_exit("da_transpose_x2y_v2")

end subroutine da_transpose_x2y_v2


subroutine da_transpose_y2x_v2 (grid)

   implicit none

   type(domain), intent(inout)               :: grid
   integer                                   :: ij, i, j, k

   if (trace_use) call da_trace_entry("da_transpose_y2x_v2")
   
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/XPOSE_V2_y2x.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
  call trans_x2y ( ntasks_y, local_communicator_y, 0, 8, 4, DATA_ORDER_XYZ , &
                   grid%xp%v2x, &  ! variable in X decomp
                   grid%sd31, grid%ed31, grid%sd32, grid%ed32, grid%sd33, grid%ed33, &
                   grid%sp31x, grid%ep31x, grid%sp32x, grid%ep32x, grid%sp33x, grid%ep33x, &
                   grid%sm31x, grid%em31x, grid%sm32x, grid%em32x, grid%sm33x, grid%em33x, &
                   grid%xp%v2y, &  ! variable in Y decomp
                   grid%sp31y, grid%ep31y, grid%sp32y, grid%ep32y, grid%sp33y, grid%ep33y, &
                   grid%sm31y, grid%em31y, grid%sm32y, grid%em32y, grid%sm33y, grid%em33y ) 
!ENDOFREGISTRYGENERATEDINCLUDE

   if (trace_use) call da_trace_exit("da_transpose_y2x_v2")

end subroutine da_transpose_y2x_v2



subroutine da_cv_to_global(cv_size, cv_size_global, x, grid, mzs, xg )


   !-----------------------------------------------------------------------
   ! Purpose: Gathers local cv-array x into domain cv-array xg(:).  
   ! Global cv-array xg will only be valid on the "monitor" task.  
   !
   ! Must be called by all MPI tasks.  
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: cv_size        ! Size of local cv-array
   integer,          intent(in)    :: cv_size_global ! Size of domain cv-array
   real,             intent(in)    :: x(1:cv_size)   ! local cv-array
   type(domain),     intent(in)    :: grid
   integer,          intent(in)    :: mzs(:)  ! mz for each variable
                                         ! (to identify empty or 2D arrays)
   real,             intent(inout) :: xg(1:cv_size_global) ! global cv-array


   ! Local declarations
   type (vp_type) :: vv_x    ! Grdipt/EOF cv_array (local)
   type (vp_type) :: vv_xg   ! Grdipt/EOF cv_array (global)
   type (xpose_type) :: xpg  ! global dimensions
   integer :: ids, ide, jds, jde, kds, &
              ims, ime, jms, jme, kms, &
              ips, ipe, jps, jpe, kps 

   integer :: n

   if (trace_use) call da_trace_entry("da_cv_to_global")      

   !
   ! Gather to mimic serial summation order.  
   !

   ! k?e varies with variable v1 - v5
   ids = ids; ide = ide; jds = jds; jde = jde; kds = kds
   ims = ims; ime = ime; jms = jms; jme = jme; kms = grid%xp%kms
   ips = grid%xp%ips; ipe = grid%xp%ipe; jps = grid%xp%jps; jpe = grid%xp%jpe; kps = grid%xp%kps

   ! TOdo:  encapsulate this crap!  
   ! allocate vv_x
   allocate(vv_x%v1(ims:ime,jms:jme,mzs(1)))
   allocate(vv_x%v2(ims:ime,jms:jme,mzs(2)))
   allocate(vv_x%v3(ims:ime,jms:jme,mzs(3)))
   allocate(vv_x%v4(ims:ime,jms:jme,mzs(4)))
   allocate(vv_x%v5(ims:ime,jms:jme,mzs(5)))
   allocate(vv_x%alpha(ims_int:ime_int,jms_int:jme_int,kms_int:kme_int,mzs(7)))

   call da_cv_to_vv (cv_size, x, mzs, vv_x )

   if (rootproc) then
      ! allocate vv_xg
      allocate(vv_xg%v1(ids:ide,jds:jde,mzs(1)))
      allocate(vv_xg%v2(ids:ide,jds:jde,mzs(2)))
      allocate(vv_xg%v3(ids:ide,jds:jde,mzs(3)))
      allocate(vv_xg%v4(ids:ide,jds:jde,mzs(4)))
      allocate(vv_xg%v5(ids:ide,jds:jde,mzs(5)))
!     allocate(vv_xg%alpha(ids:ide,jds:jde,kds:kde,mzs(7)))
      allocate(vv_xg%alpha(ids_int:ide_int,jds_int:jde_int,kds_int:kde_int,mzs(7)))
   else
      ! Allocate dummy array for non-monitor process to keep Fortran
      ! CICO happy...
      allocate(vv_xg%v1(1,1,1))
      allocate(vv_xg%v2(1,1,1))
      allocate(vv_xg%v3(1,1,1))
      allocate(vv_xg%v4(1,1,1))
      allocate(vv_xg%v5(1,1,1))
      allocate(vv_xg%alpha(1,1,1,1))
   end if

   ! TOdo:  encapsulate this crap!  
   ! gather to global data structures
   ! it is possible to gather straight into a globally-sized cv-array
   ! TOdo:  add this optimization later
   call da_patch_to_global(grid, vv_x%v1,    vv_xg%v1,    mzs(1))
   call da_patch_to_global(grid, vv_x%v2,    vv_xg%v2,    mzs(2))
   call da_patch_to_global(grid, vv_x%v3,    vv_xg%v3,    mzs(3))
   call da_patch_to_global(grid, vv_x%v4,    vv_xg%v4,    mzs(4))
   call da_patch_to_global(grid, vv_x%v5,    vv_xg%v5,    mzs(5))
   if ( mzs(7) > 0 ) then
      do n = 1, mzs(7) ! Ensemble size
         if ( anal_type_hybrid_dual_res ) then
            call da_patch_to_global_dual_res(grid%intermediate_grid, vv_x%alpha(:,:,:,n), vv_xg%alpha(:,:,:,n), mzs(10))
         else
            call da_patch_to_global(grid, vv_x%alpha(:,:,:,n), vv_xg%alpha(:,:,:,n), mzs(6))
         endif
      end do
   end if
   ! deallocate vv_x
   deallocate (vv_x%v1, vv_x%v2, vv_x%v3, vv_x%v4, vv_x%v5, vv_x%alpha)
   if (rootproc) then
      ! finally, collapse data back into a globally-sized cv-array
      xpg%ids = ids; xpg%ide = ide
      xpg%ims = ids; xpg%ime = ide
      xpg%its = ids; xpg%ite = ide
      xpg%jds = jds; xpg%jde = jde
      xpg%jms = jds; xpg%jme = jde
      xpg%jts = jds; xpg%jte = jde
      xpg%kds = kds; xpg%kde = kde
      xpg%kms = kds; xpg%kme = kde
      xpg%kts = kds; xpg%kte = kde
      call da_vv_to_cv(vv_xg, xpg, mzs, cv_size_global, xg)
   end if
   ! deallocate vv_xg
   deallocate(vv_xg%v1, vv_xg%v2, vv_xg%v3, vv_xg%v4, vv_xg%v5, vv_xg%alpha)
   if (trace_use) call da_trace_exit("da_cv_to_global")    


end subroutine da_cv_to_global


subroutine da_patch_to_global_2d (grid, vlocal, vglobal)

   !---------------------------------------------------------------------
   ! Purpose: Gathers local 2D array vlocal into global array vglobal. 
   !
   ! Must be called by all MPI tasks.
   !---------------------------------------------------------------------  

   implicit none

   type(domain), intent(in)  :: grid
   real,         intent(in)  :: vlocal(:,:)
   real,         intent(out) :: vglobal(:,:)

   real, allocatable :: vlocal3d(:,:,:), vglobal3d(:,:,:)

   if (trace_use_frequent) call da_trace_entry("da_patch_to_global_2d")

   allocate(vlocal3d (ims:ime, jms:jme, 1:1))
   allocate(vglobal3d(ids:ide, jds:jde, 1:1))

   vlocal3d(:,:,1) = vlocal(:,:)
   call da_patch_to_global_3d(grid, vlocal3d, vglobal3d, 1)
   if (rootproc) then
      vglobal(:,:) = vglobal3d(:,:,1)
   end if

   deallocate(vlocal3d)
   deallocate(vglobal3d)

   if (trace_use_frequent) call da_trace_exit("da_patch_to_global_2d")

end subroutine da_patch_to_global_2d


subroutine da_patch_to_global_3d(grid, vlocal, vglobal, mz)

   !----------------------------------------------------------------------
   ! Purpose: Gathers local 3D array vlocal into global array vglobal.  
   ! Assumes that "k" is not decomposed.  End indices in the "k" dimension 
   ! are inferred from mz, which can be less than kde.  
   !
   ! Must be called by all MPI tasks.  
   !----------------------------------------------------------------------

   implicit none

   type(domain),      intent(in)  :: grid
   real,              intent(in)  :: vlocal(ims:ime,jms:jme,kms:kme)
   real,              intent(out) :: vglobal(ids:ide,jds:jde,kds:kde)
   integer, optional, intent(in)  :: mz


   integer :: local_kde
   integer :: local_kme
   integer :: local_kpe

   if (trace_use_frequent) call da_trace_entry("da_patch_to_global_3d")


   if (present(mz)) then
      local_kde = kds + mz - 1
      local_kme = local_kde
      local_kpe = local_kde
   else
      local_kde = kde
      local_kme = kme
      local_kpe = kpe
   end if

   if (local_kde > 0) then
      call wrf_patch_to_global_real (vlocal, vglobal, grid%xp%domdesc, &
         trim(grid_stagger), trim(grid_ordering), &
         ids, ide, jds, jde, kds, local_kde,  &
         ims, ime, jms, jme, kms, local_kme,  &
         ips, ipe, jps, jpe, kps, local_kpe)
   end if

   if (trace_use_frequent) call da_trace_exit("da_patch_to_global_3d")

end subroutine da_patch_to_global_3d


subroutine da_res_generic_set_info( res_generic, proc_domain, &
                                      obs_global_index)

   !---------------------------------------------------------------------------
   ! Purpose:  Eliminate redundant code for many obs types.  
   !
   ! Method:   Set common fields in res_generic.
   !---------------------------------------------------------------------------
   
   implicit none

   type(residual_generic_type), intent(inout) :: res_generic  ! generic type
   logical,                      intent(in  ) :: proc_domain
   integer,                      intent(in  ) :: obs_global_index

   res_generic%proc_domain = proc_domain
   res_generic%obs_global_index = obs_global_index
end subroutine da_res_generic_set_info


subroutine da_res_sound_create_template( template)

   !---------------------------------------------------------------------------
   ! Purpose:  Return storage template for specific type stored as 
   !           residual_generic_type.
   !---------------------------------------------------------------------------

   implicit none

   type(residual_template_type), intent(inout) :: template

   ! only vector data for this type
   template%lbnd = 1
   template%ubnd = 4

end subroutine da_res_sound_create_template


subroutine da_res_sound_to_generic( res_specific, iv_num_levels, &
                                      res_generic)

   !---------------------------------------------------------------------------
   ! Purpose:  Eliminate redundant code for many obs types.  
   !
   ! Method:   Specific type res_specific is translated to generic type 
   !           res_generic.  Pointer manipulation is used for vector data, no 
   !           vector data is copied.  Scalar data is copied.  This routine 
   !           allocates memory for res_generic.  The caller must ensure that 
   !           memory is deallocated later.  
   !           iv_num_levels is used as a sanity check and is ignored for 
   !           generic types that do not contain vector data. 
   !---------------------------------------------------------------------------

   implicit none

   type(residual_sound_type),   intent(in  ) :: res_specific ! specific type
   integer,                      intent(in  ) :: iv_num_levels ! levels
   type(residual_generic_type), intent(inout) :: res_generic  ! generic type
   ! Local declarations
   type(residual_template_type) :: template
   integer :: n

   call da_res_sound_create_template( template)
   allocate( res_generic%values(template%lbnd:template%ubnd))
   ! only vector data for this type
   ! store references to vector data
   res_generic%values(1)%ptr => res_specific%u
   res_generic%values(2)%ptr => res_specific%v
   res_generic%values(3)%ptr => res_specific%t
   res_generic%values(4)%ptr => res_specific%q
   ! ASSERTION
   do n=1,4
      !TBH:  NOTE:  We could handle iv_num_levels < size(res_generic%values(n)%ptr) 
      !TBH:         HOWEVER, we would have to add "num_levels" state to 
      !TBH:         residual_generic_type AND send this around.  Better to fix 
      !TBH:         allocation in specific types to avoid wasting memory!  
      if (size(res_generic%values(n)%ptr) /= iv_num_levels) then
         call da_error("da_generic_methods.inc",78, &
           (/'residual_sound_to_generic:  mismatched levels'/))
      end if
   end do

end subroutine da_res_sound_to_generic


subroutine da_res_sound_from_generic( res_generic, res_specific)

   !---------------------------------------------------------------------------
   ! Purpose:  Eliminate redundant code for many obs types.  
   !
   ! Method:   Generic type res_generic is translated to specific type 
   !           res_specific.  Pointer manipulation is used for vector data, no 
   !           vector data is copied.  Scalar data is copied. 
   !---------------------------------------------------------------------------

   implicit none

   type(residual_generic_type), intent(in  ) :: res_generic  ! generic type
   type(residual_sound_type),   intent(inout) :: res_specific ! specific type

   ! only vector data for this type
   ! store references to vector data
   res_specific%u => res_generic%values(1)%ptr
   res_specific%v => res_generic%values(2)%ptr
   res_specific%t => res_generic%values(3)%ptr
   res_specific%q => res_generic%values(4)%ptr

end subroutine da_res_sound_from_generic


subroutine da_res_synop_create_template( template)

   !---------------------------------------------------------------------------
   ! Purpose:  Return storage template for specific type stored as 
   !           residual_generic_type.  
   !---------------------------------------------------------------------------

   implicit none

   type(residual_template_type), intent(inout) :: template

   ! only scalar data for this type
   template%lbnd = 0
   template%ubnd = 0

end subroutine da_res_synop_create_template


subroutine da_res_synop_to_generic( res_specific, iv_num_levels, &
                                      res_generic)

   !---------------------------------------------------------------------------
   ! Purpose:  Eliminate redundant code for many obs types.  
   !
   ! Method:   Specific type res_specific is translated to generic type 
   !           res_generic.  Pointer manipulation is used for vector data, no 
   !           vector data is copied.  Scalar data is copied.  This routine 
   !           allocates memory for res_generic.  The caller must ensure that 
   !           memory is deallocated later.  
   !           iv_num_levels is used as a sanity check and is ignored for 
   !           generic types that do not contain vector data.
   !---------------------------------------------------------------------------

   implicit none

   type(residual_synop_type),   intent(in  ) :: res_specific ! specific type
   integer,                      intent(in  ) :: iv_num_levels ! levels
   type(residual_generic_type), intent(inout) :: res_generic  ! generic type
   ! Local declarations
   type(residual_template_type) :: template

   ! use to avoid compiler warning
   if (iv_num_levels==0) then
   end if

   call da_res_synop_create_template( template)
   allocate( res_generic%values(template%lbnd:template%ubnd))
   ! only scalar data for this type
   allocate( res_generic%values(0)%ptr(5))
   ! copy scalar data
   res_generic%values(0)%ptr(1) = res_specific%u
   res_generic%values(0)%ptr(2) = res_specific%v
   res_generic%values(0)%ptr(3) = res_specific%t
   res_generic%values(0)%ptr(4) = res_specific%p
   res_generic%values(0)%ptr(5) = res_specific%q

end subroutine da_res_synop_to_generic


subroutine da_res_synop_from_generic( res_generic, res_specific)

   !---------------------------------------------------------------------------
   ! Purpose:  Eliminate redundant code for many obs types.  
   !
   ! Method:   Generic type res_generic is translated to specific type 
   !           res_specific.  Pointer manipulation is used for vector data, no 
   !           vector data is copied.  Scalar data is copied.
   !---------------------------------------------------------------------------

   implicit none

   type(residual_generic_type), intent(in  ) :: res_generic  ! generic type
   type(residual_synop_type),   intent(inout) :: res_specific ! specific type

   ! only scalar data for this type
   ! copy scalar data
   res_specific%u = res_generic%values(0)%ptr(1)
   res_specific%v = res_generic%values(0)%ptr(2)
   res_specific%t = res_generic%values(0)%ptr(3)
   res_specific%p = res_generic%values(0)%ptr(4)
   res_specific%q = res_generic%values(0)%ptr(5)

end subroutine da_res_synop_from_generic


subroutine da_y_facade_create( slice, num_obs, num_obs_glo, template)

   !---------------------------------------------------------------------------
   ! Purpose:  Create a y_facade_type containing specified number of 
   !           residual_generic_type objects.  
   !
   ! Method:   Allocate memory and call residual_generic_create.  
   !           Call y_facade_free() to deallocate memory allocated here.  
   !           If template is not present, residual_generic_type objects will be 
   !           left uninitialized.  If template is present, then 
   !           residual_generic_type objects will be allocated but no values 
   !           will be copied into them.  
   !
   !---------------------------------------------------------------------------

   implicit none

   type(y_facade_type),          intent(inout) :: slice ! selected residual obs
   integer,                       intent(in  ) :: num_obs
   integer,                       intent(in  ) :: num_obs_glo
   type(residual_template_type), optional, intent(in  ) :: template
   ! Local declarations
   integer :: n

   slice%num_obs     = num_obs
   slice%num_obs_glo = num_obs_glo
   allocate( slice%obs(slice%num_obs))
   do n=1, slice%num_obs
      call da_res_generic_create( slice%obs(n), template=template)
   end do

end subroutine da_y_facade_create


subroutine da_res_generic_create( res_generic, template)

   !---------------------------------------------------------------------------
   ! Purpose:  Create a residual_generic_type object.  This must be called 
   !           before any other operation is performed on a 
   !           residual_generic_type object.  Do not pass an already-allocated 
   !           object into this routine as res_generic or you will cause a 
   !           memory leak!  
   !
   ! Method:   If template is not present, create an empty object.  
   !           If template is present, create an uninitialized object with 
   !             space allocated to match template.  Caller is responsible 
   !             for deallocation via call to residual_generic_free().  Note 
   !             that this is *NOT* a copy-constructor because values are not 
   !             copied.  Also, proc_domain and obs_global_index fields are 
   !             not copied from the template.  Finally, memory is not allocated 
   !             for vector or scalar data, these pointers 
   !            (res_generic%values(:)%ptr) are nullified.
   !---------------------------------------------------------------------------

   implicit none

   type(residual_generic_type), intent(inout) :: res_generic  ! generic type
   type(residual_template_type), optional, intent(in  ) :: template
   ! Local declarations
   integer :: i

   nullify( res_generic%values)
   if (present( template)) then
      allocate( res_generic%values(template%lbnd:template%ubnd))
      do i=template%lbnd, template%ubnd
         nullify( res_generic%values(i)%ptr)
      end do
   end if
end subroutine da_res_generic_create


subroutine da_res_generic_alloc_and_set( res_generic, val_index, values)

   ! TOdo:  replace this with a full-blown copy-constructor!  
   !---------------------------------------------------------------------------
   ! Purpose:  Allocates and initializes one of the values(either scalar or 
   !           vector) in an already-created residual_generic_type object.  
   !
   ! Method:   Allocate pointer and copy data from values.  Note that a call
   !           to residual_generic_free() will deallocate scalar data but leave
   !           vector data alone.  It is still the callers responsibility to
   !           deallocate pointers to vector data.  In this particular case, 
   !           any vector data allocated here is later passed by reference to 
   !           a global specific object via call to residual_*_from_generic() 
   !          (called from y_type_insert_sound_global()) and later explictly 
   !           deallocated via call to free_global_*().
   !---------------------------------------------------------------------------

   implicit none

   type(residual_generic_type), intent(inout) :: res_generic  ! generic type
   integer,                      intent(in  ) :: val_index    ! which value
   real,                         intent(in  ) :: values(:)    ! values

   if (( val_index < LBOUND(res_generic%values,1)) .OR. &
       ( val_index > UBOUND(res_generic%values,1))) then
      call da_error("da_generic_methods.inc",292, &
        (/'residual_generic_alloc_and_set:  bad val_index'/))
   end if
   allocate( res_generic%values(val_index)%ptr(size(values)))
   res_generic%values(val_index)%ptr(:) = values(:)

end subroutine da_res_generic_alloc_and_set


logical function da_res_generic_has_vector( res_generic)

   !---------------------------------------------------------------------------
   ! Purpose:  Returns .true. iff res_generic stores vector values
   !---------------------------------------------------------------------------

   implicit none

   type(residual_generic_type), intent(in) :: res_generic  ! generic type
   da_res_generic_has_vector =( UBOUND(res_generic%values,1) > 0)

end function da_res_generic_has_vector


logical function da_res_generic_has_scalar( res_generic)

   !---------------------------------------------------------------------------
   ! Purpose:  Returns .true. iff res_generic stores scalar values. 
   !---------------------------------------------------------------------------

   implicit none

   type(residual_generic_type), intent(in) :: res_generic  ! generic type
   da_res_generic_has_scalar =( LBOUND(res_generic%values,1) == 0)

end function da_res_generic_has_scalar


subroutine da_res_generic_free( res_generic)

   !---------------------------------------------------------------------------
   ! Purpose:  Free memory for residual_generic_type object.  
   !
   ! Method: pointers to vector values are assumed to point to memory allocated 
   !         by others and are not deallocated.  
   !         This will fail if called for a residual_generic_type object that 
   !         has never been created.  
   !---------------------------------------------------------------------------

   implicit none

   type(residual_generic_type), intent(inout) :: res_generic  ! generic type
   ! Local declarations
   logical :: oktofreevalues

   oktofreevalues = .false.
   if (da_res_generic_has_scalar( res_generic)) then
      ! free scalar memory
      deallocate( res_generic%values(0)%ptr)
      oktofreevalues = .true.
   end if
   if (da_res_generic_has_vector( res_generic)) then
      oktofreevalues = .true.
   end if
   if (oktofreevalues) then
      deallocate( res_generic%values)
   end if

end subroutine da_res_generic_free


subroutine da_y_facade_free( slice)

   !---------------------------------------------------------------------------
   ! Purpose:  Free memory for y_facade_type object.  
   !
   ! Method:   Brute force.  May want to make this smarter... 
   !---------------------------------------------------------------------------

   implicit none

   type(y_facade_type), intent(inout) :: slice
   ! Local declarations
   integer :: n
   if (associated( slice%obs)) then
      do n=1, size(slice%obs)
         call da_res_generic_free( slice%obs(n))
      end do
      deallocate( slice%obs)
   end if
end subroutine da_y_facade_free


subroutine da_deallocate_global_sonde_sfc(iv_glob, re_glob, jo_grad_y_glob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout) :: iv_glob        ! Innovation vector
   type (y_type),  intent(inout) :: re_glob        ! residual vector
   type (y_type),  intent(inout) :: jo_grad_y_glob ! Grad_y(Jo)

   if (trace_use_dull) call da_trace_entry("da_deallocate_global_sonde_sfc")

   deallocate(iv_glob%sonde_sfc)
   deallocate(re_glob%sonde_sfc)
   deallocate(jo_grad_y_glob%sonde_sfc)

   if (trace_use_dull) call da_trace_exit("da_deallocate_global_sonde_sfc")
   
end subroutine da_deallocate_global_sonde_sfc


subroutine da_deallocate_global_sound (iv_glob, re_glob, jo_grad_y_glob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv_glob        ! Innovation vector
   type(y_type),  intent(inout) :: re_glob        ! residual vector
   type(y_type),  intent(inout) :: jo_grad_y_glob ! Grad_y(Jo)

   integer :: n

   if (trace_use_dull) call da_trace_entry("da_deallocate_global_sound")

   deallocate(iv_glob%sound)
   do n=1,size(re_glob%sound)
      deallocate (re_glob%sound(n)%u)
      deallocate (re_glob%sound(n)%v)
      deallocate (re_glob%sound(n)%t)
      deallocate (re_glob%sound(n)%q)
   end do
   deallocate(re_glob%sound)
   do n=1,size(jo_grad_y_glob%sound)
      deallocate (jo_grad_y_glob%sound(n)%u)
      deallocate (jo_grad_y_glob%sound(n)%v)
      deallocate (jo_grad_y_glob%sound(n)%t)
      deallocate (jo_grad_y_glob%sound(n)%q)
   end do
   deallocate(jo_grad_y_glob%sound)

   if (trace_use_dull) call da_trace_exit("da_deallocate_global_sound")

end subroutine da_deallocate_global_sound


subroutine da_deallocate_global_synop(iv_glob, re_glob, jo_grad_y_glob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout) :: iv_glob        ! Innovation vector
   type (y_type),  intent(inout) :: re_glob        ! residual vector
   type (y_type),  intent(inout) :: jo_grad_y_glob ! Grad_y(Jo)

   if (trace_use_dull) call da_trace_entry("da_deallocate_global_synop")

   deallocate(iv_glob%synop)
   deallocate(re_glob%synop)
   deallocate(jo_grad_y_glob%synop)

   if (trace_use_dull) call da_trace_exit("da_deallocate_global_synop")

end subroutine da_deallocate_global_synop


!
! WRFVAR generic type macro file
!
! This file is used to generate a series of simple boiler-plate calls 
! to support residual generic types for bitwise-exact testing.  
! It contains M4 macros and then
! a series of invocations of the macros to generate the subroutine
! definitions, which are then included in another source file.  
!

! $1 = specific ob name, $2 = specific ob type, $3 = ob counter



















!--- sound sound num_sound

SUBROUTINE da_y_type_ex_sound( iv, re, slice )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Extract all sound obs from y and place them in generic 
!           object slice.  
!           Call da_y_facade_free() to deallocate memory allocated here.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (iv_type),       INTENT(IN   ) :: iv     ! Innovation vector
   type (y_type),        INTENT(IN   ) :: re     ! all residual obs
   type (y_facade_type), INTENT(INOUT) :: slice  ! selected residual obs
   ! Local declarations
   INTEGER :: n

   CALL da_y_facade_create( slice, iv%info(sound)%nlocal, iv%info(sound)%ntotal )
   DO n=1, slice%num_obs
stop
!     CALL da_res_generic_set_info( slice%obs(n),                     &
!                                     iv%sound(n)%loc%proc_domain,      &
!                                     iv%sound(n)%loc%obs_global_index )
!     CALL da_res_sound_to_generic( re%sound(n), iv%sound(n)%info%levels, &
!                                     slice%obs(n) )
   ENDDO

END SUBROUTINE da_y_type_ex_sound  


!--- sound sound num_sound

SUBROUTINE da_y_type_ins_sound_global( slice_glob, re_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert obs from generic object slice_glob into 
!           globally-scoped y_type re_glob.  re_glob is 
!           allocated minimally.  Caller must deallocate.  
!           Memory for global object slice_glob is deallocated here.  
!           Do not use slice_glob after this call.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(INOUT) :: slice_glob ! generic
   type (y_type),        INTENT(INOUT) :: re_glob    ! selected residual obs
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize obs
   ! deallocation is done in free_global_sound()
   ALLOCATE( re_glob%sound(slice_glob%num_obs) )
   DO n=1, slice_glob%num_obs
     CALL da_res_sound_from_generic( slice_glob%obs(n), re_glob%sound(n) )
   ENDDO
   re_glob%nlocal(sound) = slice_glob%num_obs  ! duplication!
   CALL da_y_facade_free( slice_glob )

END SUBROUTINE da_y_type_ins_sound_global 


!--- sound sound num_sound

SUBROUTINE da_iv_type_ins_sound_global( slice_glob, iv_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert meta-data from generic object slice_glob into 
!           globally-scoped iv_type iv_glob.  iv_glob is 
!           allocated minimally.  Caller must deallocate.  
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(IN   ) :: slice_glob ! selected residual obs
   type (iv_type),       INTENT(INOUT) :: iv_glob    ! partial Innovation vector
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize needed bits of iv_glob (ugly)
   iv_glob%info(sound)%nlocal  = slice_glob%num_obs
   iv_glob%info(sound)%ntotal = slice_glob%num_obs_glo
   ! deallocation is done in free_global_sound()
   ALLOCATE( iv_glob%sound(iv_glob%info(sound)%nlocal) )
   DO n=1, iv_glob%info(sound)%nlocal
stop
!     iv_glob%sound(n)%loc%proc_domain = slice_glob%obs(n)%proc_domain
!     iv_glob%sound(n)%loc%obs_global_index = &
!                                        slice_glob%obs(n)%obs_global_index
!     IF ( da_res_generic_has_vector( slice_glob%obs(n) ) ) THEN
!       iv_glob%sound(n)%info%levels = SIZE(slice_glob%obs(n)%values(1)%ptr)
!     ENDIF
   ENDDO

END SUBROUTINE da_iv_type_ins_sound_global  


!--- sound sound num_sound

!------------------------------------------------------------------------------
! PURPOSE:  Collect local arrays of residual_sound_type objects into 
!           global arrays in serial-code storage order.  This is used to 
!           perform bitwise-exact floating-point summations in 
!           serial-code-order during bitwise-exact testing of 
!           distributed-memory parallel configurations.  
!
! METHOD:   Indices stowed away during Read_Obs() are used to restore serial 
!           storage order.  Memory for global objects is allocated here.  
!           Global objects are minimally allocated to save memory.  
!           Memory is deallocated in free_global_sound().  
!------------------------------------------------------------------------------
  SUBROUTINE da_to_global_sound( iv,      re,      jo_grad_y, &
                                  iv_glob, re_glob, jo_grad_y_glob )

    IMPLICIT NONE

    ! task-local objects
    type (iv_type), INTENT( IN) :: iv             ! Innovation vector
    type (y_type),  INTENT( IN) :: re             ! residual vector
    type (y_type),  INTENT( IN) :: jo_grad_y      ! Grad_y(Jo)
    ! task-global objects
    type (iv_type), INTENT(OUT) :: iv_glob        ! Innovation vector
    type (y_type),  INTENT(OUT) :: re_glob        ! residual vector
    type (y_type),  INTENT(OUT) :: jo_grad_y_glob ! Grad_y(Jo)

    ! Local declarations
    type (y_facade_type) :: re_slice, re_glob_slice
    type (y_facade_type) :: jo_grad_y_slice, jo_grad_y_glob_slice
    type (residual_template_type) :: template  ! allocation info

    ! create process-local generic objects from specific objects
    CALL da_y_type_ex_sound( iv, re,        re_slice )
    CALL da_y_type_ex_sound( iv, jo_grad_y, jo_grad_y_slice )

    ! create global versions of generic objects from process-local objects
    ! and destroy process-local generic objects
    CALL da_res_sound_create_template( template )  ! use template in case 
                                                     ! some tasks have no obs
    CALL da_y_facade_to_global( re_slice,        template, re_glob_slice )
    CALL da_y_facade_to_global( jo_grad_y_slice, template, jo_grad_y_glob_slice )

    ! create global versions of specific objects
    ! and destroy global versions of generic objects
    ! iv first
    CALL da_iv_type_ins_sound_global( re_glob_slice, iv_glob )
    ! then y_types
    CALL da_y_type_ins_sound_global( re_glob_slice,        re_glob )
    CALL da_y_type_ins_sound_global( jo_grad_y_glob_slice, jo_grad_y_glob )
    ! global versions of specific objects are destroyed in 
    ! free_global_sound()

    RETURN

  END SUBROUTINE da_to_global_sound  


!--- sonde_sfc synop num_sound

SUBROUTINE da_y_type_ex_sonde_sfc( iv, re, slice )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Extract all sonde_sfc obs from y and place them in generic 
!           object slice.  
!           Call da_y_facade_free() to deallocate memory allocated here.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (iv_type),       INTENT(IN   ) :: iv     ! Innovation vector
   type (y_type),        INTENT(IN   ) :: re     ! all residual obs
   type (y_facade_type), INTENT(INOUT) :: slice  ! selected residual obs
   ! Local declarations
   INTEGER :: n

   CALL da_y_facade_create( slice, iv%info(sonde_sfc)%nlocal, iv%info(sonde_sfc)%ntotal )
   DO n=1, slice%num_obs
stop
!     CALL da_res_generic_set_info( slice%obs(n),                     &
!                                     iv%sonde_sfc(n)%loc%proc_domain,      &
!                                     iv%sonde_sfc(n)%loc%obs_global_index )
!     CALL da_res_synop_to_generic( re%sonde_sfc(n), iv%sonde_sfc(n)%info%levels, &
!                                     slice%obs(n) )
   ENDDO

END SUBROUTINE da_y_type_ex_sonde_sfc  


!--- sonde_sfc synop num_sound

SUBROUTINE da_y_type_ins_sonde_sfc_global( slice_glob, re_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert obs from generic object slice_glob into 
!           globally-scoped y_type re_glob.  re_glob is 
!           allocated minimally.  Caller must deallocate.  
!           Memory for global object slice_glob is deallocated here.  
!           Do not use slice_glob after this call.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(INOUT) :: slice_glob ! generic
   type (y_type),        INTENT(INOUT) :: re_glob    ! selected residual obs
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize obs
   ! deallocation is done in free_global_sonde_sfc()
   ALLOCATE( re_glob%sonde_sfc(slice_glob%num_obs) )
   DO n=1, slice_glob%num_obs
     CALL da_res_synop_from_generic( slice_glob%obs(n), re_glob%sonde_sfc(n) )
   ENDDO
   re_glob%nlocal(sonde_sfc) = slice_glob%num_obs  ! duplication!
   CALL da_y_facade_free( slice_glob )

END SUBROUTINE da_y_type_ins_sonde_sfc_global 


!--- sonde_sfc synop num_sound

SUBROUTINE da_iv_type_ins_sonde_sfc_global( slice_glob, iv_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert meta-data from generic object slice_glob into 
!           globally-scoped iv_type iv_glob.  iv_glob is 
!           allocated minimally.  Caller must deallocate.  
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(IN   ) :: slice_glob ! selected residual obs
   type (iv_type),       INTENT(INOUT) :: iv_glob    ! partial Innovation vector
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize needed bits of iv_glob (ugly)
   iv_glob%info(sonde_sfc)%nlocal  = slice_glob%num_obs
   iv_glob%info(sonde_sfc)%ntotal = slice_glob%num_obs_glo
   ! deallocation is done in free_global_sonde_sfc()
   ALLOCATE( iv_glob%sonde_sfc(iv_glob%info(sonde_sfc)%nlocal) )
   DO n=1, iv_glob%info(sonde_sfc)%nlocal
stop
!     iv_glob%sonde_sfc(n)%loc%proc_domain = slice_glob%obs(n)%proc_domain
!     iv_glob%sonde_sfc(n)%loc%obs_global_index = &
!                                        slice_glob%obs(n)%obs_global_index
!     IF ( da_res_generic_has_vector( slice_glob%obs(n) ) ) THEN
!       iv_glob%sonde_sfc(n)%info%levels = SIZE(slice_glob%obs(n)%values(1)%ptr)
!     ENDIF
   ENDDO

END SUBROUTINE da_iv_type_ins_sonde_sfc_global  


!--- sonde_sfc synop num_sound

!------------------------------------------------------------------------------
! PURPOSE:  Collect local arrays of residual_synop_type objects into 
!           global arrays in serial-code storage order.  This is used to 
!           perform bitwise-exact floating-point summations in 
!           serial-code-order during bitwise-exact testing of 
!           distributed-memory parallel configurations.  
!
! METHOD:   Indices stowed away during Read_Obs() are used to restore serial 
!           storage order.  Memory for global objects is allocated here.  
!           Global objects are minimally allocated to save memory.  
!           Memory is deallocated in free_global_sonde_sfc().  
!------------------------------------------------------------------------------
  SUBROUTINE da_to_global_sonde_sfc( iv,      re,      jo_grad_y, &
                                  iv_glob, re_glob, jo_grad_y_glob )

    IMPLICIT NONE

    ! task-local objects
    type (iv_type), INTENT( IN) :: iv             ! Innovation vector
    type (y_type),  INTENT( IN) :: re             ! residual vector
    type (y_type),  INTENT( IN) :: jo_grad_y      ! Grad_y(Jo)
    ! task-global objects
    type (iv_type), INTENT(OUT) :: iv_glob        ! Innovation vector
    type (y_type),  INTENT(OUT) :: re_glob        ! residual vector
    type (y_type),  INTENT(OUT) :: jo_grad_y_glob ! Grad_y(Jo)

    ! Local declarations
    type (y_facade_type) :: re_slice, re_glob_slice
    type (y_facade_type) :: jo_grad_y_slice, jo_grad_y_glob_slice
    type (residual_template_type) :: template  ! allocation info

    ! create process-local generic objects from specific objects
    CALL da_y_type_ex_sonde_sfc( iv, re,        re_slice )
    CALL da_y_type_ex_sonde_sfc( iv, jo_grad_y, jo_grad_y_slice )

    ! create global versions of generic objects from process-local objects
    ! and destroy process-local generic objects
    CALL da_res_synop_create_template( template )  ! use template in case 
                                                     ! some tasks have no obs
    CALL da_y_facade_to_global( re_slice,        template, re_glob_slice )
    CALL da_y_facade_to_global( jo_grad_y_slice, template, jo_grad_y_glob_slice )

    ! create global versions of specific objects
    ! and destroy global versions of generic objects
    ! iv first
    CALL da_iv_type_ins_sonde_sfc_global( re_glob_slice, iv_glob )
    ! then y_types
    CALL da_y_type_ins_sonde_sfc_global( re_glob_slice,        re_glob )
    CALL da_y_type_ins_sonde_sfc_global( jo_grad_y_glob_slice, jo_grad_y_glob )
    ! global versions of specific objects are destroyed in 
    ! free_global_sonde_sfc()

    RETURN

  END SUBROUTINE da_to_global_sonde_sfc  


!--- synop synop num_synop

SUBROUTINE da_y_type_ex_synop( iv, re, slice )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Extract all synop obs from y and place them in generic 
!           object slice.  
!           Call da_y_facade_free() to deallocate memory allocated here.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (iv_type),       INTENT(IN   ) :: iv     ! Innovation vector
   type (y_type),        INTENT(IN   ) :: re     ! all residual obs
   type (y_facade_type), INTENT(INOUT) :: slice  ! selected residual obs
   ! Local declarations
   INTEGER :: n

   CALL da_y_facade_create( slice, iv%info(synop)%nlocal, iv%info(synop)%ntotal )
   DO n=1, slice%num_obs
stop
!     CALL da_res_generic_set_info( slice%obs(n),                     &
!                                     iv%synop(n)%loc%proc_domain,      &
!                                     iv%synop(n)%loc%obs_global_index )
!     CALL da_res_synop_to_generic( re%synop(n), iv%synop(n)%info%levels, &
!                                     slice%obs(n) )
   ENDDO

END SUBROUTINE da_y_type_ex_synop  


!--- synop synop num_synop

SUBROUTINE da_y_type_ins_synop_global( slice_glob, re_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert obs from generic object slice_glob into 
!           globally-scoped y_type re_glob.  re_glob is 
!           allocated minimally.  Caller must deallocate.  
!           Memory for global object slice_glob is deallocated here.  
!           Do not use slice_glob after this call.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(INOUT) :: slice_glob ! generic
   type (y_type),        INTENT(INOUT) :: re_glob    ! selected residual obs
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize obs
   ! deallocation is done in free_global_synop()
   ALLOCATE( re_glob%synop(slice_glob%num_obs) )
   DO n=1, slice_glob%num_obs
     CALL da_res_synop_from_generic( slice_glob%obs(n), re_glob%synop(n) )
   ENDDO
   re_glob%nlocal(synop) = slice_glob%num_obs  ! duplication!
   CALL da_y_facade_free( slice_glob )

END SUBROUTINE da_y_type_ins_synop_global 


!--- synop synop num_synop

SUBROUTINE da_iv_type_ins_synop_global( slice_glob, iv_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert meta-data from generic object slice_glob into 
!           globally-scoped iv_type iv_glob.  iv_glob is 
!           allocated minimally.  Caller must deallocate.  
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(IN   ) :: slice_glob ! selected residual obs
   type (iv_type),       INTENT(INOUT) :: iv_glob    ! partial Innovation vector
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize needed bits of iv_glob (ugly)
   iv_glob%info(synop)%nlocal  = slice_glob%num_obs
   iv_glob%info(synop)%ntotal = slice_glob%num_obs_glo
   ! deallocation is done in free_global_synop()
   ALLOCATE( iv_glob%synop(iv_glob%info(synop)%nlocal) )
   DO n=1, iv_glob%info(synop)%nlocal
stop
!     iv_glob%synop(n)%loc%proc_domain = slice_glob%obs(n)%proc_domain
!     iv_glob%synop(n)%loc%obs_global_index = &
!                                        slice_glob%obs(n)%obs_global_index
!     IF ( da_res_generic_has_vector( slice_glob%obs(n) ) ) THEN
!       iv_glob%synop(n)%info%levels = SIZE(slice_glob%obs(n)%values(1)%ptr)
!     ENDIF
   ENDDO

END SUBROUTINE da_iv_type_ins_synop_global  


!--- synop synop num_synop

!------------------------------------------------------------------------------
! PURPOSE:  Collect local arrays of residual_synop_type objects into 
!           global arrays in serial-code storage order.  This is used to 
!           perform bitwise-exact floating-point summations in 
!           serial-code-order during bitwise-exact testing of 
!           distributed-memory parallel configurations.  
!
! METHOD:   Indices stowed away during Read_Obs() are used to restore serial 
!           storage order.  Memory for global objects is allocated here.  
!           Global objects are minimally allocated to save memory.  
!           Memory is deallocated in free_global_synop().  
!------------------------------------------------------------------------------
  SUBROUTINE da_to_global_synop( iv,      re,      jo_grad_y, &
                                  iv_glob, re_glob, jo_grad_y_glob )

    IMPLICIT NONE

    ! task-local objects
    type (iv_type), INTENT( IN) :: iv             ! Innovation vector
    type (y_type),  INTENT( IN) :: re             ! residual vector
    type (y_type),  INTENT( IN) :: jo_grad_y      ! Grad_y(Jo)
    ! task-global objects
    type (iv_type), INTENT(OUT) :: iv_glob        ! Innovation vector
    type (y_type),  INTENT(OUT) :: re_glob        ! residual vector
    type (y_type),  INTENT(OUT) :: jo_grad_y_glob ! Grad_y(Jo)

    ! Local declarations
    type (y_facade_type) :: re_slice, re_glob_slice
    type (y_facade_type) :: jo_grad_y_slice, jo_grad_y_glob_slice
    type (residual_template_type) :: template  ! allocation info

    ! create process-local generic objects from specific objects
    CALL da_y_type_ex_synop( iv, re,        re_slice )
    CALL da_y_type_ex_synop( iv, jo_grad_y, jo_grad_y_slice )

    ! create global versions of generic objects from process-local objects
    ! and destroy process-local generic objects
    CALL da_res_synop_create_template( template )  ! use template in case 
                                                     ! some tasks have no obs
    CALL da_y_facade_to_global( re_slice,        template, re_glob_slice )
    CALL da_y_facade_to_global( jo_grad_y_slice, template, jo_grad_y_glob_slice )

    ! create global versions of specific objects
    ! and destroy global versions of generic objects
    ! iv first
    CALL da_iv_type_ins_synop_global( re_glob_slice, iv_glob )
    ! then y_types
    CALL da_y_type_ins_synop_global( re_glob_slice,        re_glob )
    CALL da_y_type_ins_synop_global( jo_grad_y_glob_slice, jo_grad_y_glob )
    ! global versions of specific objects are destroyed in 
    ! free_global_synop()

    RETURN

  END SUBROUTINE da_to_global_synop  

subroutine da_y_facade_to_global( re_slice, template, re_glob_slice )

   !---------------------------------------------------------------------------
   ! Purpose:  Collect a local y_facade_type object into a global y_facade_type 
   !           object while preserving serial-code storage order.  This is 
   !           used to perform bitwise-exact floating-point summations in 
   !           serial-code-order during bitwise-exact testing of 
   !           distributed-memory parallel configurations.  
   !
   ! Method:   Indices stowed away during Read_Obs() are used to restore serial 
   !           storage order.  Memory for global objects is allocated here.  
   !           Global objects are minimally allocated to save memory.  
   !           Memory must be deallocated by the caller via a call to 
   !           da_y_facade_free().  
   !           Memory for local object re_slice is deallocated here.  Do not 
   !           use re_slice after this call.  
   !           The "template" argument is needed because some tasks may not 
   !           have any local obs.  
   !
   ! PERFORMANCE NOTE:   This implementation is NOT very efficient.  Speed it 
   !                     up if testing is too slow.  
   !---------------------------------------------------------------------------

   implicit none

   ! task-local objects  (really intent(in   ))
   type (y_facade_type),          intent(inout) :: re_slice      ! residual vector
   type (residual_template_type), intent(in)    :: template  ! allocation info
   ! task-global objects (really intent(  out))
   type (y_facade_type),          intent(inout) :: re_glob_slice ! residual vector

   ! Local declarations
   integer                      :: n, k, serial_n, serial_numvals
   integer                      :: proc
   integer                      :: num_obs_send
   integer                      :: buf_idx
   integer                      :: num_obs_all
   integer                      :: num_recv_all
   integer                      :: obs_global_index(re_slice%num_obs)
   integer                      :: num_values(re_slice%num_obs)
   integer                      :: num_send_buf_lcl
   integer,allocatable          :: num_obs_send_all(:)
   integer,allocatable          :: obs_global_index_all(:)
   integer,allocatable          :: obs_global_index_inv(:)
   integer,allocatable          :: obs_start_all(:)  ! start index of each obs
   integer,pointer              :: num_values_all(:)
   integer,allocatable          :: num_send_buf_all(:)
   integer,allocatable          :: recv_counts(:)
   integer,allocatable          :: recv_displs(:)
   real,allocatable             :: re_vals_lcl(:)
   real,allocatable             :: re_vals_all(:)
   integer                      :: i

   if (trace_use) call da_trace_entry("da_y_facade_to_global")

   ! todo:  later, upgrade from "allgather" to "gather-broadcast"

   ! collect information about all observations
   num_obs_send = 0
   obs_global_index = -1
   do n=1, re_slice%num_obs
      if (re_slice%obs(n)%proc_domain) then
         num_obs_send = num_obs_send + 1
         obs_global_index(num_obs_send) = re_slice%obs(n)%obs_global_index
      end if
   end do
   do n=1, num_obs_send
      if (obs_global_index(n) < 0) then
         call da_error("da_y_facade_to_global.inc",70, &
            (/'ASSERTION ERROR:  bad local value of obs_global_index'/))
      end if
   end do

   ! exchange num_obs_send and obs_global_index
   allocate (num_obs_send_all(0:num_procs-1))
   allocate (num_send_buf_all(0:num_procs-1))
   allocate (recv_counts(0:num_procs-1))
   allocate (recv_displs(0:num_procs-1))

   ! gather num_obs_send
   call mpi_allgather( num_obs_send, 1, mpi_integer, num_obs_send_all, 1, &
      mpi_integer, comm, ierr )
   num_obs_all = sum( num_obs_send_all )
   if ( num_obs_all /= re_slice%num_obs_glo ) then
      call da_error ("da_y_facade_to_global.inc",86, &
         (/'ASSERTION ERROR:  inconsistent count of sound obs'/))
   end if
   ! set up to gather obs_global_index
   recv_counts = num_obs_send_all
   recv_displs(0) = 0
   do proc=1, num_procs-1
      recv_displs(proc) = recv_displs(proc-1) + recv_counts(proc-1)
   end do
   allocate (num_values_all(num_obs_all))
   allocate (obs_global_index_all(num_obs_all))
   allocate (obs_global_index_inv(num_obs_all))
   allocate (obs_start_all(num_obs_all))

   ! gather obs_global_index
   call mpi_allgatherv( obs_global_index, num_obs_send, mpi_integer,    &
      obs_global_index_all, recv_counts, recv_displs, &
      mpi_integer, comm, ierr )

   ! invert obs_global_index_all
   ! initialize to "invalid" value
   obs_global_index_inv = -1
   do n=1, num_obs_all
      if ( ( obs_global_index_all(n) <  1 ) .OR. &
           ( obs_global_index_all(n) > size(obs_global_index_inv) ) ) then
         call da_error ("da_y_facade_to_global.inc",111, &
            (/'ASSERTION ERROR:  obs_global_index_all(n) out of range'/))
      end if
      if ( obs_global_index_inv(obs_global_index_all(n)) /= -1 ) then
         call da_error ("da_y_facade_to_global.inc",115, &
            (/'ASSERTION ERROR:  obs owned by more than one task'/))
      end if
      obs_global_index_inv(obs_global_index_all(n)) = n
   end do
   do n=1, num_obs_all
      if ( obs_global_index_inv(obs_global_index_all(n)) == -1 ) then
         call da_error ("da_y_facade_to_global.inc",122, &
            (/'ASSERTION ERROR:  obs not owned by any task'/))
      end if
   end do

   ! Create re_glob_slice and populate with residual_generic_type objects 
   ! allocated to match template.  
   call da_y_facade_create( re_glob_slice, num_obs_all, num_obs_all, template=template )

   ! NOTE:  This i loop should be inside the n loops.  
   ! Ideally, residual_generic class should know how to pack/unpack 
   ! (serialize/deserialize) itself for arbitrary I/O or communications (MPI or 
   ! otherwise) that clients may choose to implement.  Below are a possible set 
   ! of new methods for residual_generic_type:  
   !  residual_generic_pack_counts( res_generic, (out)rcount, (out)icount )
   !  [client allocates ibuf(icount) and rbuf(rcount)]
   !  residual_generic_pack( res_generic, (inout)rstart, (inout)rbuf, &
   !                         (inout)istart, (inout)ibuf )
   !  [client MPI communications:   ibuf->ibufg  rbuf->rbufg]
   !  residual_generic_unpack_counts( ibufg, (out)rstarts, (out)rcounts )
   !  residual_generic_unpack( (out)res_generic, rstarts(n), rcounts(n), rbufg )
   ! TOdo:  
   !  1) Design serialization for residual_generic_type.  
   !  2) Implement new methods.  
   !  3) Refactor below...  
   !  4) Optimize performance.  
   ! At the moment this refactoring is relatively low-priority since 
   ! da_y_facade_to_global() is already well-encapsulated and peformance is not 
   ! (yet) a concern for testing.  
   ! Loop over vector and scalar elements
   
   do i=template%lbnd,template%ubnd
      num_obs_send = 0
      num_values = 0
      num_send_buf_lcl = 0
      ! collect information to allocate buffers
      do n=1, re_slice%num_obs
         if ( re_slice%obs(n)%proc_domain ) then
            num_obs_send = num_obs_send + 1
            ! track number of scalars/levels per obs element
            num_values(num_obs_send) = size(re_slice%obs(n)%values(i)%ptr)
            ! track total buffer size
            num_send_buf_lcl = num_send_buf_lcl + num_values(num_obs_send)
         end if
      end do
      ! gather num_send_buf_lcl
      call mpi_allgather( num_send_buf_lcl, 1, mpi_integer, &
                          num_send_buf_all, 1,              &
                          mpi_integer, comm, ierr )
      ! gather num_values
      recv_counts = num_obs_send_all
      recv_displs(0) = 0
      do proc=1, num_procs-1
         recv_displs(proc) = recv_displs(proc-1) + recv_counts(proc-1)
      end do
      ! num_values
      call mpi_allgatherv( num_values, num_obs_send, mpi_integer,    &
                           num_values_all, recv_counts, recv_displs, &
                           mpi_integer, comm, ierr )
      ! set up to gather local arrays
      ! compute start index of each obs in gathered buffers
      obs_start_all(1) = 1
      do n=2, num_obs_all
         obs_start_all(n) = obs_start_all(n-1) + num_values_all(n-1)
      end do
      ! finish setting up to gather local arrays
      recv_counts = num_send_buf_all
      recv_displs(0) = 0
      do proc=1, num_procs-1
         recv_displs(proc) = recv_displs(proc-1) + recv_counts(proc-1)
      end do
      num_recv_all = sum( recv_counts )
      ! allocate and initialize send buffer
      allocate( re_vals_lcl(num_send_buf_lcl) )
      buf_idx = 0
      do n=1, re_slice%num_obs
         if ( re_slice%obs(n)%proc_domain ) then
            do k=1, size(re_slice%obs(n)%values(i)%ptr)
               buf_idx = buf_idx + 1
               re_vals_lcl(buf_idx) = re_slice%obs(n)%values(i)%ptr(k)
            end do
         end if
      end do
      if ( buf_idx /= num_send_buf_lcl ) then
         call da_error ("da_y_facade_to_global.inc",206, &
            (/'ASSERTION ERROR:  buf_idx /= num_send_buf_lcl'/))
      end if
      ! allocate receive buffer
      allocate( re_vals_all(num_recv_all) )
      ! finally, actually gather the values
      call mpi_allgatherv( re_vals_lcl, num_send_buf_lcl, true_mpi_real,   &
                           re_vals_all, recv_counts, recv_displs, &
                           true_mpi_real, comm, ierr )
      ! initialize re_glob_slice
      ! NOTE:  The i loop should be inside this n loop, see comment above...  
      ! TOdo:  Refactor...  
      do n=1, re_glob_slice%num_obs
         serial_n = obs_global_index_inv(n)
         serial_numvals = num_values_all(serial_n)
         buf_idx = obs_start_all(serial_n)
         ! note that we only collected obs someone owns
         re_glob_slice%obs(n)%proc_domain = .true.
         re_glob_slice%obs(n)%obs_global_index = obs_global_index_all(serial_n)
         call da_res_generic_alloc_and_set (re_glob_slice%obs(n), i, &
            re_vals_all(buf_idx:(buf_idx+serial_numvals-1)))
      end do
      ! deallocate communication buffers, etc.
      deallocate (re_vals_all)
      deallocate (re_vals_lcl) 
   end do  ! end of i loop

   call da_y_facade_free (re_slice)

   deallocate (num_values_all, obs_global_index_all, obs_global_index_inv, obs_start_all)
   deallocate (num_obs_send_all, num_send_buf_all, recv_counts, recv_displs)

   if (trace_use) call da_trace_exit("da_y_facade_to_global")


end subroutine da_y_facade_to_global


subroutine da_system(cmd)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   character (len=*) , intent(in) :: cmd

   call system(cmd)

end subroutine da_system



subroutine da_patch_to_global_dual_res(grid, vlocal, vglobal, mz)

   !----------------------------------------------------------------------
   ! Purpose: Gathers local 3D array vlocal into global array vglobal.  
   ! Assumes that "k" is not decomposed.  End indices in the "k" dimension 
   ! are inferred from mz, which can be less than kde.  
   !
   ! Must be called by all MPI tasks.  
   !----------------------------------------------------------------------

   implicit none

   type(domain),      intent(in)  :: grid
   real,              intent(in)  :: vlocal(ims_int:ime_int,jms_int:jme_int,kms_int:kme_int)
   real,              intent(out) :: vglobal(ids_int:ide_int,jds_int:jde_int,kds_int:kde_int)
   integer, optional, intent(in)  :: mz


   integer :: local_kde
   integer :: local_kme
   integer :: local_kpe

   if (trace_use_frequent) call da_trace_entry("da_patch_to_global_dual_res")


   if (present(mz)) then
      local_kde = kds_int + mz - 1
      local_kme = local_kde
      local_kpe = local_kde
   else
      local_kde = kde_int
      local_kme = kme_int
      local_kpe = kpe_int
   end if

   if (local_kde > 0) then
      call wrf_patch_to_global_real (vlocal, vglobal, grid%xp%domdesc, &
         trim(grid_stagger), trim(grid_ordering), &
         ids_int, ide_int, jds_int, jde_int, kds_int, local_kde,  &
         ims_int, ime_int, jms_int, jme_int, kms_int, local_kme,  &
         ips_int, ipe_int, jps_int, jpe_int, kps_int, local_kpe)
   end if

   if (trace_use_frequent) call da_trace_exit("da_patch_to_global_dual_res")

end subroutine da_patch_to_global_dual_res



subroutine da_proc_stats_combine(proc_ave, proc_err, proc_min, proc_max, &
   nobs_min, nobs_max, klev_min, klev_max)

   !---------------------------------------------------------------------------
   !  Purpose: Do MPI reduction operations across processors to get the average, 
   !           rms error, minimum, and maximum values for an observation field.
   !           These are stored only on the root processor, i.e., processor 0.
   !           (In this way, we do not have to do all-to-all communication.)
   !---------------------------------------------------------------------------

   implicit none

   real,      intent(inout)      :: proc_ave       ! Processor average.
   real,      intent(inout)      :: proc_err       ! Processor rms error.
   real,      intent(inout)      :: proc_min       ! Processor minumum.
   real,      intent(inout)      :: proc_max       ! Processor maximum.
   integer,   intent(inout)      :: nobs_min       ! Obs number of minimum.
   integer,   intent(inout)      :: nobs_max       ! Obs number of maximum.
   integer,   intent(inout)      :: klev_min       ! Level of minimum.
   integer,   intent(inout)      :: klev_max       ! Level of maximum.

   real    :: average            ! Global average.
   real    :: rms_err            ! Global rms_error.
   real    :: in(2)              ! mpi_reduce input value with processor rank.
   real    :: out(2)             ! mpi_reduce output min/max with associated processor.
   integer :: proc_id            ! Id of processor with max or min value.
   integer :: status(mpi_status_size) ! MPI status.
   real    :: buf(1)


   if (trace_use_frequent) call da_trace_entry("da_proc_stats_combine")

   ! Sum average and rms error over all processors and store on monitor processor.
   call mpi_reduce(proc_ave, average, 1, true_mpi_real, mpi_sum, root, comm, ierr)
   call mpi_reduce(proc_err, rms_err, 1, true_mpi_real, mpi_sum, root, comm, ierr)

   if (rootproc) then
      proc_ave = average
      proc_err = rms_err
   end if

   ! Get minimum value and associated processor index.
   in(1) = proc_min
   in(2) = myproc
   call mpi_reduce(in, out, 1, mpi_2double_precision, mpi_minloc, root, comm, ierr)

   if (myproc == root) then
      proc_min = out(1)
      proc_id = inT(out(2))
      buf(1) = real(proc_id)
   end if

   call wrf_dm_bcast_real (buf, 1)
   proc_id=int(buf(1))

   ! Get obs number and k-level where minimum occurs.
   if (proc_id .ne. root) then
      if (rootproc) then
         call mpi_recv(nobs_min, 1, mpi_integer, proc_id, 10, comm, STATUS, ierr)
         call mpi_recv(klev_min, 1, mpi_integer, proc_id, 11, comm, STATUS, ierr)
      else if (myproc == proc_id) then
         call mpi_send(nobs_min, 1, mpi_integer, root, 10, comm, ierr)
         call mpi_send(klev_min, 1, mpi_integer, root, 11, comm, ierr)
      end if
   end if

   ! Get maximum value and associated processor index.
   in(1) = proc_max
   in(2) = myproc

   call mpi_reduce(in, out, 1, mpi_2double_precision, mpi_maxloc, root, comm, ierr)

   if (rootproc) then
      proc_max = out(1)
      proc_id = int(out(2))
      buf(1) = real(proc_id)
   end if

   call wrf_dm_bcast_real (buf, 1)
   proc_id=int(buf(1))

   ! Get obs number and k-level where maximum occurs.
   if (proc_id .ne. root) then
      if (rootproc) then
         call mpi_recv(nobs_max, 1, mpi_integer, proc_id, 10, comm, STATUS, ierr)
         call mpi_recv(klev_max, 1, mpi_integer, proc_id, 11, comm, STATUS, ierr)
      else if (myproc == proc_id) then
         call mpi_send(nobs_max, 1, mpi_integer, root, 10, comm, ierr)
         call mpi_send(klev_max, 1, mpi_integer, root, 11, comm, ierr)
      end if
   end if

   if (trace_use_frequent) call da_trace_exit("da_proc_stats_combine")

end subroutine da_proc_stats_combine


subroutine da_proc_maxmin_combine(n, max, min)

   !---------------------------------------------------------------------------
   !  Purpose: Do MPI reduction operations across processors to get the minimum
   !           and maximum values for an observation field of length n. The
   !           i, j location of the minimum and maximum, for each n, is also
   !           communicated.
   !           The return values are stored only on the root processor, i.e., 
   !           processor 0.  (In this way, we do not have to do all-to-all 
   !           communication.)
   !---------------------------------------------------------------------------
   
   implicit none

   integer,                  intent(in)     :: n       ! Length of input fields.
   type (maxmin_field_type), intent(inout)  :: max(n)  ! Max values over proc.
   type (maxmin_field_type), intent(inout)  :: min(n)  ! Min values over proc.

   real    :: in(2*n)            ! mpi_reduce input value with processor rank.
   real    :: out(2*n)           ! mpi_reduce output min/max with associated processor.
   integer :: i                  ! Loop counter.
   integer :: proc_id(n)         ! Id of processor with max or min value.
   integer :: status(mpi_status_size) ! MPI status.


   if (trace_use_frequent) call da_trace_entry("da_proc_maxmin_combine")

   ! Get minimum value and associated processor index.
   do i = 1, n
      in(2*i-1) = min(i)%value
      in(2*i)   = myproc
   end do

   call mpi_reduce(in, out, n, mpi_2double_precision, mpi_minloc, root, comm, ierr)

   if (rootproc) then
      do i = 1, n
         min(i)%value = out(2*i-1)
         proc_id(i)   = int(out(2*i))
      end do
   end if

   call wrf_dm_bcast_integer (proc_id, n)

   ! Get i and j where minimum occurs.
   do i = 1, n
      if (proc_id(i) .ne. 0) then
         if (rootproc) then
            call mpi_recv(min(i)%i, 1, mpi_integer, proc_id(i), 10, comm, STATUS, ierr)
            call mpi_recv(min(i)%j, 1, mpi_integer, proc_id(i), 11, comm, STATUS, ierr)
         else if (myproc == proc_id(i)) then
            call mpi_send(min(i)%i, 1, mpi_integer, root, 10, comm, ierr)
            call mpi_send(min(i)%j, 1, mpi_integer, root, 11, comm, ierr)
         end if
      end if
   end do

   ! Get maximum value and associated processor index.
   do i = 1, n
      in(2*i-1) = max(i)%value
      in(2*i)   = myproc
   end do
   call mpi_reduce(in, out, n, mpi_2double_precision, mpi_maxloc, root, comm, ierr)

   if (rootproc) then
      do i = 1, n
         max(i)%value = out(2*i-1)
         proc_id(i)   = int(out(2*i))
      end do
   end if

   call wrf_dm_bcast_integer (proc_id, n)

   ! Get i and j where maximum occurs.
   do i = 1, n
      if (proc_id(i) .ne. root) then
         if (rootproc) then
            call mpi_recv(max(i)%i, 1, mpi_integer, proc_id(i), 10, comm, STATUS, ierr)
            call mpi_recv(max(i)%j, 1, mpi_integer, proc_id(i), 11, comm, STATUS, ierr)
         else if (myproc == proc_id(i)) then
            call mpi_send(max(i)%i, 1, mpi_integer, root, 10, comm, ierr)
            call mpi_send(max(i)%j, 1, mpi_integer, root, 11, comm, ierr)
         end if
      end if
   end do

   if (trace_use_frequent) call da_trace_exit("da_proc_maxmin_combine")

end subroutine da_proc_maxmin_combine



end module da_par_util

