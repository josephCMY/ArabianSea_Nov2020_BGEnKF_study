












module da_recursive_filter

   
   
   

   use module_domain, only : domain,  vp_type
  
   use da_control, only : ims,ime,jms,jme,kms,kme,ids,ide,jds,jde,kds,kde, &
      rf_passes, its,ite,jts,jte,kts,kte,vert_corr, trace_use, vert_corr_1,&
      trace_use_dull, cv_size

   use da_control, only : its_int, ite_int, jts_int, jte_int, kts_int, kte_int, & 
                          ims_int, ime_int, jms_int, jme_int, kms_int, kme_int, &
                          ids_int, ide_int, jds_int, jde_int, kds_int, kde_int 

   use da_define_structures, only : be_type 
   use da_par_util, only : da_transpose_z2y, da_transpose_x2y, &
      da_transpose_y2z, da_transpose_y2x, da_transpose_x2z, &
      da_transpose_z2x, da_vv_to_cv,da_cv_to_vv
   use da_tracing, only : da_trace_entry, da_trace_exit

   use da_rfz_cv3, only : da_rfz
   use da_rf_cv3, only : smoothx, smoothy

   implicit none

   contains

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


subroutine da_calculate_rf_factors (rf_lengthscale, rf_alpha, rf_scale_factor)

   !---------------------------------------------------------------------------
   ! Purpose: Calculate:
   !          1) Alpha value for recursive filter.
   !          2) Turning conditions appropriate for B=UU^T RF.
   !          3) Generic (depends only on rf_passes) rescaling factor.
   !          4) Grid-dependent (hopefully temporary) scaling factor - removed.
   !---------------------------------------------------------------------------

   implicit none
   
   real*8, intent(in) :: rf_lengthscale(:)      ! Non-dim. R.F. lengthscale.
   real*8,intent(out) :: rf_alpha(:)            ! RF alpha factor.
   real*8,intent(out) :: rf_scale_factor(:)     ! Variance scaling factor.
      
   integer, parameter :: n = 500                ! 2n +1 = # pts in delta func.
   integer            :: kz                     ! 3rd array dimension.
   integer            :: pass                   ! Pass of recursive filter.
   integer            :: k                      ! Loop counter
   ! integer          :: nn                       ! Loop counter
      
   real               :: rf_e                   ! RF E factor.      

   real*8, dimension(-n:n) :: field_in, field_out  ! Working field.
   ! real, dimension(-n:n) :: field_out1  ! Working field.

   if (trace_use_dull) call da_trace_entry("da_calculate_rf_factors")

   !-------------------------------------------------------------------------
   ! [1.0]: Initialise:
   !-------------------------------------------------------------------------  

   kz = size(rf_scale_factor)
   
   rf_scale_factor(:) = 0.0
   
   do k = 1, kz

      !-------------------------------------------------------------------------
      ! [2.0]: Calculate RF alpha:
      !-------------------------------------------------------------------------  

      rf_e = 0.25 * rf_passes / (rf_lengthscale(k) * rf_lengthscale(k))
      rf_alpha(k) = 1.0 + rf_e - sqrt(rf_e * (rf_e + 2.0))

      !-------------------------------------------------------------------------
      ! [3.0]: Calculate rescaling factor:
      !-------------------------------------------------------------------------

      ! [3.1]: Calculate generic rescaling (normalise zero distance to 1):
      ! For rf_passes=2 (SOAR) = 4*rf_lengthscale.
      ! For rf_passes=infinity (Gaussian) = sqrt(8*pi)*rf_lengthscale.

      field_in(-n:n) = 0.0
      field_in(0) = 1.0
      field_out(-n:n) = field_in(-n:n)

      do pass = 1, rf_passes / 2
         call da_recursive_filter_1d_adj(pass, rf_alpha(k), field_out, 2*n+1)
      end do

      do pass = 1, rf_passes / 2
         call da_recursive_filter_1d(pass, rf_alpha(k), field_out, 2*n+1)
      end do

      rf_scale_factor(k) = 1.0 / field_out(0)

      ! Uncomment the following to test equivalence of UU^T and RF:
      ! write(unit=stdout,fmt='(A,f15.5)') &
      !    ' RF Scaling Factor = ', 1.0 / field_out(0)
      ! field_out1(-n:n) = field_in(-n:n)
      ! do pass = 1, rf_passes
      !    call da_recursive_filter_1d(pass, rf_alpha(k), field_out1, 2*n+1)
      ! end do

      ! do nn = -n, n
      !    write(unit=stdout,fmt='(2i5,4f12.5)')k, nn, field_in(nn), &
      !                             field_out(nn) / field_out(0), &
      !                             field_out1(nn) / field_out1(0), &
      !                             exp(-0.125*(real(nn)/rf_lengthscale(k))**2)
      ! end do

   end do ! End loop over k

   if (trace_use_dull) call da_trace_exit("da_calculate_rf_factors")
   
end subroutine da_calculate_rf_factors
      

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


subroutine da_recursive_filter_1d_adj(pass, alpha, field, n)

   !---------------------------------------------------------------------------
   ! Purpose: Perform one pass of recursive filter on 1D array - adjoint.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: pass           ! Current pass of filter.
   real*8, intent(in)     :: alpha          ! Alpha coefficient for RF.
   real*8, intent(inout)  :: field(:)       ! Array to be filtered.
   integer, intent(in)    :: n              ! Size of field array.
   
   integer                :: j              ! Loop counter.
   real                   :: one_alpha      ! 1 - alpha.
   real                   :: a(1:n)         ! Input field.
   real                   :: b(1:n)         ! Field after left-right pass.
   real                   :: c(1:n)         ! Field after right-left pass.

   if (trace_use_dull) call da_trace_entry("da_recursive_filter_1d_adj")

   !-------------------------------------------------------------------------
   ! [1.0] Initialise:
   !-------------------------------------------------------------------------

   one_alpha = 1.0 - alpha

   !-------------------------------------------------------------------------
   ! [4.0] Copy and tidy up:
   !-------------------------------------------------------------------------

   c(1:n) = field(1:n)
   
   !-------------------------------------------------------------------------
   ! [3.0] Perform left-moving filter:
   !-------------------------------------------------------------------------

   ! [3.2] Perform pass left to right:

   b(1:n) = 0.0   

   do j = 1, n-1
      c(j+1) = c(j+1) + alpha * c(j)
      b(j) = one_alpha * c(j)
   end do

   ! use turning conditions as in the appendix of Hayden & Purser (1995):

   if (pass == 1) then
      b(n) = b(n) + c(n) / (1.0 + alpha)
   else
      b(n) = b(n) + one_alpha * c(n) / (1.0 - alpha**2)**2
      b(n-1) = b(n-1) - one_alpha * alpha**3 * c(n) / (1.0 - alpha**2)**2
   end if

   !-------------------------------------------------------------------------
   ! [2.0] Perform right-moving filter:
   !-------------------------------------------------------------------------

   a(1:n) = 0.0

   ! [2.2] Perform pass left to right:

   do j = n, 2, -1
      b(j-1) = b(j-1) + alpha * b(j)
      a(j) = a(j) + one_alpha * b(j)
   end do

   ! use turning conditions as in the appendix of Hayden & Purser (1995):

   if (pass == 1) then
      a(1) = a(1) + one_alpha * b(1)
   else if (pass == 2) then
      a(1) = a(1) + b(1) / (1.0 + alpha)
   else
      a(1) = a(1) + one_alpha * b(1) / (1.0 - alpha**2)**2
      a(2) = a(2) - one_alpha * alpha**3 * b(1) / (1.0 - alpha**2)**2
   end if

   field(1:n) = a(1:n)

   if (trace_use_dull) call da_trace_exit("da_recursive_filter_1d_adj")
   
end subroutine da_recursive_filter_1d_adj


subroutine da_transform_through_rf(grid,mz, rf_alpha, val,field, scaling)

   !---------------------------------------------------------------------------
   ! Purpose: Control variable transform through the recursive filter.  
   !
   ! Method: 1) Transform to nondimensional v_hat space.
   !         2) Perform recursive filter in non-dimensional space (ds = 1).
   !         3) Scale recursive filter output.
   !         4) Transform filtered field to dimensional space.
   !         5) Optionally scale by background error.
   !---------------------------------------------------------------------------

   implicit none

   type(domain),     intent(inout) :: grid
   integer,          intent(in)    :: mz                ! Vertical truncation.
   real*8,           intent(in)    :: rf_alpha(mz)      ! RF scale parameter. 
   real*8,           intent(in)    :: val(jds:jde,mz)   ! Error standard deviation.
   real,             intent(inout) :: field(ims:ime,jms:jme,kms:kme) ! Field to be transformed. 
      
   integer :: rf_passes_over_two ! rf_passes / 2
   integer :: i, j, m, n, pass, ij   ! Loop counters.
   real    :: p_x(ims:ime,jms:jme)! sqrt(Grid box area).
   real*8  :: val_j(grid%xp%jtsy:grid%xp%jtey)
   real*8  :: val_i(grid%xp%itsx:grid%xp%itex)
   
   logical, optional, intent(in) :: scaling

   !----------------------------------------------------------------------
   ! [1.0]: Initialise:
   !----------------------------------------------------------------------  

   if (trace_use_dull) call da_trace_entry("da_transform_through_rf")

   rf_passes_over_two = rf_passes / 2
   
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, i, j )
   do ij = 1 , grid%num_tiles
      ! [1.1] Define inner product (square root of grid box area):
      do j = grid%j_start(ij), grid%j_end(ij)
         do i = its, ite
            p_x(i,j) = sqrt(grid%xb%grid_box_area(i,j))
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, m, i, j )
   do ij = 1 , grid%num_tiles
      do m = 1, mz
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its, ite
               grid%xp%v1z(i,j,m) = 0.0
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

      ! [1.2] Transform to nondimensional v_hat space:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, m, i, j )
   do ij = 1 , grid%num_tiles
      do m = 1, mz
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its, ite
               grid%xp%v1z(i,j,m) = field(i,j,m) / p_x(i,j)
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   !-------------------------------------------------------------------------
   ! [2.0]: Perform 1D recursive filter in x-direction:
   !-------------------------------------------------------------------------

   ! [2.1] Apply (i',j',k -> i,j',k') (grid%xp%v1z -> grid%xp%v1x)
   ! convert from vertical column to x-stripe

   call da_transpose_z2x (grid)

   ! [2.2] Apply 1D filter in x direction:

   n = grid%xp%itex - grid%xp%itsx + 1
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( m, j, pass, val_i, i )
   do m = grid%xp%ktsx, min(grid%xp%ktex,mz)
      do j = grid%xp%jtsx, grid%xp%jtex
         do i = grid%xp%itsx, grid%xp%itex
            val_i(i) = grid%xp%v1x(i,j,m)
         end do
         do pass = 1, rf_passes_over_two
            call da_recursive_filter_1d(pass, rf_alpha(m), val_i, n)
         end do
         do i = grid%xp%itsx, grid%xp%itex
            grid%xp%v1x(i,j,m) = val_i(i)
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   !-------------------------------------------------------------------------
   ! [3.0]: Perform 1D recursive filter in y-direction:
   !-------------------------------------------------------------------------

   ! [3.1] Apply (i, j' ,k' -> i', j ,k') (grid%xp%v1x -> grid%xp%v1y)
   ! convert from vertical column to y-stripe

   call da_transpose_x2y (grid)

   ! [3.2] Apply 1D filter in y direction:

   n = grid%xp%jtey - grid%xp%jtsy + 1
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( m, i, pass, val_j, j )
   do m = grid%xp%ktsy, min(grid%xp%ktey,mz)
      do i = grid%xp%itsy, grid%xp%itey
         do j = grid%xp%jtsy, grid%xp%jtey
            val_j(j) = grid%xp%v1y(i,j,m)
         end do
         do pass = 1, rf_passes_over_two
            call da_recursive_filter_1d(pass, rf_alpha(m), val_j, n)
         end do
         do j = grid%xp%jtsy, grid%xp%jtey
            grid%xp%v1y(i,j,m) = val_j(j)
         end do
      end do
   end do
   !$OMP END PARALLEL DO
   
   !-------------------------------------------------------------------------
   ! [4.0]: Perform 1D recursive filter in y-direction:
   !-------------------------------------------------------------------------

   ! [4.1] Apply (i',j,k' -> i',j',k) (grid%xp%v1y -> grid%xp%v1z)
   ! convert from y-stripe to vertical column.
   
   call da_transpose_y2z (grid)

   ! [4.2] Transform filtered field to dimensional space:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, m, i, j )
   do ij = 1 , grid%num_tiles
      do m = 1, mz
         do j = grid%j_start(ij), grid%j_end(ij)
   !     do j = jts, jte
            do i = its, ite
               field(i,j,m) = grid%xp%v1z(i,j,m) * p_x(i,j)
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   ! [4.3] Optionally scale by background error:

   ! be_s % val = Gridpoint standard deviation - only required for
   ! vert_corr = vert_corr_1 as scaling is performed in vertical transform
   ! for vert_corr = vert_corr_2:

   if (vert_corr == vert_corr_1 .or. present(scaling)) then
      if (scaling .or. vert_corr == vert_corr_1) then
         do m = 1, mz
            do i = its, ite
               field(i,jts:jte,m) = field(i,jts:jte,m) * val(jts:jte,m)
            end do
         end do
      end if
   end if
   if (trace_use_dull) call da_trace_exit("da_transform_through_rf")

end subroutine da_transform_through_rf


subroutine da_transform_through_rf_adj(grid, mz,rf_alpha, val, field, scaling)

   !---------------------------------------------------------------------------
   ! Purpose: Control variable transform through the recursive filter.
   !
   ! Method: 1) Transform to nondimensional v_hat space.
   !         2) Perform recursive filter in non-dimensional space (ds = 1).
   !         3) Scale recursive filter output.
   !         4) Transform filtered field to dimensional space.
   !         5) Optionally scale by background error.
   !---------------------------------------------------------------------------

   implicit none

   type(domain),     intent(inout) :: grid
   integer,          intent(in)    :: mz                             ! Vertical truncation.
   real*8,           intent(in)    :: rf_alpha(mz)                   ! RF scale parameter. 
   real*8,           intent(in)    :: val(jds:jde,mz)                ! Error standard deviation.
   real,             intent(inout) :: field(ims:ime,jms:jme,kms:kme) ! Field to be transformed. 
      
   integer :: rf_passes_over_two    ! rf_passes / 2
   integer :: i, j, m, n, pass, ij  ! Loop counters.
   real    :: p_x(ims:ime,jms:jme)  ! sqrt(Grid box area).
   real*8  :: val_j(grid%xp%jtsy:grid%xp%jtey)
   real*8  :: val_i(grid%xp%itsx:grid%xp%itex)

   logical, optional, intent(in) :: scaling

   !-------------------------------------------------------------------------
   ! [1.0]: Initialise:
   !-------------------------------------------------------------------------  

   if (trace_use_dull) call da_trace_entry("da_transform_through_rf_adj")

   rf_passes_over_two = rf_passes / 2

   ! [1.1] Define inner product (square root of grid box area):
   !$OMP PARALLEL DO &
   !$OMP PRIVATE (ij, i, j)
   do ij = 1 , grid%num_tiles
      do j = grid%j_start(ij), grid%j_end(ij)
         do i = its, ite
            p_x(i,j) = sqrt(grid%xb%grid_box_area(i,j))
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   !-------------------------------------------------------------------------
   ! [4.0]: Perform 1D recursive filter in y-direction:
   !-------------------------------------------------------------------------
   
   ! [4.3] Optionally scale by background error:
   ! be_s % val = Gridpoint standard deviation - only required for
   ! vert_corr = vert_corr_1 as scaling is performed in vertical transform
   ! for vert_corr = vert_corr_2:

   if (vert_corr == vert_corr_1 .or. (present(scaling))) then
      if (scaling .or. vert_corr == vert_corr_1) then
         do m = 1, mz
            do i = its, ite
               field(i,jts:jte,m) = field(i,jts:jte,m) * val(jts:jte,m)
            end do
         end do
      end if
   end if

   ! [4.2] Transform filtered field to dimensional space:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (ij ,m, j, i)
   do ij = 1 , grid%num_tiles
      do m = 1, mz
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its, ite
               grid%xp%v1z(i,j,m) = field(i,j,m) * p_x(i,j)
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   ! [4.1] Apply (i',j',k -> i',j,k') (grid%xp%v1z -> grid%xp%v1y)
   ! convert vertical column to y-stripe

   call da_transpose_z2y (grid)

   !-------------------------------------------------------------------------
   ! [3.0]: Perform 1D recursive filter in y-direction:
   !-------------------------------------------------------------------------

   !  [3.2] Apply 1D filter in y direction:

   n=grid%xp%jtey-grid%xp%jtsy+1
   !$OMP PARALLEL DO &
   !$OMP PRIVATE (m, i, val_j, pass, j)
   do m = grid%xp%ktsy, min(grid%xp%ktey, mz)
      do i = grid%xp%itsy, grid%xp%itey
         do j = grid%xp%jtsy, grid%xp%jtey
            val_j(j) = grid%xp%v1y(i,j,m)
         end do
         do pass = rf_passes_over_two, 1, -1
            call da_recursive_filter_1d_adj(pass, rf_alpha(m), val_j, n)
         end do
         do j = grid%xp%jtsy, grid%xp%jtey
            grid%xp%v1y(i,j,m) = val_j(j)
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   ! [3.1] Apply (i',j,k' -> i,j',k') (grid%xp%v1y -> grid%xp%v1x)
   ! convert from y-stripe to x-stripe

   call da_transpose_y2x (grid)

   !-------------------------------------------------------------------------
   ! [2.0]: Perform 1D recursive filter in x-direction:
   !-------------------------------------------------------------------------

   ! [2.2] Apply 1D filter in x direction:

   n = grid%xp%itex-grid%xp%itsx+1

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( m, j, pass, i, val_i)
   do m = grid%xp%ktsx, min(grid%xp%ktex,mz)
      do j = grid%xp%jtsx, grid%xp%jtex
         do i = grid%xp%itsx, grid%xp%itex
            val_i(i) = grid%xp%v1x(i,j,m)
         end do
         do pass = rf_passes_over_two, 1, -1
            call da_recursive_filter_1d_adj(pass, rf_alpha(m), val_i, n)
         end do
         do i = grid%xp%itsx, grid%xp%itex
            grid%xp%v1x(i,j,m) = val_i(i) 
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   ! [2.1] Apply (i,j',k' -> i',j',k) (grid%xp%v1x -> grid%xp%v1z)
   ! convert from x-stripe to vertical column

   call da_transpose_x2z (grid)

   !-------------------------------------------------------------------------
   ! [1.0]: Initialise:
   !-------------------------------------------------------------------------  

   ! [1.2] Transform to nondimensional v_hat space:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (ij ,m, i, j)
   do ij = 1 , grid%num_tiles
      do m = 1, mz
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its, ite
               field(i,j,m) = grid%xp%v1z(i,j,m) / p_x(i,j)
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   if (trace_use_dull) call da_trace_exit("da_transform_through_rf_adj")

end subroutine da_transform_through_rf_adj



SUBROUTINE da_apply_rf_1v( be, vp, grid, nv )

   IMPLICIT NONE

   TYPE (be_type), INTENT(IN)       :: be     ! Background error structure.
   type (domain) , intent(inout)    :: grid   ! Dimensions and xpose buffers.

   integer, intent(in) :: nv                        ! # of var.

   real, dimension(ims:ime, jms:jme, kms:kme), INTENT(INOUT) :: vp   ! working array

   integer                 :: in, jn, kn, k
!-------------------------------------------------------------------------
!  [1.0] Initialise:
!-------------------------------------------------------------------------

   in=ite-its+1
   jn=jte-jts+1
   kn=kte-kts+1

   call da_rfz(vp(ims:ime,jms:jme,kms:kme),in,jn,kn,be%ndeg,&
     be%vz(kts:kte,its:ite,jts:jte,nv),be%be,be%table,be%nta,be%swidth, &
                                    ids,ide, jds,jde, kds,kde,  &
                                    ims,ime, jms,jme, kms,kme,  &
                                    its,ite, jts,jte, kts,kte )

!-------------------------------------------------------------------------
!  [2.0]: Perform 1D recursive filter in y-x direction:
!-------------------------------------------------------------------------

   do k = kts,kte
      grid%xp % v1z(its:ite,jts:jte,k) = vp(its:ite,jts:jte,k)
   end do

   call da_transpose_z2x ( grid )

   in=grid%xp%ipex-grid%xp%ipsx
   jn=grid%xp%jpex-grid%xp%jpsx

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( k )
   do k=grid%xp%kpsx,grid%xp%kpex
     call smoothx(in,jn, &
          grid%xp% v1x(grid%xp%ipsx:grid%xp%ipex,grid%xp%jpsx:grid%xp%jpex,k),&
          be%slix(grid%xp%ipsx:grid%xp%ipex,grid%xp%jpsx:grid%xp%jpex,k,nv),  &
          be%ndeg,be%be,be%nta,be%swidth,be%table)
   enddo
   !$OMP END PARALLEL DO


   call da_transpose_x2y ( grid )

   in=grid%xp%ipey-grid%xp%ipsy
   jn=grid%xp%jpey-grid%xp%jpsy

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( k )
   do k=grid%xp%kpsy,grid%xp%kpey
     call smoothy(in,jn, &
          grid%xp%v1y(grid%xp%ipsy:grid%xp%ipey,grid%xp%jpsy:grid%xp%jpey,k),&
          be%sljy(grid%xp%ipsy:grid%xp%ipey,grid%xp%jpsy:grid%xp%jpey,k,nv), &
          be%ndeg,be%be,be%nta,be%swidth,be%table)
   enddo
   !$OMP END PARALLEL DO


   call da_transpose_y2z ( grid )

   do k = kts,kte
      vp(its:ite,jts:jte,k)= grid%xp % v1z(its:ite,jts:jte,k)
   end do

END SUBROUTINE da_apply_rf_1v

SUBROUTINE da_apply_rf_1v_adj( be, vp, grid, nv )

   IMPLICIT NONE

   TYPE (be_type), INTENT(IN)       :: be     ! Background error structure.
   type (domain) , intent(inout)    :: grid   ! Dimensions and xpose buffers.

   integer, intent(in) :: nv                        ! # of var.

   real, dimension(ims:ime, jms:jme, kms:kme), INTENT(INOUT) :: vp   ! working array

   integer             :: in, jn, kn, k

!-------------------------------------------------------------------------
!  [1.0] Perform 1D recursive filter in x-y direction:
!-------------------------------------------------------------------------

   do k = kts,kte
      grid%xp % v1z(its:ite,jts:jte,k) = vp(its:ite,jts:jte,k)
   end do

   call da_transpose_z2y ( grid )

   in=grid%xp%ipey-grid%xp%ipsy
   jn=grid%xp%jpey-grid%xp%jpsy

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( k )
   do k=grid%xp%kpsy,grid%xp%kpey
      call smoothy(in,jn,&
               grid%xp % v1y(grid%xp%ipsy:grid%xp%ipey,grid%xp%jpsy:grid%xp%jpey,k),&
               be%sljy(grid%xp%ipsy:grid%xp%ipey,grid%xp%jpsy:grid%xp%jpey,k,nv),   &
               be%ndeg,be%be,be%nta,be%swidth,be%table)
   enddo
   !$OMP END PARALLEL DO

   call da_transpose_y2x ( grid )

   in=grid%xp%ipex-grid%xp%ipsx
   jn=grid%xp%jpex-grid%xp%jpsx

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( k )
   do k=grid%xp%kpsx,grid%xp%kpex
      call smoothx(in,jn, &
               grid%xp % v1x(grid%xp%ipsx:grid%xp%ipex,grid%xp%jpsx:grid%xp%jpex,k),&
               be%slix(grid%xp%ipsx:grid%xp%ipex,grid%xp%jpsx:grid%xp%jpex,k,nv),   &
               be%ndeg,be%be,be%nta,be%swidth,be%table)
   enddo
   !$OMP END PARALLEL DO

   call da_transpose_x2z ( grid )

   do k = kts,kte
      vp(its:ite,jts:jte,k)= grid%xp % v1z(its:ite,jts:jte,k)
   end do

!-------------------------------------------------------------------------
!  [2.0]: Perform 1D recursive filter in z direction:
!-------------------------------------------------------------------------

   in=ite-its+1
   jn=jte-jts+1
   kn=kte-kts+1

   call da_rfz(vp(ims:ime,jms:jme,kms:kme),in,jn,kn,be%ndeg,&
     be%vz(kts:kte,its:ite,jts:jte,nv),be%be,be%table,be%nta,be%swidth,&
                                    ids,ide, jds,jde, kds,kde,  &
                                    ims,ime, jms,jme, kms,kme,  &
                                    its,ite, jts,jte, kts,kte )
END SUBROUTINE da_apply_rf_1v_adj

SUBROUTINE da_apply_rf( be, vp, grid )

   IMPLICIT NONE

   TYPE (be_type), INTENT(IN)       :: be   ! Background error structure.
   TYPE (vp_type), INTENT(INOUT)    :: vp   ! working array
   type (domain) , intent(inout)    :: grid   ! Dimensions and xpose buffers.

   integer :: in, jn

   integer :: i, j, k

!-------------------------------------------------------------------------

   call da_apply_rf_1v( be, vp%v1, grid, 1)

!-------------------------------------------------------------------------

   call da_apply_rf_1v( be, vp%v2, grid, 2)

!-------------------------------------------------------------------------

   call da_apply_rf_1v( be, vp%v3, grid, 3)

!-------------------------------------------------------------------------

   call da_apply_rf_1v( be, vp%v4, grid, 4)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!  [2.0]: Perform 1D recursive filter in y-x direction:
!-------------------------------------------------------------------------

   grid%xp%v1z(its:ite,jts:jte,1) = vp%v5(its:ite,jts:jte,1)

   call da_transpose_z2x ( grid )

   in=grid%xp%ipex-grid%xp%ipsx
   jn=grid%xp%jpex-grid%xp%jpsx

   if ( LBOUND(grid%xp%v1x,3) == 1 ) then
      call smoothx(in,jn,&
           grid%xp%v1x(grid%xp%ipsx:grid%xp%ipex,grid%xp%jpsx:grid%xp%jpex,1),&
           be%slipx(grid%xp%ipsx:grid%xp%ipex,grid%xp%jpsx:grid%xp%jpex),&
           be%ndeg,be%be,be%nta,be%swidth,be%table)
   endif

   call da_transpose_x2y ( grid )

   in=grid%xp%ipey-grid%xp%ipsy
   jn=grid%xp%jpey-grid%xp%jpsy

   if ( LBOUND(grid%xp%v1y,3) == 1 ) then
      call smoothy(in,jn, &
           grid%xp%v1y(grid%xp%ipsy:grid%xp%ipey,grid%xp%jpsy:grid%xp%jpey,1),&
           be%sljpy(grid%xp%ipsy:grid%xp%ipey,grid%xp%jpsy:grid%xp%jpey),&
           be%ndeg,be%be,be%nta,be%swidth,be%table)
   endif

   call da_transpose_y2z ( grid )

   vp%v5(its:ite,jts:jte,1)= grid%xp % v1z(its:ite,jts:jte,1)

!-------------------------------------------------------------------------

END SUBROUTINE da_apply_rf          

SUBROUTINE da_apply_rf_adj( be, vp , grid )

   IMPLICIT NONE

   TYPE (be_type), INTENT(IN)       :: be     ! Background error structure.
   TYPE (vp_type), INTENT(INOUT)    :: vp     ! working array
   type (domain) , intent(inout)    :: grid   ! Dimensions and xpose buffers.

   integer             :: i, j, k
   integer             :: in, jn, kn

!-------------------------------------------------------------------------

   call da_apply_rf_1v_adj( be, vp%v1, grid, 1 )
!-------------------------------------------------------------------------

   call da_apply_rf_1v_adj( be, vp%v2, grid, 2 )

!-------------------------------------------------------------------------

   call da_apply_rf_1v_adj( be, vp%v3, grid, 3 )

!-------------------------------------------------------------------------

   call da_apply_rf_1v_adj( be, vp%v4, grid, 4 )

!-------------------------------------------------------------------------
!  [1.0] Perform 1D recursive filter in x-y direction:
!-------------------------------------------------------------------------

   grid%xp % v1z(its:ite,jts:jte,1) = vp%v5(its:ite,jts:jte,1)

   call da_transpose_z2y ( grid )

   in=grid%xp%ipey-grid%xp%ipsy
   jn=grid%xp%jpey-grid%xp%jpsy
 
   if ( LBOUND(grid%xp%v1y,3) == 1 ) then
      call smoothy(in,jn, &
           grid%xp%v1y(grid%xp%ipsy:grid%xp%ipey,grid%xp%jpsy:grid%xp%jpey,1),&
           be%sljpy(grid%xp%ipsy:grid%xp%ipey,grid%xp%jpsy:grid%xp%jpey),&
           be%ndeg,be%be,be%nta,be%swidth,be%table)
   endif

   call da_transpose_y2x ( grid )

   in=grid%xp%ipex-grid%xp%ipsx
   jn=grid%xp%jpex-grid%xp%jpsx

   if ( LBOUND(grid%xp%v1x,3) == 1 ) then
      call smoothx(in,jn,&
           grid%xp%v1x(grid%xp%ipsx:grid%xp%ipex,grid%xp%jpsx:grid%xp%jpex,1),&
           be%slipx(grid%xp%ipsx:grid%xp%ipex,grid%xp%jpsx:grid%xp%jpex),&
           be%ndeg,be%be,be%nta,be%swidth,be%table)
   endif

   call da_transpose_x2z ( grid )

   vp%v5(its:ite,jts:jte,1)= grid%xp % v1z(its:ite,jts:jte,1)

!-------------------------------------------------------------------------

END SUBROUTINE da_apply_rf_adj      


subroutine da_transform_through_rf_dual_res(grid,mz, rf_alpha, val,field, scaling)

   !---------------------------------------------------------------------------
   ! Purpose: Control variable transform through the recursive filter.  
   !
   ! Method: 1) Transform to nondimensional v_hat space.
   !         2) Perform recursive filter in non-dimensional space (ds = 1).
   !         3) Scale recursive filter output.
   !         4) Transform filtered field to dimensional space.
   !         5) Optionally scale by background error.
   !---------------------------------------------------------------------------

   implicit none

   type(domain),     intent(inout) :: grid
   integer,          intent(in)    :: mz                ! Vertical truncation.
   real*8,           intent(in)    :: rf_alpha(mz)      ! RF scale parameter. 
   real*8,           intent(in)    :: val(jds_int:jde_int,mz)   ! Error standard deviation.
!  real*8,           intent(in)    :: val(1:jde_int-jds_int+1,mz)   ! Error standard deviation.
   real,             intent(inout) :: field(ims_int:ime_int,jms_int:jme_int,kms_int:kme_int) ! Field to be transformed.
      
   integer :: rf_passes_over_two ! rf_passes / 2
   integer :: i, j, m, n, pass, ij   ! Loop counters.
   real    :: p_x(ims_int:ime_int,jms_int:jme_int)! sqrt(Grid box area).
   real*8  :: val_j(grid%xp%jtsy:grid%xp%jtey)
   real*8  :: val_i(grid%xp%itsx:grid%xp%itex)
   
   logical, optional, intent(in) :: scaling

   !----------------------------------------------------------------------
   ! [1.0]: Initialise:
   !----------------------------------------------------------------------  

   if (trace_use_dull) call da_trace_entry("da_transform_through_rf")

   rf_passes_over_two = rf_passes / 2
   
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, i, j )
   do ij = 1 , grid%num_tiles
      ! [1.1] Define inner product (square root of grid box area):
      do j = grid%j_start(ij), grid%j_end(ij)
         do i = its_int, ite_int
            p_x(i,j) = sqrt(grid%xb%grid_box_area(i,j))
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, m, i, j )
   do ij = 1 , grid%num_tiles
      do m = 1, mz
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its_int, ite_int
               grid%xp%v1z(i,j,m) = 0.0
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

      ! [1.2] Transform to nondimensional v_hat space:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, m, i, j )
   do ij = 1 , grid%num_tiles
      do m = 1, mz
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its_int, ite_int
               grid%xp%v1z(i,j,m) = field(i,j,m) / p_x(i,j)
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   !-------------------------------------------------------------------------
   ! [2.0]: Perform 1D recursive filter in x-direction:
   !-------------------------------------------------------------------------

   ! [2.1] Apply (i',j',k -> i,j',k') (grid%xp%v1z -> grid%xp%v1x)
   ! convert from vertical column to x-stripe

   call da_transpose_z2x (grid)

   ! [2.2] Apply 1D filter in x direction:

   n = grid%xp%itex - grid%xp%itsx + 1
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( m, j, pass, val_i, i )
   do m = grid%xp%ktsx, min(grid%xp%ktex,mz)
      do j = grid%xp%jtsx, grid%xp%jtex
         do i = grid%xp%itsx, grid%xp%itex
            val_i(i) = grid%xp%v1x(i,j,m)
         end do
         do pass = 1, rf_passes_over_two
            call da_recursive_filter_1d(pass, rf_alpha(m), val_i, n)
         end do
         do i = grid%xp%itsx, grid%xp%itex
            grid%xp%v1x(i,j,m) = val_i(i)
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   !-------------------------------------------------------------------------
   ! [3.0]: Perform 1D recursive filter in y-direction:
   !-------------------------------------------------------------------------

   ! [3.1] Apply (i, j' ,k' -> i', j ,k') (grid%xp%v1x -> grid%xp%v1y)
   ! convert from vertical column to y-stripe

   call da_transpose_x2y (grid)

   ! [3.2] Apply 1D filter in y direction:

   n = grid%xp%jtey - grid%xp%jtsy + 1
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( m, i, pass, val_j, j )
   do m = grid%xp%ktsy, min(grid%xp%ktey,mz)
      do i = grid%xp%itsy, grid%xp%itey
         do j = grid%xp%jtsy, grid%xp%jtey
            val_j(j) = grid%xp%v1y(i,j,m)
         end do
         do pass = 1, rf_passes_over_two
            call da_recursive_filter_1d(pass, rf_alpha(m), val_j, n)
         end do
         do j = grid%xp%jtsy, grid%xp%jtey
            grid%xp%v1y(i,j,m) = val_j(j)
         end do
      end do
   end do
   !$OMP END PARALLEL DO
   
   !-------------------------------------------------------------------------
   ! [4.0]: Perform 1D recursive filter in y-direction:
   !-------------------------------------------------------------------------

   ! [4.1] Apply (i',j,k' -> i',j',k) (grid%xp%v1y -> grid%xp%v1z)
   ! convert from y-stripe to vertical column.
   
   call da_transpose_y2z (grid)

   ! [4.2] Transform filtered field to dimensional space:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, m, i, j )
   do ij = 1 , grid%num_tiles
      do m = 1, mz
         do j = grid%j_start(ij), grid%j_end(ij)
   !     do j = jts, jte
            do i = its_int, ite_int
               field(i,j,m) = grid%xp%v1z(i,j,m) * p_x(i,j)
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   ! [4.3] Optionally scale by background error:

   ! be_s % val = Gridpoint standard deviation - only required for
   ! vert_corr = vert_corr_1 as scaling is performed in vertical transform
   ! for vert_corr = vert_corr_2:

   if (vert_corr == vert_corr_1 .or. present(scaling)) then
      if (scaling .or. vert_corr == vert_corr_1) then
         do m = 1, mz
            do i = its_int, ite_int
               field(i,jts_int:jte_int,m) = field(i,jts_int:jte_int,m) * val(jts_int:jte_int,m)
            end do
         end do
      end if
   end if
   if (trace_use_dull) call da_trace_exit("da_transform_through_rf")

end subroutine da_transform_through_rf_dual_res


subroutine da_transform_through_rf_adj_dual_res(grid, mz,rf_alpha, val, field, scaling)

   !---------------------------------------------------------------------------
   ! Purpose: Control variable transform through the recursive filter.
   !
   ! Method: 1) Transform to nondimensional v_hat space.
   !         2) Perform recursive filter in non-dimensional space (ds = 1).
   !         3) Scale recursive filter output.
   !         4) Transform filtered field to dimensional space.
   !         5) Optionally scale by background error.
   !---------------------------------------------------------------------------

   implicit none

   type(domain),     intent(inout) :: grid
   integer,          intent(in)    :: mz                             ! Vertical truncation.
   real*8,           intent(in)    :: rf_alpha(mz)                   ! RF scale parameter. 
!  real*8,           intent(in)    :: val(1:jde_int-jds_int+1,mz)                ! Error standard deviation.
   real*8,           intent(in)    :: val(jds_int:jde_int,mz)   ! Error standard deviation.
   real,             intent(inout) :: field(ims_int:ime_int,jms_int:jme_int,kms_int:kme_int) ! Field to be transformed.
      
   integer :: rf_passes_over_two    ! rf_passes / 2
   integer :: i, j, m, n, pass, ij  ! Loop counters.
   real    :: p_x(ims_int:ime_int,jms_int:jme_int)  ! sqrt(Grid box area).
   real*8  :: val_j(grid%xp%jtsy:grid%xp%jtey)
   real*8  :: val_i(grid%xp%itsx:grid%xp%itex)

   logical, optional, intent(in) :: scaling

   !-------------------------------------------------------------------------
   ! [1.0]: Initialise:
   !-------------------------------------------------------------------------  

   if (trace_use_dull) call da_trace_entry("da_transform_through_rf_adj")

   rf_passes_over_two = rf_passes / 2

   ! [1.1] Define inner product (square root of grid box area):
   !$OMP PARALLEL DO &
   !$OMP PRIVATE (ij, i, j)

   do ij = 1 , grid%num_tiles
      do j = grid%j_start(ij), grid%j_end(ij)
         do i = its_int, ite_int
            p_x(i,j) = sqrt(grid%xb%grid_box_area(i,j))
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   !-------------------------------------------------------------------------
   ! [4.0]: Perform 1D recursive filter in y-direction:
   !-------------------------------------------------------------------------
   
   ! [4.3] Optionally scale by background error:
   ! be_s % val = Gridpoint standard deviation - only required for
   ! vert_corr = vert_corr_1 as scaling is performed in vertical transform
   ! for vert_corr = vert_corr_2:

   if (vert_corr == vert_corr_1 .or. (present(scaling))) then
      if (scaling .or. vert_corr == vert_corr_1) then
         do m = 1, mz
            do i = its_int, ite_int
               field(i,jts_int:jte_int,m) = field(i,jts_int:jte_int,m) * val(jts_int:jte_int,m)
            end do
         end do
      end if
   end if

   ! [4.2] Transform filtered field to dimensional space:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (ij ,m, j, i)
   do ij = 1 , grid%num_tiles
      do m = 1, mz
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its_int, ite_int
               grid%xp%v1z(i,j,m) = field(i,j,m) * p_x(i,j)
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   ! [4.1] Apply (i',j',k -> i',j,k') (grid%xp%v1z -> grid%xp%v1y)
   ! convert vertical column to y-stripe

   call da_transpose_z2y (grid)

   !-------------------------------------------------------------------------
   ! [3.0]: Perform 1D recursive filter in y-direction:
   !-------------------------------------------------------------------------

   !  [3.2] Apply 1D filter in y direction:

   n=grid%xp%jtey-grid%xp%jtsy+1
   !$OMP PARALLEL DO &
   !$OMP PRIVATE (m, i, val_j, pass, j)
   do m = grid%xp%ktsy, min(grid%xp%ktey, mz)
      do i = grid%xp%itsy, grid%xp%itey
         do j = grid%xp%jtsy, grid%xp%jtey
            val_j(j) = grid%xp%v1y(i,j,m)
         end do
         do pass = rf_passes_over_two, 1, -1
            call da_recursive_filter_1d_adj(pass, rf_alpha(m), val_j, n)
         end do
         do j = grid%xp%jtsy, grid%xp%jtey
            grid%xp%v1y(i,j,m) = val_j(j)
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   ! [3.1] Apply (i',j,k' -> i,j',k') (grid%xp%v1y -> grid%xp%v1x)
   ! convert from y-stripe to x-stripe

   call da_transpose_y2x (grid)

   !-------------------------------------------------------------------------
   ! [2.0]: Perform 1D recursive filter in x-direction:
   !-------------------------------------------------------------------------

   ! [2.2] Apply 1D filter in x direction:

   n = grid%xp%itex-grid%xp%itsx+1

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( m, j, pass, i, val_i)
   do m = grid%xp%ktsx, min(grid%xp%ktex,mz)
      do j = grid%xp%jtsx, grid%xp%jtex
         do i = grid%xp%itsx, grid%xp%itex
            val_i(i) = grid%xp%v1x(i,j,m)
         end do
         do pass = rf_passes_over_two, 1, -1
            call da_recursive_filter_1d_adj(pass, rf_alpha(m), val_i, n)
         end do
         do i = grid%xp%itsx, grid%xp%itex
            grid%xp%v1x(i,j,m) = val_i(i) 
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   ! [2.1] Apply (i,j',k' -> i',j',k) (grid%xp%v1x -> grid%xp%v1z)
   ! convert from x-stripe to vertical column

   call da_transpose_x2z (grid)

   !-------------------------------------------------------------------------
   ! [1.0]: Initialise:
   !-------------------------------------------------------------------------  

   ! [1.2] Transform to nondimensional v_hat space:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (ij ,m, i, j)
   do ij = 1 , grid%num_tiles
      do m = 1, mz
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its_int, ite_int
               field(i,j,m) = grid%xp%v1z(i,j,m) / p_x(i,j)
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   if (trace_use_dull) call da_trace_exit("da_transform_through_rf_adj")

end subroutine da_transform_through_rf_adj_dual_res



end module da_recursive_filter

