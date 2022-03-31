












module da_dynamics

   
   
   

   use da_control, only : ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte, &
      ids,ide,jds,jde,kds,kde,ips,ipe,jps,jpe,kps,kpe,gamma, gravity,global,test_transforms, &
      fg_format, fg_format_wrf_arw_regional,fg_format_wrf_nmm_regional, &
      fg_format_wrf_arw_global, fg_format_kma_global, balance_geo,balance_geocyc, &
      balance_type, balance_cyc, gravity, convert_fd2uv, trace_use
   use module_domain, only : domain,xb_type
   use module_dm, only : local_communicator, &
      ntasks_x, ntasks_y, data_order_xy, mytask, &
      ntasks
   use module_comm_dm, only : halo_2d_work_sub, halo_wpec_sub, halo_wpec_adj_sub

   use da_define_structures, only : xbx_type
   use da_reporting, only : message, da_error   
   use da_ffts, only : da_solve_poissoneqn_fst, da_solve_poissoneqn_fst_adj
   use da_tracing, only : da_trace_entry, da_trace_exit
   use da_tools, only : da_set_boundary_3d

   implicit none

   contains

subroutine da_balance_cycloterm_adj (rho, ub, vb, u, v, coefx, coefy, term_x, term_y)

   !---------------------------------------------------------------------------
   ! Purpose: Adjoint of da_balance_cycloterm
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(in)    :: rho(ims:ime,jms:jme)    ! Density.
   real, intent(in)    :: ub(ims:ime,jms:jme)     ! Background u wind
   real, intent(in)    :: vb(ims:ime,jms:jme)     ! Background u wind
   real, intent(in)    :: term_x(ims:ime,jms:jme) ! x component of term
   real, intent(in)    :: term_y(ims:ime,jms:jme) ! y component of term
   real, intent(in)    :: coefx(ims:ime,jms:jme)
   real, intent(in)    :: coefy(ims:ime,jms:jme)  ! Mulplicative coeff. 

   real, intent(inout) :: u(ims:ime,jms:jme)      ! u wind increment
   real, intent(inout) :: v(ims:ime,jms:jme)      ! v wind increment

   integer              :: i, j                         ! Loop counters.
   integer              :: is, ie                       ! 1st dim. end points.
   integer              :: js, je                       ! 2nd dim. end points.
   real                 :: data(ims:ime,jms:jme)        ! Work array.

   real                 :: var, varb, uar

   if (trace_use) call da_trace_entry("da_balance_cycloterm_adj")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   ! Computation to check for edge of domain:
   is = its; ie = ite; js = jts; je = jte
   if (.not. global .and. its == ids) is = ids+1
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds) js = jds+1
   if (jte == jde) je = jde-1

   !---------------------------------------------------------------------------
   ! [3.0] Calculate term_y = rho M ( u'dv/dx + v'dv/dy + udv'/dx + vdv'/dy ):
   !---------------------------------------------------------------------------

   ! [3.7] Multiply by rho and add to term_y

   data(its:ite,jts:jte) = rho(its:ite,jts:jte) * term_y(its:ite,jts:jte)

   if (.NOT. global) then      

      ! [3.6] Corner points:

      if (its == ids .AND. jts == jds ) then
         data(its,jts+1) = data(its,jts+1) + 0.5 * data(its,jts)
         data(its+1,jts) = data(its+1,jts) + 0.5 * data(its,jts)
      end if

      if (ite == ide .AND. jts == jds ) then
         data(ite-1,jts) = data(ite-1,jts) + 0.5 * data(ite,jts)
         data(ite,jts+1) = data(ite,jts+1) + 0.5 * data(ite,jts)
      end if

      if (its == ids .AND. jte == jde ) then
         data(its,jde-1) = data(its,jde-1) + 0.5 * data(its,jde)
         data(its+1,jde) = data(its+1,jde) + 0.5 * data(its,jde)
      end if

      if (ite == ide .AND. jte == jde ) then 
         data(ite-1,jte) = data(ite-1,jte) + 0.5 * data(ite,jte)
         data(ite,jte-1) = data(ite,jte-1) + 0.5 * data(ite,jte)
      end if

      ! [3.5] Right boundaries:

      if (jte == jde ) then
         j = jte

         do i = is, ie
            varb = 3.0*vb(i,j)-4.0*vb(i,j-1)+vb(i,j-2)

            var = coefy(i,j)* vb(i,j) * data(i,j)
            uar = coefx(i,j)* data(i,j) * ub(i,j)

            u(i,j)   = u(i,j) + coefx(i,j)*data(i,j) * ( vb(i+1,j) - vb(i-1,j) )
            v(i,j)   = v(i,j) + coefy(i,j)*data(i,j) * varb
            v(i+1,j) = v(i+1,j) + uar                             
            v(i-1,j) = v(i-1,j) - uar                             

            v(i,j) = v(i,j) + 3.0*var
            v(i,j-1) = v(i,j-1) -4.0*var
            v(i,j-2) = v(i,j-2) + var
         end do
      end if

      ! [3.4] Left boundaries:

      if (jts == jds ) then
         j = jts

         do i = is, ie
            varb = -3.0*vb(i,j)+4.0*vb(i,j+1)-vb(i,j+2)

            var = coefy(i,j)*vb(i,j) * data(i,j)
            uar = coefx(i,j)*ub(i,j) * data(i,j)

            v(i,j)   = v(i,j) + coefy(i,j)*data(i,j) * varb
            v(i+1,j) = v(i+1,j) + uar                           
            v(i-1,j) = v(i-1,j) - uar                           

            v(i,j) = v(i,j) - 3.0*var
            v(i,j+1) = v(i,j+1) +4.0*var
            v(i,j+2) = v(i,j+2) - var
         end do
      end if

      ! [3.3] Top boundaries:

      if (ite == ide ) then
         i = ite

         do j = js, je
            varb = 3.0*vb(i,j)-4.0*vb(i-1,j)+vb(i-2,j)

            var = coefx(i,j)* ub(i,j) * data(i,j)
            uar = coefy(i,j)* vb(i,j) * data(i,j)

            u(i,j) = u(i,j) + coefx(i,j)*data(i,j) * varb
            v(i,j) = v(i,j) + coefy(i,j)*data(i,j) * ( vb(i,j+1) - vb(i,j-1) )
            v(i,j+1) = v(i,j+1) + uar                           
            v(i,j-1) = v(i,j-1) - uar                           

            v(i,j) = v(i,j) + 3.0*var
            v(i-1,j) = v(i-1,j) -4.0**var
            v(i-2,j) = v(i-2,j) + var
         end do
      end if

      ! [3.2] Bottom boundaries:

      if (its == ids ) then
         i = its

         do j = js, je
            varb = -3.0*vb(i,j)+4.0*vb(i+1,j)-vb(i+2,j)

            var = coefx(i,j)* ub(i,j) * data(i,j)
            uar = coefy(i,j)* vb(i,j) * data(i,j)

            u(i,j) = u(i,j) + coefx(i,j)*data(i,j) * varb
            v(i,j) = v(i,j) + coefy(i,j)*data(i,j) * ( vb(i,j+1) - vb(i,j-1) )
            v(i,j+1) = v(i,j+1) + uar                             
            v(i,j-1) = v(i,j-1) - uar                             

            v(i,j) = v(i,j) - 3.0*var
            v(i+1,j) = v(i+1,j) +4.0**var
            v(i+2,j) = v(i+2,j) - var
         end do
      end if
   end if ! not global
   
   !  [3.1] Interior points:

   do j = je, js, -1
      do i = ie, is, -1
         uar = coefx(i,j) * ub(i,j) * data(i,j)  
         var = coefy(i,j) * vb(i,j) * data(i,j)  

         u(i,j) = u(i,j) + coefx(i,j)*data(i,j)*( vb(i+1,j) - vb(i-1,j) ) 
         v(i,j) = v(i,j) + coefy(i,j)*data(i,j)*( vb(i,j+1) - vb(i,j-1) ) 
         v(i+1,j) = v(i+1,j) + uar                 
         v(i-1,j) = v(i-1,j) - uar                 
         v(i,j+1) = v(i,j+1) + var                  
         v(i,j-1) = v(i,j-1) - var                  
      end do
   end do
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate term_x = rho M ( u'du/dx + v'du/dy + udu'/dx + vdu'/dy ):
   !---------------------------------------------------------------------------

   ! [2.7] Multiply by rho and add to term_x:

   data(its:ite,jts:jte) = rho(its:ite,jts:jte) * term_x(its:ite,jts:jte)

   if( .NOT. global) then
      ! [2.6] Corner points:

      if (its == ids .AND. jts == jds ) then
         data(its,jts+1) = data(its,jts+1) + 0.5 * data(its,jts)
         data(its+1,jts) = data(its+1,jts) + 0.5 * data(its,jts)
      end if

      if (ite == ide .AND. jts == jds ) then
         data(ite-1,jts) = data(ite-1,jts) + 0.5 * data(ite,jts)
         data(ite,jts+1) = data(ite,jts+1) + 0.5 * data(ite,jts)
      end if

      if (its == ids .AND. jte == jde ) then
         data(its,jde-1) = data(its,jde-1) + 0.5 * data(its,jde)
         data(its+1,jde) = data(its+1,jde) + 0.5 * data(its,jde)
      end if

      if (ite == ide .AND. jte == jde ) then 
         data(ite-1,jte) = data(ite-1,jte) + 0.5 * data(ite,jte)
         data(ite,jte-1) = data(ite,jte-1) + 0.5 * data(ite,jte)
      end if

      ! [2.5] Right boundaries:

      if (jte == jde ) then
         j = jte

         do i = is, ie
            varb = 3.0*ub(i,j)-4.0*ub(i,j-1)+ub(i,j-2)
            var  = coefy(i,j) * vb(i,j) * data(i,j)
            uar  = coefx(i,j) * ub(i,j) * data(i,j)

            u(i+1,j) = u(i+1,j) + uar                   
            u(i-1,j) = u(i-1,j) - uar                   
            v(i,j) = v(i,j) + coefy(i,j)*data(i,j) * varb
            u(i,j) = u(i,j) + coefx(i,j)*data(i,j) * ( ub(i+1,j) - ub(i-1,j) )

            u(i,j) = u(i,j) + 3.0*var
            u(i,j-1) = u(i,j-1) -4.0*var
            u(i,j-2) = u(i,j-2) + var
         end do
      end if

      ! [2.4] Left boundaries:

      if (jts == jds ) then
         j = jts

         do i = is, ie
            varb = -3.0*ub(i,j)+4.0*ub(i,j+1)-ub(i,j+2)
            var = coefy(i,j)*vb(i,j) * data(i,j)
            uar = coefx(i,j)*ub(i,j) * data(i,j)

            u(i+1,j) = u(i+1,j) + uar                 
            u(i-1,j) = u(i-1,j) - uar                 
            v(i,j) = v(i,j) + coefy(i,j)*data(i,j) * varb
            u(i,j) = u(i,j) + coefx(i,j)*data(i,j) * ( ub(i+1,j) - ub(i-1,j) )

            u(i,j) = u(i,j) - 3.0*var
            u(i,j+1) = u(i,j+1) +4.0*var
            u(i,j+2) = u(i,j+2) - var
         end do
      end if

      ! [2.3] Top boundaries:

      if (ite == ide ) then
         i = ite

         do j = js, je
            varb = 3.0*ub(i,j)-4.0*ub(i-1,j)+ub(i-2,j)
            var = coefx(i,j)*ub(i,j) * data(i,j)
            uar = coefy(i,j)*vb(i,j) * data(i,j)

            u(i,j+1) = u(i,j+1) + uar                  
            u(i,j-1) = u(i,j-1) - uar                  
            v(i,j) = v(i,j) + coefy(i,j)*data(i,j) * ( ub(i,j+1) - ub(i,j-1) )
            u(i,j) = u(i,j) + coefx(i,j)*data(i,j) * varb

            u(i,j)   = u(i,j) + 3.0*var
            u(i-1,j) =  u(i-1,j) - 4.0*var
            u(i-2,j) =  u(i-2,j) + var
         end do
      end if

      ! [2.2] Bottom boundaries:

      if (its == ids ) then
         i = its

         do j = js, je
            varb = -3.0*ub(i,j)+4.0*ub(i+1,j)-ub(i+2,j)
            var = coefy(i,j)*ub(i,j) * data(i,j)
            uar = coefy(i,j)*vb(i,j) * data(i,j)

            u(i,j+1) = u(i,j+1) + uar                  
            u(i,j-1) = u(i,j-1) - uar                  
            v(i,j) = v(i,j) + coefy(i,j)*data(i,j) * ( ub(i,j+1) - ub(i,j-1) )
            u(i,j) = u(i,j) + coefx(i,j)*data(i,j) * varb

            u(i,j) = u(i,j) - 3.0*var
            u(i+1,j) =  u(i+1,j) + 4.0*var
            u(i+2,j) =  u(i+2,j) - var
         end do
      end if
   end if ! not global

   ! [2.1] Interior points:

   do j = je, js, -1
      do i = ie, is, -1
         uar = coefx(i,j) * ub(i,j) * data(i,j)
         var = coefy(i,j) * vb(i,j) * data(i,j)

         u(i,j) = u(i,j) + coefx(i,j)*( ub(i+1,j) - ub(i-1,j) ) * data(i,j)
         v(i,j) = v(i,j) + coefy(i,j)*( ub(i,j+1) - ub(i,j-1) ) * data(i,j)
         u(i+1,j) = u(i+1,j) + uar                 
         u(i-1,j) = u(i-1,j) - uar                 
         u(i,j+1) = u(i,j+1) + var                   
         u(i,j-1) = u(i,j-1) - var  
      end do
   end do

   if (trace_use) call da_trace_exit("da_balance_cycloterm_adj")

end subroutine da_balance_cycloterm_adj


subroutine da_balance_cycloterm_lin( rho, ub, vb, u, v, coefx, coefy, term_x, term_y)
                                     
   !---------------------------------------------------------------------------
   ! Purpose: Calculates linearised cyclostrophic term in balance equation.
   !
   ! Method:  Term is rho (u.grad) u on a single level.
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(in)    :: rho(ims:ime,jms:jme)    ! Density.
   real, intent(in)    :: ub(ims:ime,jms:jme)     ! Background u wind
   real, intent(in)    :: vb(ims:ime,jms:jme)     ! Background u wind
   real, intent(in)    :: u(ims:ime,jms:jme)      ! u wind increment
   real, intent(in)    :: v(ims:ime,jms:jme)      ! v wind increment
   real, intent(in)    :: coefx(ims:ime,jms:jme)  ! Multiplicative const.
   real, intent(in)    :: coefy(ims:ime,jms:jme)
   real, intent(inout) :: term_x(ims:ime,jms:jme) ! x component of term.
   real, intent(inout) :: term_y(ims:ime,jms:jme) ! y component of term.

   integer :: i, j                         ! Loop counters.
   integer :: is, ie                       ! 1st dim. end points.
   integer :: js, je                       ! 2nd dim. end points.
   
   real    :: data(ims:ime,jms:jme)        ! Temporary storage.

   real    :: varb, var

   if (trace_use) call da_trace_entry("da_balance_cycloterm_lin")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   ! Computation to check for edge of domain:
   is = its; ie = ite; js = jts; je = jte
   if (.not. global .and. its == ids) is = ids+1
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds) js = jds+1; if (jte == jde) je = jde-1

   !---------------------------------------------------------------------------
   ! [2.0] Calculate term_x = rho M ( u'du/dx + v'du/dy + udu'/dx + vdu'/dy):
   !---------------------------------------------------------------------------

   ! [2.1] Interior points:

   do j = js, je
      do i = is, ie
         data(i,j) = coefx(i,j)*u(i,j) * ( ub(i+1,j) - ub(i-1,j)) + &
                     coefy(i,j)*v(i,j) * ( ub(i,j+1) - ub(i,j-1)) + &
                     coefx(i,j)*ub(i,j) * ( u(i+1,j) - u(i-1,j)) + &
                     coefy(i,j)*vb(i,j) * ( u(i,j+1) - u(i,j-1))
      end do
   end do

   if (.NOT. global) then ! For global only interior points needed  

      ! [2.2] Bottom boundaries:

      if (its == ids) then
         i = its

         do j = js, je 
            var  = -3.0*u (i,j)+4.0*u(i+1,j)-u(i+2,j)
            varb = -3.0*ub(i,j)+4.0*ub(i+1,j)-ub(i+2,j)

            data(i,j) = coefx(i,j)* u(i,j) * varb + &
                        coefy(i,j)* v(i,j) * ( ub(i,j+1) - ub(i,j-1)) + &
                        coefx(i,j)*ub(i,j) * var + &
                        coefy(i,j)*vb(i,j) * ( u(i,j+1) - u(i,j-1))
         end do
      end if

      ! [2.3] Top boundaries:

      if (ite == ide) then
         i = ite

         do j = js, je
            var  = 3.0*u (i,j)-4.0*u (i-1,j)+u (i-2,j)
            varb = 3.0*ub(i,j)-4.0*ub(i-1,j)+ub(i-2,j)

            data(i,j) = coefx(i,j)*u(i,j) * varb + &
                        coefy(i,j)*v(i,j) * ( ub(i,j+1) - ub(i,j-1)) + &
                        coefx(i,j)*ub(i,j) * var + &
                        coefy(i,j)*vb(i,j) * ( u(i,j+1) - u(i,j-1))
         end do
      end if

      ! [2.4] Left boundaries:

      if (jts == jds) then
         j = jts

         do i = is, ie
            var  = -3.0*u (i,j)+4.0*u (i,j+1)-u (i,j+2)
            varb = -3.0*ub(i,j)+4.0*ub(i,j+1)-ub(i,j+2)

            data(i,j) = coefx(i,j)*u(i,j) * ( ub(i+1,j) - ub(i-1,j)) + &
                        coefy(i,j)*v(i,j) * varb + &
                        coefx(i,j)*ub(i,j) * ( u(i+1,j) - u(i-1,j)) + &
                        coefy(i,j)*vb(i,j) * var
         end do
      end if

      ! [2.5] Right boundaries:

      if (jte == jde) then
         j = jte

         do i = is, ie
            var  = 3.0*u (i,j)-4.0*u (i,j-1)+u (i,j-2)
            varb = 3.0*ub(i,j)-4.0*ub(i,j-1)+ub(i,j-2)

            data(i,j) = coefx(i,j)*u(i,j) * ( ub(i+1,j) - ub(i-1,j)) + &
                        coefy(i,j)*v(i,j) * varb + &
                        coefx(i,j)*ub(i,j) * ( u(i+1,j) - u(i-1,j)) + &
                        coefy(i,j)*vb(i,j) * var
         end do
      end if

      ! [2.6] Corner points:

      if (its == ids .AND. jts == jds) then
         data(its,jts) = 0.5 * ( data(its,jts+1) + data(its+1,jts))
      end if

      if (ite == ide .AND. jts == jds) then
         data(ite,jts) = 0.5 * ( data(ite-1,jts) + data(ite,jts+1))
      end if

      if (its == ids .AND. jte == jde) then
         data(its,jde) = 0.5 * ( data(its,jde-1) + data(its+1,jde))
      end if

      if (ite == ide .AND. jte == jde) then 
         data(ite,jte) = 0.5 * ( data(ite-1,jte) + data(ite,jte-1))
      end if
   end if ! not global
      
   ! [2.7] Multiply by rho and add to term_x:

   term_x(its:ite,jts:jte) = rho(its:ite,jts:jte) * data(its:ite,jts:jte) + &
                             term_x(its:ite,jts:jte)

   !---------------------------------------------------------------------------
   ! [3.0] Calculate term_y = rho M ( u'dv/dx + v'dv/dy + udv'/dx + vdv'/dy):
   !---------------------------------------------------------------------------

   ! [3.1] Interior points:

   do j = js, je
      do i = is, ie
         data(i,j) = coefx(i,j)*u(i,j) * ( vb(i+1,j) - vb(i-1,j)) + &
                     coefy(i,j)*v(i,j) * ( vb(i,j+1) - vb(i,j-1)) + &
                     coefx(i,j)*ub(i,j) * ( v(i+1,j) - v(i-1,j)) + &
                     coefy(i,j)*vb(i,j) * ( v(i,j+1) - v(i,j-1))
      end do
   end do
   
   if (.NOT. global) then      ! For global only interior points needed  
      ! [3.2] Bottom boundaries:

      if (its == ids) then
         i = its

         do j = js, je
            var  = -3.0*v (i,j)+4.0*v (i+1,j)-v (i+2,j)
            varb = -3.0*vb(i,j)+4.0*vb(i+1,j)-vb(i+2,j)

            data(i,j) = coefx(i,j)*u(i,j) * varb + &
                        coefy(i,j)*v(i,j) * ( vb(i,j+1) - vb(i,j-1)) + &
                        coefx(i,j)*ub(i,j) * var + &
                        coefy(i,j)*vb(i,j) * ( v(i,j+1) - v(i,j-1))
         end do
      end if

      !  [3.3] Top boundaries:

      if (ite == ide) then
         i = ite

         do j = js, je
            var  = 3.0*v (i,j)-4.0*v (i-1,j)+v (i-2,j)
            varb = 3.0*vb(i,j)-4.0*vb(i-1,j)+vb(i-2,j)

            data(i,j) = coefx(i,j)*u(i,j) * varb + &
                        coefy(i,j)*v(i,j) * ( vb(i,j+1) - vb(i,j-1)) + &
                        coefx(i,j)*ub(i,j) * var + &
                        coefy(i,j)*vb(i,j) * ( v(i,j+1) - v(i,j-1))
         end do
      end if

      !  [3.4] Left boundaries:

      if (jts == jds) then
         j = jts

         do i = is, ie
            varb = -3.0*vb(i,j)+4.0*vb(i,j+1)-vb(i,j+2)
            var  = -3.0*v (i,j)+4.0*v (i,j+1)-v (i,j+2)

            data(i,j) = coefx(i,j)*u(i,j) * ( vb(i+1,j) - vb(i-1,j)) + &
                        coefy(i,j)*v(i,j) * varb + &
                        coefx(i,j)*ub(i,j) * ( v(i+1,j) - v(i-1,j)) + &
                        coefy(i,j)*vb(i,j) * var
         end do
      end if

      !  [3.5] Right boundaries:

      if (jte == jde) then
         j = jte

         do i = is, ie
            varb = 3.0*vb(i,j)-4.0*vb(i,j-1)+vb(i,j-2)
            var  = 3.0*v (i,j)-4.0*v (i,j-1)+v (i,j-2)

            data(i,j) = coefx(i,j)*u(i,j) * ( vb(i+1,j) - vb(i-1,j)) + &
                        coefy(i,j)*v(i,j) * varb + &
                        coefx(i,j)*ub(i,j) * ( v(i+1,j) - v(i-1,j)) + &
                        coefy(i,j)*vb(i,j) * var
         end do
      end if

      !  [3.6] Corner points:

      if (its == ids .AND. jts == jds) then
         data(its,jts) = 0.5 * ( data(its,jts+1) + data(its+1,jts))
      end if

      if (ite == ide .AND. jts == jds) then
         data(ite,jts) = 0.5 * ( data(ite-1,jts) + data(ite,jts+1))
      end if

      if (its == ids .AND. jte == jde) then
         data(its,jde) = 0.5 * ( data(its,jde-1) + data(its+1,jde))
      end if

      if (ite == ide .AND. jte == jde) then 
         data(ite,jte) = 0.5 * ( data(ite-1,jte) + data(ite,jte-1))
      end if
   end if ! not global

   ! [3.7] Multiply by  rho and add to term_y

   term_y(its:ite,jts:jte) = rho(its:ite,jts:jte) * data(its:ite,jts:jte) + term_y(its:ite,jts:jte)

   if (trace_use) call da_trace_exit("da_balance_cycloterm_lin")

end subroutine da_balance_cycloterm_lin


subroutine da_balance_cycloterm (xb, k, term_x, term_y)

   !---------------------------------------------------------------------------
   !  Purpose: Calculates cyclostrophic term in balance equation.
   !
   !  Method:  Term is rho (u.grad) u on a single level.
   !---------------------------------------------------------------------------

   implicit none
   
   type(xb_type), intent(in)    :: xb           ! First guess structure.
   integer,       intent(in)    :: k            ! Model level.
   real,          intent(inout) :: term_x(:,:)  ! x component of term.
   real,          intent(inout) :: term_y(:,:)  ! y component of term.

   integer :: i, j                         ! Loop counters.
   integer :: is, ie                       ! 1st dim. end points.
   integer :: js, je                       ! 2nd dim. end points.
   
   real    :: data(ims:ime,jms:jme)        ! Temporary storage.

   real    :: varb

   if (trace_use) call da_trace_entry("da_balance_cycloterm")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------
   
   ! Computation to check for edge of domain:
   is = its; ie = ite; js = jts; je = jte
   if (.not. global .and. its==ids) is = ids+1
   if (.not. global .and. ite==ide) ie = ide-1
   if (jts==jds) js = jds+1
   if (jte==jde) je = jde-1
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate term_x = rho M (u.du/dx + v.du/dy):
   !---------------------------------------------------------------------------

   ! [2.1] Interior points:

   do j = js, je
      do i = is, ie
         data(i,j) = xb%u(i,j,k) * xb%coefx(i,j)*(xb%u(i+1,j,k) - &
            xb%u(i-1,j,k)) + xb%v(i,j,k) * xb%coefy(i,j)*(xb%u(i,j+1,k) - &
            xb%u(i,j-1,k))
      end do
   end do
   
   if (.NOT. global) then ! For global only interior points

      ! [2.2] Bottom boundaries:

      if (its==ids) then
         i = its

         do j = js, je 
            varb = -3.0*xb%u(i,j,k)+4.0*xb%u(i+1,j,k)-xb%u(i+2,j,k)

            data(i,j) = xb%u(i,j,k) * xb%coefx(i,j) * varb + &
               xb%v(i,j,k) * xb%coefy(i,j) * (xb%u(i,j+1,k) - xb%u(i,j-1,k))
         end do
      end if

      ! [2.3] Top boundaries:

      if (ite==ide) then
         i = ite

         do j = js, je
            varb = 3.0*xb%u(i,j,k)-4.0*xb%u(i-1,j,k)+xb%u(i-2,j,k)

            data(i,j) = xb%u(i,j,k) * xb%coefx(i,j) * varb + &
               xb%v(i,j,k) * xb%coefy(i,j) * (xb%u(i,j+1,k) - xb%u(i,j-1,k))
         end do
      end if

      ! [2.4] Left boundaries:

      if (jts==jds) then
         j = jts

         do i = is, ie
            varb = -3.0*xb%u(i,j,k)+4.0*xb%u(i,j+1,k)-xb%u(i,j+2,k)

            data(i,j) = xb%u(i,j,k) * xb%coefx(i,j) * (xb%u(i+1,j,k) - &
               xb%u(i-1,j,k)) + xb%v(i,j,k) * xb%coefy(i,j) * varb
         end do
      end if

      ! [2.5] Right boundaries:

      if (jte==jde) then
         j = jte

         do i = is, ie
            varb = 3.0*xb%u(i,j,k)-4.0*xb%u(i,j-1,k)+xb%u(i,j-2,k)

            data(i,j) = xb%u(i,j,k) * xb%coefx(i,j) * (xb%u(i+1,j,k) - &
               xb%u(i-1,j,k)) + xb%v(i,j,k) * xb%coefy(i,j) * varb
         end do
      end if

      ! [2.6] Corner points:

      if (its==ids .AND. jts==jds) then
         data(its,jts) = 0.5 * (data(its,jts+1) + data(its+1,jts))
      end if

      if (ite==ide .AND. jts==jds) then
         data(ite,jts) = 0.5 * (data(ite-1,jts) + data(ite,jts+1))
      end if

      if (its==ids .AND. jte==jde) then
         data(its,jde) = 0.5 * (data(its,jde-1) + data(its+1,jde))
      end if

      if (ite==ide .AND. jte==jde) then 
         data(ite,jte) = 0.5 * (data(ite-1,jte) + data(ite,jte-1))
      end if
   end if ! not global

   !  [2.7] Multiply by rho  and add to term_x:
      
   term_x(its:ite,jts:jte) = xb%rho(its:ite,jts:jte,k)*data(its:ite,jts:jte) + term_x(its:ite,jts:jte)

   !---------------------------------------------------------------------------
   ! [3.0] Calculate term_y = rho M (u.dv/dx + v.dv/dy):
   !---------------------------------------------------------------------------

   ! [3.1] Interior points:

   do j = js, je
      do i = is, ie
         data(i,j) = xb%u(i,j,k) * xb%coefx(i,j)*(xb%v(i+1,j,k) - xb%v(i-1,j,k)) + &
                     xb%v(i,j,k) * xb%coefy(i,j)*(xb%v(i,j+1,k) - xb%v(i,j-1,k))
      end do
   end do
   
   
   if (.NOT. global) then ! For global only interior points

      ! [3.2] Bottom boundaries:

      if (its==ids) then
         i = its

         do j = js, je
            varb = -3.0*xb%v(i,j,k)+4.0*xb%v(i+1,j,k)-xb%v(i+2,j,k)

            data(i,j) = xb%u(i,j,k) * xb%coefx(i,j)* varb + &
                        xb%v(i,j,k) * xb%coefy(i,j)* (xb%v(i,j+1,k) - xb%v(i,j-1,k))
         end do
      end if

      !  [3.3] Top boundaries:

      if (ite==ide) then
         i = ite

         do j = js, je
            varb = 3.0*xb%v(i,j,k)-4.0*xb%v(i-1,j,k)+xb%v(i-2,j,k)

            data(i,j) = xb%u(i,j,k) * xb%coefx(i,j)* varb + &
                        xb%v(i,j,k) * xb%coefy(i,j)* (xb%v(i,j+1,k) - xb%v(i,j-1,k))
         end do
      end if

      ! [3.4] Left boundaries:

      if (jts==jds) then
         j = jts

         do i = is, ie
            varb = -3.0*xb%v(i,j,k)+4.0*xb%v(i,j+1,k)-xb%v(i,j+2,k)

            data(i,j) = xb%u(i,j,k) * xb%coefx(i,j)* (xb%v(i+1,j,k) - xb%v(i-1,j,k)) + &
                        xb%v(i,j,k) * xb%coefy(i,j)* varb
         end do
      end if

      ! [3.5] Right boundaries:

      if (jte==jde) then
         j = jte

         do i = is, ie
            varb = 3.0*xb%v(i,j,k)-4.0*xb%v(i,j+1,k)+xb%v(i,j+2,k)

            data(i,j) = xb%u(i,j,k) * xb%coefx(i,j)* (xb%v(i+1,j,k) - xb%v(i-1,j,k)) + &
                        xb%v(i,j,k) * xb%coefy(i,j)* varb
         end do
      end if

      ! [3.6] Corner points:

      if (its==ids .AND. jts==jds) then
         data(its,jts) = 0.5 * (data(its,jts+1) + data(its+1,jts))
      end if

      if (ite==ide .AND. jts==jds) then
         data(ite,jts) = 0.5 * (data(ite-1,jts) + data(ite,jts+1))
      end if

      if (its==ids .AND. jte==jde) then
         data(its,jde) = 0.5 * (data(its,jde-1) + data(its+1,jde))
      end if

      if (ite==ide .AND. jte==jde) then 
         data(ite,jte) = 0.5 * (data(ite-1,jte) + data(ite,jte-1))
      end if
   end if ! not global

   ! [3.7] Multiply by rho and add to term_y

   term_y(its:ite,jts:jte)=xb%rho(its:ite,jts:jte,k)* data(its:ite,jts:jte) + &
      term_y(its:ite,jts:jte)

   if (trace_use) call da_trace_exit("da_balance_cycloterm")

end subroutine da_balance_cycloterm


subroutine da_balance_equation_adj(grid, xbx, phi_b, u, v)

   !---------------------------------------------------------------------------
   ! Purpose: Adjoint of da_balance_equation
   !---------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)               :: grid

   type (xbx_type),intent(in) :: xbx              ! Header & non-gridded vars.
   real, intent(in)    :: phi_b(ims:ime,jms:jme,kms:kme) ! Balanced mass increment.
   real, intent(inout) :: u(ims:ime,jms:jme,kms:kme) ! u wind comp. (dot pts)
   real, intent(inout) :: v(ims:ime,jms:jme,kms:kme) ! v wind comp. (dot pts)

   integer             :: i, j, k          ! Loop counters.
   integer             :: is, ie           ! 1st dim. end points.
   integer             :: js, je           ! 2nd dim. end points.

   real, dimension(ims:ime,jms:jme) :: coefx, &   ! Multiplicative coefficient.
                                       coefy, &   ! Multiplicative coefficient.
                                       term_x, &  ! Balance eqn x term
                                       term_y     ! Balance eqn y term

   real, dimension(ims:ime,jms:jme,kms:kme) :: del2phi_b  ! Del**2 phi_b/M**2
   real                       :: coeff1, coeff2   ! Multiplicative coefficient.

   if (trace_use) call da_trace_entry("da_balance_equation_adj")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   ! Computation to check for edge of domain:
   is = its-1; ie = ite+1; js = jts-1; je = jte+1
   if (.not. global .and. its == ids) is = ids+1
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds ) js = jds+1
   if (jte == jde ) je = jde-1

   if (fg_format == fg_format_kma_global) then
      coefx(is:ie,js:je) = grid%xb%coefx(is:ie,js:je)
      coefy(is:ie,js:je) = grid%xb%coefy(is:ie,js:je)
   else if (fg_format == fg_format_wrf_arw_global) then
      write (unit=message(1),fmt='(A,I3)') ' needs work for fg_format  = ',fg_format
      call da_error("da_balance_equation_adj.inc",46,message(1:1))
   else if (fg_format == fg_format_wrf_arw_regional) then
      coefx(is:ie,js:je) = grid%xb%coefz(is:ie,js:je)
      coefy(is:ie,js:je) = coefx(is:ie,js:je)
   else if (fg_format == fg_format_wrf_nmm_regional) then
      write (unit=message(1),fmt='(A,I3)') ' needs work for fg_format  = ',fg_format
      call da_error("da_balance_equation_adj.inc",52,message(1:1))
   else
      write (unit=message(1),fmt='(A,I3)') ' Wrong FG_FORMAT = ',fg_format
      call da_error("da_balance_equation_adj.inc",55,message(1:1))
   end if

   ! [1.1] Multiplicative coefficent for conversion RHS->Del**2 phi_b/M**2:

   del2phi_b(:,:,:) = 0.0

   !---------------------------------------------------------------------------
   ! [3.0] Solve Del**2 phi_b = RHS for phi_b:
   !---------------------------------------------------------------------------

   call da_solve_poissoneqn_fst_adj(grid,xbx, phi_b, del2phi_b)

   !---------------------------------------------------------------------------
   ! [2.0] Calculate RHS of balance equation in gridpt space:
   !---------------------------------------------------------------------------

   do k = kts, kte
      ! [2.4] Del**2 Phi_b boundary conditions (null as zero boundary conditions):

      ! [2.3] Take divergence to get Del**2 phi_b/M**2:

      term_x(ims:ime,jms:jme) = 0.0
      term_y(ims:ime,jms:jme) = 0.0

      do j = je, js, -1
         do i = ie, is, -1
            coeff1 = coefx(i,j) * del2phi_b(i,j,k)
            coeff2 = coefy(i,j) * del2phi_b(i,j,k)

            term_x(i+1,j) = term_x(i+1,j) - coeff1
            term_x(i-1,j) = term_x(i-1,j) + coeff1
            term_y(i,j+1) = term_y(i,j+1) - coeff2
            term_y(i,j-1) = term_y(i,j-1) + coeff2
         end do
      end do

      ! [2.2] Include cyclostrophic terms in balance eqn if requested:

      if (balance_type == balance_cyc .OR. balance_type == balance_geocyc ) then
         call da_balance_cycloterm_adj (grid%xb%rho(:,:,k),grid%xb%u(:,:,k),&
            grid%xb%v(:,:,k), u(:,:,k), v(:,:,k), grid%xb%coefx(:,:), grid%xb%coefy(:,:),&
            term_x(:,:), term_y(:,:))
      end if

      
      ! [2.1] Calculate geostrophic terms in balance eqn:
 
      if (balance_type == balance_geo .OR. balance_type == balance_geocyc ) then
         call da_balance_geoterm_adj (grid%xb%cori, grid%xb%rho(:,:,k), term_x, term_y, &
            u(:,:,k), v(:,:,k))
      end if
   end do

   if (trace_use) call da_trace_exit("da_balance_equation_adj")

end subroutine da_balance_equation_adj


subroutine da_balance_equation_lin(grid, xbx, u, v, phi_b)

   !---------------------------------------------------------------------------
   !  Purpose: Calculates balanced mass increment phi_b from wind increments.
   !
   !  Method:  Solve grad**2 Phi_b = - div[ k x rho f u + rho (u.grad ) u ] by
   !           1) Calculate RHS of balance equation in gridpt space.
   !           2) Solve Del**2 phi_b = RHS for phi_b using double (FCT).
   !
   !---------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid
   type(xbx_type), intent(in)    :: xbx                            ! Header & non-gridded vars. 
   real,           intent(in)    :: u(ims:ime,jms:jme,kms:kme)     ! u wind comp.
   real,           intent(in)    :: v(ims:ime,jms:jme,kms:kme)     ! v wind comp.
   real,           intent(out)   :: phi_b(ims:ime,jms:jme,kms:kme) ! Balanced mass increment.

   integer :: i, j, k                 ! Loop counters.
   integer :: is, ie                  ! 1st dim. end points.
   integer :: js, je                  ! 2nd dim. end points.

   real    :: coefx(ims:ime,jms:jme)  ! Multiplicative coefficient.
   real    :: coefy(ims:ime,jms:jme)  ! Multiplicative coefficient.
   real    :: term_x(ims:ime,jms:jme) ! Balance eqn x term
   real    :: term_y(ims:ime,jms:jme) ! Balance eqn y term
   
   real    :: del2phi_b(ims:ime,jms:jme,kms:kme)  ! Del**2 phi_b/M**2

   if (trace_use) call da_trace_entry("da_balance_equation_lin")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise iand set Multipilactive constants
   !---------------------------------------------------------------------------

   ! Computation to check for edge of domain:

   is = its; ie = ite; js = jts; je = jte
   if (.not.global .and. its == ids) is = ids+1
   if (.not.global .and. ite == ide) ie = ide-1
   if (jts == jds ) js = jds+1; if (jte == jde) je = jde-1

   if (fg_format == fg_format_kma_global) then
      coefx = grid%xb%coefx
      coefy = grid%xb%coefy
   else if( fg_format == fg_format_wrf_arw_regional) then
      coefx = grid%xb%coefz
      coefy = coefx   
   else if (fg_format == fg_format_wrf_arw_global) then
      write (unit=message(1),fmt='(A,I3)') ' needs work for fg_format_wrf_arw_global  = ',fg_format
      call da_error("da_balance_equation_lin.inc",52,message(1:1))
   else if (fg_format == fg_format_wrf_nmm_regional) then
      write (unit=message(1),fmt='(A,I3)') ' needs work for fg_format_wrf_nmm_regional = ',fg_format
      call da_error("da_balance_equation_lin.inc",55,message(1:1))
   else
      write(unit=message(1),fmt='(A,I3)') 'Wrong FG_FORMAT = ',fg_format
      call da_error("da_balance_equation_lin.inc",58,message(1:1))
   end if

   ! [1.1] Multiplicative coefficent for conversion RHS->Del**2 phi_b/M**2:

   do k = kts, kte
      term_x(ims:ime,jms:jme) = 0.0
      term_y(ims:ime,jms:jme) = 0.0

      !---------------------------------------------------------------------------
      ! [2.0] Calculate RHS of balance equation in gridpt space:
      !---------------------------------------------------------------------------

      ! [2.1] Include geostrophic terms in balance eqn if requested:

      if (balance_type == balance_geo .OR. balance_type == balance_geocyc ) then
         call da_balance_geoterm_lin (grid%xb % cori, grid%xb % rho(:,:,k), u(:,:,k), v(:,:,k), &
            term_x, term_y)
      end if
      
      ! [2.2] Include cyclostrophic terms in balance eqn if requested:

      if (balance_type == balance_cyc .OR. balance_type == balance_geocyc ) then
         call da_balance_cycloterm_lin (grid%xb % rho(:,:,k), grid%xb % u(:,:,k),   &
            grid%xb % v(:,:,k), u(:,:,k), v(:,:,k), grid%xb % coefx(:,:), grid%xb % coefy(:,:), &
            term_x(:,:), term_y(:,:))
      end if
      
      ! [2.3] Take divergence to get Del**2 phi_b/M**2:

      do j = js, je
         do i = is, ie
          del2phi_b(i,j,k) = -coefx(i,j)*( term_x(i+1,j) - term_x(i-1,j)) &
                             -coefy(i,j)*( term_y(i,j+1) - term_y(i,j-1)) 
         end do
      end do

      ! [2.4] Del**2 Phi_b  boundary conditions:

      if (.not. global .and. its == ids ) del2phi_b(its,js:je,k) = 0.0
      if (.not. global .and. ite == ide ) del2phi_b(ite,js:je,k) = 0.0
      if (jts == jds ) del2phi_b(is:ie,jts,k) = 0.0
      if (jte == jde ) del2phi_b(is:ie,jte,k) = 0.0

      if (.not. global .and. (its == ids .and. jts == jds) ) del2phi_b(its,jts,k) = 0.0
      if (.not. global .and. (its == ids .and. jte == jde) ) del2phi_b(its,jte,k) = 0.0
      if (.not. global .and. (ite == ide .and. jts == jds) ) del2phi_b(ite,jts,k) = 0.0
      if (.not. global .and. (ite == ide .and. jte == jde) ) del2phi_b(ite,jte,k) = 0.0     
   end do

   !------------------------------------------------------------------------------
   !  [3.0] Solve Del**2 phi_b = RHS for phi_b:
   !------------------------------------------------------------------------------

   call da_solve_poissoneqn_fst(grid, xbx, del2phi_b, phi_b) 

   if (trace_use) call da_trace_exit("da_balance_equation_lin")

end subroutine da_balance_equation_lin


subroutine da_balance_geoterm_adj( cori, rho, term_x, term_y, u, v)
 
   !---------------------------------------------------------------------------
   ! Purpose: Adjoint of da_balance_geoterm.
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(in)    :: cori(ims:ime,jms:jme)   ! Coriolis factor.
   real, intent(in)    :: rho(ims:ime,jms:jme)    ! Density
   real, intent(in)    :: term_x(ims:ime,jms:jme) ! x component of term.
   real, intent(in)    :: term_y(ims:ime,jms:jme) ! y component of term.
   real, intent(inout) :: u(ims:ime,jms:jme)      ! u wind increment
   real, intent(inout) :: v(ims:ime,jms:jme)      ! v wind increment

   if (trace_use) call da_trace_entry("da_balance_geoterm_adj")

   !---------------------------------------------------------------------------
   ! [2.0] Calculate term_y = f rho u~:
   !---------------------------------------------------------------------------

   u(its:ite,jts:jte) = u(its:ite,jts:jte) + rho(its:ite,jts:jte) * cori(its:ite,jts:jte) &
      * term_y(its:ite,jts:jte)

   !---------------------------------------------------------------------------
   ! [1.0] Calculate term_x = -f rho v~:
   !---------------------------------------------------------------------------

   v(its:ite,jts:jte) = v(its:ite,jts:jte) - rho(its:ite,jts:jte) * cori(its:ite,jts:jte) &
      * term_x(its:ite,jts:jte)

   if (trace_use) call da_trace_exit("da_balance_geoterm_adj")

end subroutine da_balance_geoterm_adj


subroutine da_balance_geoterm_lin( cori, rho, u, v, term_x, term_y)

   !---------------------------------------------------------------------------
   !  Purpose: calculates linearised geostrophic term in balance equation.
   !
   !  method:  term is k x rho f u on a single level.
   !
   !  assumptions: Various (see documentation).
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(in)    :: cori(ims:ime,jms:jme)   ! Coriolis factor.
   real, intent(in)    :: rho(ims:ime,jms:jme)    ! Density
   real, intent(in)    :: u(ims:ime,jms:jme)      ! u wind increment
   real, intent(in)    :: v(ims:ime,jms:jme)      ! v wind increment
   real, intent(inout) :: term_x(ims:ime,jms:jme) ! x component of term.
   real, intent(inout) :: term_y(ims:ime,jms:jme) ! y component of term.

   integer :: is, js              ! i,j lower loop limits
   integer :: ie, je              ! i,j upper loop limits

   if (trace_use) call da_trace_entry("da_balance_geoterm_lin")

   ! Set loop limits

   is = its-1
   ie = ite+1
   js = jts-1
   je = jte+1
   if (its == ids ) is = ids
   if (ite == ide ) ie = ide
   if (jts == jds ) js = jds
   if (jte == jde ) je = jde

   !---------------------------------------------------------------------------
   ! [1.0] Calculate term_x = -f rho v~:
   !---------------------------------------------------------------------------

   term_x(is:ie,js:je) = term_x(is:ie,js:je) - rho(is:ie,js:je) * cori(is:ie,js:je) * v(is:ie,js:je)

   !---------------------------------------------------------------------------
   ! [2.0] Calculate term_y = f rho u~:
   !---------------------------------------------------------------------------

   term_y(is:ie,js:je) = term_y(is:ie,js:je) + rho(is:ie,js:je) * cori(is:ie,js:je) * u(is:ie,js:je)

   if (trace_use) call da_trace_exit("da_balance_geoterm_lin")

end subroutine da_balance_geoterm_lin


subroutine da_hydrostaticp_to_rho_adj (xb, rho, p) 

   !---------------------------------------------------------------------------
   !  Purpose: Adjoint of da_hydrostaticp_to_rho.
   !
   !  Method:  Standard adjoint coding.
   !
   !  Assumptions: 1) Hydrostatic pressure.
   !               2) Model level stored top down.
   !---------------------------------------------------------------------------

   implicit none
   
   type(xb_type), intent(in)    :: xb                           ! First guess structure.
   real,          intent(in)    :: rho(ims:ime,jms:jme,kms:kme) ! Density inc. (cross pts).
   real,          intent(inout) :: p(ims:ime,jms:jme,kms:kme)   ! Pressure inc. (cross pts)

   integer               :: i, j, k      ! Loop counters.
   real                  :: delta1       ! Height difference.
   real                  :: delta2       ! Height difference.
   real                  :: dPdz         ! Vertical pressure gradient.
   real                  :: temp(its:ite,jts:jte) ! Temporary array.

   if (trace_use) call da_trace_entry("da_hydrostaticp_to_rho_adj")
   
   !---------------------------------------------------------------------------
   ! [4.0] Calculate density increment on top level:
   !---------------------------------------------------------------------------

   do j = jts, jte
      do i = its, ite
      
         ! Put together to get density increment at bottom level:
         dPdz = -rho(i,j,kte) / gravity      
      
         ! dP~/dz by forwards one-sided 2nd order finite differences:
         
         delta1 = xb % h(i,j,kte) - xb % h(i,j,kte-1)
         delta2 = xb % h(i,j,kte) - xb % h(i,j,kte-2)
         
         p(i,j,k) = p(i,j,kte) + ( delta1 + delta2 ) * dPdz / (delta1 * delta2)
         p(i,j,k-1) = p(i,j,kte-1) - (delta2 / delta1) * dPdz / (delta2 - delta1)
         p(i,j,k-2) = p(i,j,kte-2) + (delta1 / delta2) * dPdz / (delta2 - delta1)
      end do
   end do

   !---------------------------------------------------------------------------
   ! [3.0] Calculate density increment on top level:
   !---------------------------------------------------------------------------

   do j = jts, jte
      do i = its, ite

         ! Put together to get density increment at top level:
         dPdz = -rho(i,j,kts) / gravity

         ! dP~/dz by backwards one-sided 2nd order finite differences:

         delta1 = xb % h(i,j,kts+1) - xb % h(i,j,kts)
         delta2 = xb % h(i,j,kts+2) - xb % h(i,j,kts)

         p(i,j,kts)   = p(i,j,k)     - (delta1 + delta2) * dPdz / (delta1 * delta2)
         p(i,j,kts+1) = p(i,j,kts+1) + (delta2 / delta1) * dPdz / (delta2 - delta1)
         p(i,j,kts+2) = p(i,j,kts+2) - (delta1 / delta2) * dPdz / (delta2 - delta1)
      end do
   end do  
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate density increments at all levels except top/bottom:
   !---------------------------------------------------------------------------

   do k = kte-1, kts+1, -1 
      temp(its:ite,jts:jte) = -rho(its:ite,jts:jte,k) / &
         ((xb%h(its:ite,jts:jte,k+1) - xb%h(its:ite,jts:jte,k-1)) * gravity)

      p(its:ite,jts:jte,k-1) = p(its:ite,jts:jte,k-1) - temp(its:ite,jts:jte)
      p(its:ite,jts:jte,k+1) = p(its:ite,jts:jte,k+1) + temp(its:ite,jts:jte)
   end do                                  

end subroutine da_hydrostaticp_to_rho_adj


subroutine da_hydrostaticp_to_rho_lin(grid, p, rho ) 

   !---------------------------------------------------------------------------
   !  Purpose: Calculates density increment from pressure increment fiteld.
   !
   !  Method:  Hydrostatic eqn => rho~ = -dP~/dz / g.
   !
   !  Assumptions: 1) Hydrostatic pressure.
   !               2) Model level stored top down.
   !
   !---------------------------------------------------------------------------

   implicit none
   
   type(domain), intent(in)  :: grid
   real,         intent(in)  :: p(ims:ime,jms:jme,kms:kme)   ! Pressure inc. (cross pts)
   real,         intent(out) :: rho(ims:ime,jms:jme,kms:kme) ! Density inc. (cross pts)

   integer :: i, j, k      ! Loop counters.
   real    :: delta1       ! Height difference.
   real    :: delta2       ! Height difference.
   real    :: dPdz         ! Vertical pressure gradient.

   if (trace_use) call da_trace_entry("da_hydrostaticp_to_rho_lin")

   !---------------------------------------------------------------------------
   ! [2.0] Calculate density increments at all levels except top/bottom:
   !---------------------------------------------------------------------------

   do k = kts+1, kte-1 
      rho(its:ite,jts:jte,k) =  -( p(its:ite,jts:jte,k+1) - p(its:ite,jts:jte,k-1) ) / &
                             ( ( grid%xb % h(its:ite,jts:jte,k+1) - &
                             grid%xb % h(its:ite,jts:jte,k-1) ) * gravity )
   end do                                  

   !---------------------------------------------------------------------------
   ! [3.0] Calculate density increment on bottom level:
   !---------------------------------------------------------------------------

   k = kts
   do j = jts, jte
      do i = its, ite
         ! dP~/dz by backwards one-sided 2nd order finite differences:
         
         delta1 = grid%xb % h(i,j,k+1) - grid%xb % h(i,j,k)
         delta2 = grid%xb % h(i,j,k+2) - grid%xb % h(i,j,k)
         dPdz = -( delta1 + delta2 ) * p(i,j,k)  / ( delta1 * delta2 ) + &
                 ( delta2 / delta1 ) * p(i,j,k+1)/ ( delta2 - delta1 ) - &
                 ( delta1 / delta2 ) * p(i,j,k+2)/ ( delta2 - delta1 )
                      
         ! Put together to get density increment at top level:
         rho(i,j,k) = -dPdz / gravity
      end do
   end do
                       
   !---------------------------------------------------------------------------
   ! [4.0] Calculate density increment on top level:
   !---------------------------------------------------------------------------

   k = kte
   do j = jts, jte
      do i = its, ite
         ! dP~/dz by forwards one-sided 2nd order finite differences:
         
         delta1 = grid%xb % h(i,j,k) - grid%xb % h(i,j,k-1)
         delta2 = grid%xb % h(i,j,k) - grid%xb % h(i,j,k-2)
         
         dPdz = ( delta1 + delta2 ) * p(i,j,k)   / ( delta1 * delta2 ) - &
                ( delta2 / delta1 ) * p(i,j,k-1) / ( delta2 - delta1 ) + &
                ( delta1 / delta2 ) * p(i,j,k-2) / ( delta2 - delta1 )

         ! Put together to get density increment at bottom level:
         rho(i,j,k) = -dPdz / gravity
      end do
   end do

   if (trace_use) call da_trace_exit("da_hydrostaticp_to_rho_lin")
   
end subroutine da_hydrostaticp_to_rho_lin


subroutine da_psichi_to_uv(psi, chi, coefx,coefy, u, v)

   !---------------------------------------------------------------------------
   !  Purpose: Calculate wind components u and v from psi and chi.
   !
   !  Method:  u = coefx * (-dpsi/dy + dchi/dx)
   !           v = coefy * ( dpsi/dx + dchi/dy)
   !
   !  Assumptions: Unstaggered grid.
   !               Lateral boundary conditions - dpsi/dn, dchi/dn = 0 (FCT)
   !               grid size may or may not be equal
   !
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !----------------------------------------------------------------------------

   implicit none
   
   real, intent(inout) :: psi(ims:ime,jms:jme,kms:kme) ! Stream function
   real, intent(inout) :: chi(ims:ime,jms:jme,kms:kme) ! Velocity potential
   real, intent(in)    :: coefx(ims:ime,jms:jme)       ! Multiplicative coeff.
   real, intent(in)    :: coefy(ims:ime,jms:jme)       ! Multiplicative coeff.
   real, intent(out)   :: u(ims:ime,jms:jme,kms:kme)   ! u wind comp.
   real, intent(out)   :: v(ims:ime,jms:jme,kms:kme)   ! v wind comp.


   integer            :: i, j, k                      ! Loop counters.
   integer            :: is, ie                       ! 1st dim. end points.
   integer            :: js, je                       ! 2nd dim. end points.

   if (trace_use) call da_trace_entry("da_psichi_to_uv")

   !---------------------------------------------------------------------------
   ! [1.0] For Global application, set Wast-Eest Periodic boundary
   !---------------------------------------------------------------------------

   ! Computation to check for edge of domain:
   is = its
   ie = ite
   js = jts
   je = jte
   if (jts == jds) js = jds+1
   if (jte == jde) je = jde-1

   if (global) then 
      call da_set_boundary_3d(psi)
      call da_set_boundary_3d(chi)
   else
      if (its == ids) is = ids+1
      if (ite == ide) ie = ide-1
   end if


   !$OMP PARALLEL DO &
   !$OMP PRIVATE (i, j, k)
   do k = kts, kte
      !------------------------------------------------------------------------
      !  [2.0] Compute u, v at interior points (2nd order central finite diffs):
      !------------------------------------------------------------------------
      do j = js, je
         do i = is, ie
            u(i,j,k) = -coefy(i,j)*(psi(i  ,j+1,k) - psi(i  ,j-1,k)) + &
                        coefx(i,j)*(chi(i+1,j  ,k) - chi(i-1,j  ,k)) 

            v(i,j,k) = coefx(i,j)*(psi(i+1,j  ,k) - psi(i-1,j  ,k))  + &
                       coefy(i,j)*(chi(i  ,j+1,k) - chi(i  ,j-1,k)) 
         end do
      end do

      if (global) cycle
      !------------------------------------------------------------------------
      ! [3.0] Compute u, v at domain boundaries:
      !------------------------------------------------------------------------

      ! [3.1] Western boundaries:

      if (its == ids) then
         i = its
         do j = js, je
            u(i,j,k) = -coefy(i,j)*(psi(i  ,j+1,k) - psi(i,j-1,k)) + &
                        coefx(i,j)*(chi(i+2,j  ,k) - chi(i,j  ,k))  

            v(i,j,k) = coefx(i,j)*(psi(i+2,j  ,k) - psi(i,j  ,k))  + &
                       coefy(i,j)*(chi(i  ,j+1,k) - chi(i,j-1,k)) 
         end do
      end if

      ! [3.2] Eastern boundaries:

      if (ite == ide) then
         i = ite
         do j = js, je
            u(i,j,k) = -coefy(i,j)*(psi(i,j+1,k) - psi(i  ,j-1,k)) + &
                        coefx(i,j)*(chi(i,j  ,k) - chi(i-2,j  ,k)) 

            v(i,j,k) = coefx(i,j)*(psi(i,j  ,k) - psi(i-2,j  ,k))+ &
                       coefy(i,j)*(chi(i,j+1,k) - chi(i  ,j-1,k)) 
         end do
      end if

      ! [3.3] Southern boundaries:

      if (jts == jds) then
         j = jts
         do i = is, ie
            u(i,j,k) = -coefy(i,j)*(psi(i  ,j+2,k) - psi(i  ,j,k)) + &
                        coefx(i,j)*(chi(i+1,j  ,k) - chi(i-1,j,k))  

            v(i,j,k) = coefx(i,j)*(psi(i+1,j  ,k) - psi(i-1,j,k)) + &
                       coefy(i,j)*(chi(i  ,j+2,k) - chi(i  ,j,k)) 
         end do
      end if

      ! [3.4] Northern boundaries:

      if (jte == jde) then
         j = jte
         do i = is, ie
            u(i,j,k) = -coefy(i,j)*(psi(i  ,j,k) - psi(i  ,j-2,k)) + &
                        coefx(i,j)*(chi(i+1,j,k) - chi(i-1,j  ,k))  

            v(i,j,k) = coefx(i,j)*(psi(i+1,j,k) - psi(i-1,j  ,k))+ &
                       coefy(i,j)*(chi(i  ,j,k) - chi(i  ,j-2,k)) 
         end do
      end if

      !------------------------------------------------------------------------
      ! [4.0] Corner points (assume average of surrounding points - poor?):
      !------------------------------------------------------------------------

      ! [4.1] Bottom-left point:

      if (its == ids .AND. jts == jds) then
         u(its,jts,k) = 0.5 * (u(its+1,jts,k) + u(its,jts+1,k))
         v(its,jts,k) = 0.5 * (v(its+1,jts,k) + v(its,jts+1,k))
      end if

      ! [4.2] Top-left point:

      if (ite == ide .AND. jts == jds) then
         u(ite,jts,k) = 0.5 * (u(ite-1,jts,k) + u(ite,jts+1,k))
         v(ite,jts,k) = 0.5 * (v(ite-1,jts,k) + v(ite,jts+1,k))
      endif

      ! [4.3] Bottom-right point:

      if (its == ids .AND. jte == jde) then
         u(its,jte,k) = 0.5 * (u(its+1,jte,k) + u(its,jte-1,k))
         v(its,jte,k) = 0.5 * (v(its+1,jte,k) + v(its,jte-1,k))
      end if

      ! [4.4] Top-right point:

      if (ite == ide .AND. jte == jde) then
         u(ite,jte,k) = 0.5 * (u(ite-1,jte,k) + u(ite,jte-1,k))
         v(ite,jte,k) = 0.5 * (v(ite-1,jte,k) + v(ite,jte-1,k))
      end if
   end do
   !$OMP END PARALLEL DO

   !---------------------------------------------------------------------------
   ! [5.0] For Global application, set Wast-Eest Periodic boundary
   !---------------------------------------------------------------------------

   if (global) then
      call da_set_boundary_3d(u)
      call da_set_boundary_3d(v)
   end if

   if (trace_use) call da_trace_exit("da_psichi_to_uv")

end subroutine da_psichi_to_uv


subroutine da_psichi_to_uv_adj(u, v, coefx, coefy, psi, chi)

   !---------------------------------------------------------------------------
   !  Purpose: Adjoint code of da_psichi_to_uv  
   !
   !  Method:  u = coefx * (-dpsi/dy + dchi/dx)
   !           v = coefy * ( dpsi/dx + dchi/dy)
   !
   !  Assumptions: Unstaggered grid.
   !               Lateral boundary conditions - dpsi/dn, dchi/dn = 0 (FCT)
   !               grid size may or may not be equal
   ! 
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(inout) :: u(ims:ime,jms:jme,kms:kme)   ! u wind comp.
   real, intent(inout) :: v(ims:ime,jms:jme,kms:kme)   ! v wind comp.
   real, intent(inout) :: psi(ims:ime,jms:jme,kms:kme) ! Stream function
   real, intent(inout) :: chi(ims:ime,jms:jme,kms:kme) ! Velocity potential
   real, intent(in)    :: coefx(ims:ime,jms:jme)       ! Multiplicative coeff.
   real, intent(in)    :: coefy(ims:ime,jms:jme)       ! Multiplicative coeff.
   
   integer :: i, j, k                      ! Loop counters.
   integer :: is, ie                       ! 1st dim. end points.
   integer :: js, je                       ! 2nd dim. end points.

   if (trace_use) call da_trace_entry("da_psichi_to_uv_adj")

   !---------------------------------------------------------------------------
   ! [5.0] For Global application, set Wast-Eest Periodic boundary
   !---------------------------------------------------------------------------

   if (global) then
      call da_set_boundary_3d(u) 
      call da_set_boundary_3d(v) 
   end if

   !---------------------------------------------------------------------------
   ! Computation to check for edge of domain:
   !---------------------------------------------------------------------------

   is = its-1
   ie = ite+1
   js = jts-1
   je = jte+1
   if (jts == jds) js = jds+1
   if (jte == jde) je = jde-1

   if (.not. global) then
      if (its == ids) is = ids+1
      if (ite == ide) ie = ide-1

      !$OMP PARALLEL DO &
      !$OMP PRIVATE ( k ) 
      do k = kts, kte

         !---------------------------------------------------------------------
         ! [4.0] Corner points (assume average of surrounding points - poor?):
         !---------------------------------------------------------------------

         ! [4.1] Bottom-left point:

         if (its == ids .AND. jts == jds) then
            u(its+1,jts,k) = u(its+1,jts,k) + 0.5 * u(its,jts,k)
            u(its,jts+1,k) = u(its,jts+1,k) + 0.5 * u(its,jts,k)
            v(its+1,jts,k) = v(its+1,jts,k) + 0.5 * v(its,jts,k)
            v(its,jts+1,k) = v(its,jts+1,k) + 0.5 * v(its,jts,k)
         end if

        ! [4.2] Top-left point:

         if (ite == ide .AND. jts == jds) then
            u(ite-1,jts,k) = u(ite-1,jts,k) + 0.5 * u(ite,jts,k)
            u(ite,jts+1,k) = u(ite,jts+1,k) + 0.5 * u(ite,jts,k)
            v(ite-1,jts,k) = v(ite-1,jts,k) + 0.5 * v(ite,jts,k)
            v(ite,jts+1,k) = v(ite,jts+1,k) + 0.5 * v(ite,jts,k)
         end if

         ! [4.3] Bottom-right point:

         if (its == ids .AND. jte == jde) then
            u(its+1,jte,k) = u(its+1,jte,k) + 0.5 * u(its,jte,k)
            u(its,jte-1,k) = u(its,jte-1,k) + 0.5 * u(its,jte,k)
            v(its+1,jte,k) = v(its+1,jte,k) + 0.5 * v(its,jte,k)
            v(its,jte-1,k) = v(its,jte-1,k) + 0.5 * v(its,jte,k)
         end if
         ! [4.4] Top-right point:

         if (ite == ide .AND. jte == jde) then
            u(ite-1,jte,k) = u(ite-1,jte,k) + 0.5 * u(ite,jte,k)
            u(ite,jte-1,k) = u(ite,jte-1,k) + 0.5 * u(ite,jte,k)
            v(ite-1,jte,k) = v(ite-1,jte,k) + 0.5 * v(ite,jte,k)
            v(ite,jte-1,k) = v(ite,jte-1,k) + 0.5 * v(ite,jte,k)
         end if
      end do
      !$OMP END PARALLEL DO
   end if


   ! [3.0] Compute u, v at domain boundaries:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( i, j, k, message ) 
   do k = kts, kte
   if (.not. global) then
         ! [3.4] Northern boundaries:

         if (jte == jde) then
            j = jte

            do i = ie, is, -1
               psi(i+1,j,k) = psi(i+1,j,k) + coefx(i,j) * v(i,j,k)
               psi(i-1,j,k) = psi(i-1,j,k) - coefx(i,j) * v(i,j,k)
               chi(i,j  ,k) = chi(i,j  ,k) + coefy(i,j) * v(i,j,k)
               chi(i,j-2,k) = chi(i,j-2,k) - coefy(i,j) * v(i,j,k)

               psi(i,j  ,k) = psi(i,j  ,k) - coefy(i,j) * u(i,j,k)
               psi(i,j-2,k) = psi(i,j-2,k) + coefy(i,j) * u(i,j,k)
               chi(i+1,j,k) = chi(i+1,j,k) + coefx(i,j) * u(i,j,k)
               chi(i-1,j,k) = chi(i-1,j,k) - coefx(i,j) * u(i,j,k)
            end do
         end if
         ! [3.3] Southern boundaries:

         if (jts == jds) then
            j = jts

            do i = ie, is, -1


               psi(i+1,j,k) = psi(i+1,j,k) + coefx(i,j) * v(i,j,k)
               psi(i-1,j,k) = psi(i-1,j,k) - coefx(i,j) * v(i,j,k)
               chi(i,j+2,k) = chi(i,j+2,k) + coefy(i,j) * v(i,j,k)
               chi(i,j  ,k) = chi(i,j  ,k) - coefy(i,j) * v(i,j,k)

               psi(i,j+2,k) = psi(i,j+2,k) - coefy(i,j) * u(i,j,k)
               psi(i,j  ,k) = psi(i,j  ,k) + coefy(i,j) * u(i,j,k)
               chi(i+1,j,k) = chi(i+1,j,k) + coefx(i,j) * u(i,j,k)
               chi(i-1,j,k) = chi(i-1,j,k) - coefx(i,j) * u(i,j,k)
            end do
         end if

         ! [3.2] Eastern boundaries:

         if (ite == ide) then
            i = ite
            do j = je, js, -1
               psi(i  ,j,k) = psi(i  ,j,k) + coefx(i,j) * v(i,j,k)
               psi(i-2,j,k) = psi(i-2,j,k) - coefx(i,j) * v(i,j,k)
               chi(i,j+1,k) = chi(i,j+1,k) + coefy(i,j) * v(i,j,k)
               chi(i,j-1,k) = chi(i,j-1,k) - coefy(i,j) * v(i,j,k)

               psi(i,j+1,k) = psi(i,j+1,k) - coefy(i,j) * u(i,j,k)
               psi(i,j-1,k) = psi(i,j-1,k) + coefy(i,j) * u(i,j,k)
               chi(i  ,j,k) = chi(i  ,j,k) + coefx(i,j) * u(i,j,k)
               chi(i-2,j,k) = chi(i-2,j,k) - coefx(i,j) * u(i,j,k)
            end do
         end if

         ! [3.1] Western boundaries:
         if (its == ids) then
            i = its

            do j = je, js, -1
               psi(i+2,j,k) = psi(i+2,j,k) + coefx(i,j) * v(i,j,k)
               psi(i  ,j,k) = psi(i  ,j,k) - coefx(i,j) * v(i,j,k)
               chi(i,j+1,k) = chi(i,j+1,k) + coefy(i,j) * v(i,j,k)
               chi(i,j-1,k) = chi(i,j-1,k) - coefy(i,j) * v(i,j,k)

               psi(i,j+1,k) = psi(i,j+1,k) - coefy(i,j) * u(i,j,k)
               psi(i,j-1,k) = psi(i,j-1,k) + coefy(i,j) * u(i,j,k)
               chi(i+2,j,k) = chi(i+2,j,k) + coefx(i,j) * u(i,j,k)
               chi(i,  j,k) = chi(i,  j,k) - coefx(i,j) * u(i,j,k)
            end do
         end if
   end if   

      !------------------------------------------------------------------------
      ! [2.0] Compute u, v at interior points (2nd order central finite diffs):
      !------------------------------------------------------------------------
      do j = je, js, -1
         do i = ie, is, -1
            psi(i+1,j,k) = psi(i+1,j,k) + coefx(i,j) * v(i,j,k)
            psi(i-1,j,k) = psi(i-1,j,k) - coefx(i,j) * v(i,j,k)
            chi(i,j+1,k) = chi(i,j+1,k) + coefy(i,j) * v(i,j,k)
            chi(i,j-1,k) = chi(i,j-1,k) - coefy(i,j) * v(i,j,k)

            psi(i,j+1,k) = psi(i,j+1,k) - coefy(i,j) * u(i,j,k)
            psi(i,j-1,k) = psi(i,j-1,k) + coefy(i,j) * u(i,j,k)
            chi(i+1,j,k) = chi(i+1,j,k) + coefx(i,j) * u(i,j,k)
            chi(i-1,j,k) = chi(i-1,j,k) - coefx(i,j) * u(i,j,k)
         end do
      end do
   end do    !  loop over levels
   !$OMP END PARALLEL DO

   !---------------------------------------------------------------------------
   ! [5.0] For Global application, set Wast-Eest Periodic boundary
   !---------------------------------------------------------------------------

   if (global) then
       call da_set_boundary_3d(psi)
       call da_set_boundary_3d(chi)
   end if

   if (trace_use) call da_trace_exit("da_psichi_to_uv_adj")

end subroutine da_psichi_to_uv_adj
subroutine da_uv_to_divergence(xb, u, v, div)

   !---------------------------------------------------------------------------
   !  Purpose: Calculate divergence on a co-ordinate surface, given an input
   !           wind field.
   !
   !                        d   U      d   V
   !           Div = m^2 *[---(---) + ---(---) ] 
   !                        dx  m      dy  M
   !---------------------------------------------------------------------------

   implicit none

   type (xb_type), intent(in)           :: xb         ! First guess structure.
   real, intent(in)   :: u(ims:ime,jms:jme,kms:kme)   ! u wind comp.
   real, intent(in)   :: v(ims:ime,jms:jme,kms:kme)   ! v wind comp.
   real, intent(inout):: div(ims:ime,jms:jme,kms:kme) ! Divergence.

   integer            :: i, j, k                      ! Loop counters.
   integer            :: is, ie                       ! 1st dim. end points.
   integer            :: js, je                       ! 2nd dim. end points.
   real               :: one_third                    ! 1/3.

   real               :: coeff, inv_2ds
   real               :: um(ims:ime,jms:jme)          ! Temp. storage of u/m.
   real               :: vm(ims:ime,jms:jme)          ! Temp. storage of v/m.

   if (trace_use) call da_trace_entry("da_uv_to_divergence")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   one_third = 1.0 / 3.0
   div = 0.0

   !---------------------------------------------------------------------------
   ! Computation to check for edge of domain:
   !---------------------------------------------------------------------------

   is = its; ie = ite; js = jts; je = jte
   if (.not. global .and. its == ids) is = ids+1  
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds) js = jds+1; if (jte == jde) je = jde-1
  
   if (.not.global) inv_2ds = 0.5 / xb%ds

   !---------------------------------------------------------------------------
   ! [2.0] Calculate divergence:
   !---------------------------------------------------------------------------

   if (global) then
      do k = kts, kte
         ! [2.1] Compute fd divergence at interior points:

         do j = js, je
            do i = is, ie
               div(i,j,k) = xb%coefx(i,j) * (u(i+1,j,k) - u(i-1,j,k)) + &
                            xb%coefy(i,j) * (v(i,j+1,k) - v(i,j-1,k))
            end do
         end do
      end do
      call da_set_boundary_3d(div)
   else
      do k = kts, kte

         um(is-1:ie+1,js-1:je+1) = u(is-1:ie+1,js-1:je+1,k) / xb%map_factor(is-1:ie+1,js-1:je+1)
         vm(is-1:ie+1,js-1:je+1) = v(is-1:ie+1,js-1:je+1,k) / xb%map_factor(is-1:ie+1,js-1:je+1)

         ! [2.1] Compute fd divergence at interior points:

         do j = js, je
            do i = is, ie
               coeff = xb%map_factor(i,j) * xb%map_factor(i,j) * inv_2ds
               div(i,j,k) = (um(i+1,j) - um(i-1,j) + vm(i,j+1) - vm(i,j-1)) * coeff
            end do
         end do

         ! [2.2] Impose zero divergence gradient condition at boundaries:

         ! [2.2.1] Bottom boundaries:

         if (its == ids) then
            i = its 
            do j = jts, jte
               div(i,j,k) = one_third * (4.0 * div(i+1,j,k) - div(i+2,j,k))
            end do
         end if

         ! [2.2.2] Top boundaries:

         if (ite == ide) then
            i = ite
            do j = jts, jte
               div(i,j,k) = one_third * (4.0 * div(i-1,j,k) - div(i-2,j,k))
            end do
         end if

         ! [2.2.3] Left boundaries:

         if (jts == jds) then
            j = jts
            do i = its, ite
               div(i,j,k) = one_third * (4.0 * div(i,j+1,k) - div(i,j+2,k))
            end do
         end if

         ! [2.2.4] right boundaries:

         if (jte == jde) then
            j = jte
            do i = its, ite
               div(i,j,k) = one_third * (4.0 * div(i,j-1,k) - div(i,j-2,k))
            end do
         end if
      end do
   end if

   if (trace_use) call da_trace_exit("da_uv_to_divergence")

end subroutine da_uv_to_divergence


subroutine da_uv_to_divergence_adj(grid, u, v, div)

   !---------------------------------------------------------------------------
   !  Purpose: Adjoint of the subroutine da_uv_to_divergence
   !                        d   U      d   V
   !           Div = m^2 *[---(---) + ---(---) ] 
   !                        dx  m      dy  M
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   real, intent(out)  :: u(ims:ime,jms:jme,kms:kme)   ! u wind comp.
   real, intent(out)  :: v(ims:ime,jms:jme,kms:kme)   ! v wind comp.
   real, intent(inout):: div(ims:ime,jms:jme,kms:kme) ! Divergence.

   integer            :: i, j, k                      ! Loop counters.
   integer            :: is, ie                       ! 1st dim. end points.
   integer            :: js, je                       ! 2nd dim. end points.
   real               :: one_third                    ! 1/3.

   real               :: coeff, inv_2ds
   real               :: um(ims:ime,jms:jme)          ! Temp. storage of u/m.
   real               :: vm(ims:ime,jms:jme)          ! Temp. storage of v/m

   if (trace_use) call da_trace_entry("da_uv_to_divergence_adj")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   one_third = 1.0 / 3.0

   is = its - 1; ie = ite + 1; js = jts - 1; je = jte + 1

   if (.not. global .and. its == ids) is = ids+1
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds) js = jds+1; if (jte == jde) je = jde-1

   !---------------------------------------------------------------------------
   ! [2.0] Calculate divergence:
   !---------------------------------------------------------------------------

   if (.not. global) then
      inv_2ds = 0.5 / grid%xb%ds
      do k =kds, kde
         um(ims:ime,jms:jme) = 0.0
         vm(ims:ime,jms:jme) = 0.0
         ! [2.2] Impose zero divergence gradient condition at boundaries:

         ! [2.2.4] Right boundaries:

         if (jte == jde) then
            j = jte
            do i = its, ite    ! This is different to original
               div(i,j-1,k)=div(i,j-1,k)+div(i,j,k)*one_third*4.0
               div(i,j-2,k)=div(i,j-2,k)-div(i,j,k)*one_third
               div(i,j,k)=0.0
            end do
         end if

         ! [2.2.3] Left boundaries:

         if (jts == jds) then
            j = jts
            do i = its, ite    ! This is different to original
               div(i,j+1,k)=div(i,j+1,k)+div(i,j,k)*one_third*4.0
               div(i,j+2,k)=div(i,j+2,k)-div(i,j,k)*one_third
               div(i,j,k)=0.0
            end do
         end if

         ! [2.2.2] Top boundaries:

         if (ite == ide) then
            i = ite
            do j = jts, jte
               div(i-1,j,k)=div(i-1,j,k)+div(i,j,k)*one_third*4.0
               div(i-2,j,k)=div(i-2,j,k)-div(i,j,k)*one_third
               div(i,j,k)=0.0
            end do
         end if

         ! [2.2.1] Bottom boundaries:

         if (its == ids) then
            i = its 
            do j = jts, jte
               div(i+1,j,k)=div(i+1,j,k)+div(i,j,k)*one_third*4.0
               div(i+2,j,k)=div(i+2,j,k)-div(i,j,k)*one_third
               div(i,j,k)=0.0
            end do
         end if

         ! [2.1] Compute fd divergence at interior points:
         ! Computation to check for edge of domain:
         ! This is only for adjoint, as we have to cross the processor boundary
         ! to get the contribution.

         grid%xp%vxy(its:ite, jts:jte) = div(its:ite, jts:jte, k)
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_2D_WORK.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_2D_WORK_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE

         div(is, js:je, k) = grid%xp%vxy(is, js:je)
         div(ie, js:je, k) = grid%xp%vxy(ie, js:je)
         div(is:ie, js, k) = grid%xp%vxy(is:ie, js)
         div(is:ie, je, k) = grid%xp%vxy(is:ie, je)

         do j = js, je
            do i = is, ie
               coeff = grid%xb%map_factor(i,j) * grid%xb%map_factor(i,j) * inv_2ds
               um(i+1,j)=um(i+1,j)+div(i,j,k)*coeff
               um(i-1,j)=um(i-1,j)-div(i,j,k)*coeff
               vm(i,j+1)=vm(i,j+1)+div(i,j,k)*coeff
               vm(i,j-1)=vm(i,j-1)-div(i,j,k)*coeff
            end do
         end do

         u(is-1:ie+1,js-1:je+1,k) = u(is-1:ie+1,js-1:je+1,k) +um(is-1:ie+1,js-1:je+1) / grid%xb%map_factor(is-1:ie+1,js-1:je+1)
         v(is-1:ie+1,js-1:je+1,k) = v(is-1:ie+1,js-1:je+1,k) +vm(is-1:ie+1,js-1:je+1) / grid%xb%map_factor(is-1:ie+1,js-1:je+1)
      end do

   else ! global
      call da_set_boundary_3d(div)

      do k =kds, kde
         !-------------------------------------------------------------------------
         ! [2.1] Compute fd divergence at interior points:
         !-------------------------------------------------------------------------

         do j = je, js, -1
            do i = ie, is, -1  
               u(i+1,j,k) = u(i+1,j,k) + grid%xb%coefx(i,j) * div(i,j,k) 
               u(i-1,j,k) = u(i-1,j,k) - grid%xb%coefx(i,j) * div(i,j,k) 
               v(i,j+1,k) = v(i,j+1,k) + grid%xb%coefy(i,j) * div(i,j,k) 
               v(i,j-1,k) = v(i,j-1,k) - grid%xb%coefy(i,j) * div(i,j,k) 
            end do
         end do
      end do
   end if

   div = 0.0

   if (trace_use) call da_trace_exit("da_uv_to_divergence_adj")

end subroutine da_uv_to_divergence_adj


subroutine da_w_adjustment_lin(xb,W_a,WZ_a)

   !---------------------------------------------------------------------------
   ! Purpose: Adjust vertical velocity increments
   !
   ! Assumptions: 1) Model level stored top down.
   !---------------------------------------------------------------------------

   implicit none

   type (xb_type), intent(in)    :: xb                ! first guess structure.
   real,           intent(out)   :: w_a(ims:ime,jms:jme,kms:kme)
   real,           intent(inout) :: wz_a(ims:ime,jms:jme,kms:kme)

   integer :: i,j,k
   real    :: wz_b(ims:ime,jms:jme,kms:kme)
   real    :: ebxl1,ebxl2
   real    :: ebxl19,ebxl29

   if (trace_use) call da_trace_entry("da_w_adjustment_lin")

   call da_wz_base(xb,WZ_b)

   do j=jts,jte
      do i=its,ite
         ebxl1=0.0
         ebxl19=0.0
         ebxl2=0.0
         ebxl29=0.0
         do k=kte,kts,-1
            ebxl1=ebxl1+wz_a(i,j,k)*(xb%hf(i,j,k)-xb%hf(i,j,k+1))
            ebxl19=ebxl19+wz_b(i,j,k)*(xb%hf(i,j,k)-xb%hf(i,j,k+1))
            ebxl2=ebxl2+wz_a(i,j,k)*(xb%hf(i,j,k)-xb%hf(i,j,k+1))   &
                  *sign(1.0,wz_b(i,j,k))
            ebxl29=ebxl29+abs(wz_b(i,j,k))*(xb%hf(i,j,k)-xb%hf(i,j,k+1))
         end do

         do k=kts,kte
            wz_a(i,j,k)=wz_a(i,j,k)*(1.-ebxl19/ebxl29*sign(1.0,wz_b(i,j,k)))-  &
                        ebxl1*abs(wz_b(i,j,k))/ebxl29+                        &
                        ebxl2*abs(wz_b(i,j,k))*ebxl19/ebxl29**2
         end do

      end do
   end do

   if (trace_use) call da_trace_exit("da_w_adjustment_lin")

end subroutine da_w_adjustment_lin


subroutine da_w_adjustment_adj(xb,WZ_a)

   !---------------------------------------------------------------------------
   ! Purpose: Adjust vertical velocity increments
   !
   ! Assumptions: 1) Model level stored top down.
   !---------------------------------------------------------------------------

   implicit none

   type (xb_type), intent(in)    :: xb                ! First guess structure.

   real, dimension(ims:ime,jms:jme,kms:kme), intent(inout) :: WZ_a

   integer :: I,J,K

   real, dimension(ims:ime,jms:jme,kms:kme)   :: WZ_b

   real :: EBXL1,EBXL2
   real :: EBXL19,EBXL29

   if (trace_use) call da_trace_entry("da_w_adjustment_adj")

   call da_wz_base(xb,WZ_b)

   do J=jts,jte
      do I=its,ite
         EBXL19=0.0
         EBXL29=0.0

         do K=kte,kts,-1
            EBXL19=EBXL19+WZ_b(I,J,K)*(xb%hf(I,J,K)-xb%hf(I,J,K+1))
            EBXL29=EBXL29+ABS(WZ_b(I,J,K))*(xb%hf(I,J,K)-xb%hf(I,J,K+1))
         end do

         EBXL1=0.0
         EBXL2=0.0

         do K=kts,kte
            EBXL1=EBXL1-WZ_a(I,J,K)*ABS(WZ_b(I,J,K))/EBXL29
            EBXL2=EBXL2-WZ_a(I,J,K)*   &
                  ABS(WZ_b(I,J,K))*(-EBXL19)/EBXL29**2
            WZ_a(I,J,K)=WZ_a(I,J,K)*(1.0-EBXL19/EBXL29   &
                                   *SIGN(1.0,WZ_b(I,J,K)))
         end do

         do K=kte,kts,-1
            WZ_a(I,J,K)=WZ_a(I,J,K)+EBXL2*(xb%hf(I,J,K)-xb%hf(I,J,K+1))   &
                          *SIGN(1.0,WZ_b(I,J,K))
            WZ_a(I,J,K)=WZ_a(I,J,K)+EBXL1*(xb%hf(I,J,K)-xb%hf(I,J,K+1))
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_w_adjustment_adj")

end subroutine da_w_adjustment_adj


subroutine da_uv_to_vorticity(xb , u, v, vor)

   !---------------------------------------------------------------------------
   !  Purpose: Calculate vorticity on a co-ordinate surface, given an input 
   !           wind field.
   !
   !  Method:  vor = m**2 (d(v/m)/dx - d(u/m)/dy)
   !---------------------------------------------------------------------------

   implicit none

   type (xb_type), intent(in)    :: xb         ! First guess structure.
   real,           intent(in)    :: u(ims:ime,jms:jme,kms:kme)   ! u wind comp.
   real,           intent(in)    :: v(ims:ime,jms:jme,kms:kme)   ! v wind comp.
   real,           intent(inout) :: vor(ims:ime,jms:jme,kms:kme) ! Vorticity.

   integer :: i, j, k                      ! Loop counters.
   integer :: is, ie                       ! 1st dim. end points.
   integer :: js, je                       ! 2nd dim. end points.
   real    :: one_third                    ! 1/3.
   real    :: inv_2ds                      ! 1/2ds.
   real    :: coeff(ims:ime,jms:jme)       ! Mult. coefficient.
   real    :: um(ims:ime,jms:jme)          ! Temp. storage of u/m.
   real    :: vm(ims:ime,jms:jme)          ! Temp. storage of v/m.

   if (trace_use) call da_trace_entry("da_uv_to_vorticity")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   one_third = 1.0 / 3.0

   vor = 0.0

   !---------------------------------------------------------------------------
   ! Computation to check for edge of domain:
   !---------------------------------------------------------------------------

   is = its
   ie = ite
   js = jts
   je = jte
   if (.not. global .and. its == ids) is = ids+1 
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds) js = jds+1
   if (jte == jde) je = jde-1

   if (.not.global) then
      inv_2ds = 0.5 / xb%ds
      coeff(its:ite,jts:jte) = xb%map_factor(its:ite,jts:jte) * &
         xb%map_factor(its:ite,jts:jte) * inv_2ds
   end if
     
   !---------------------------------------------------------------------------
   ! [2.0] Calculate vorticity:
   !---------------------------------------------------------------------------

   do k = kts, kte
      ! [2.1] Compute finite difference vorticity at interior points:

      if (global) then
         do j = js, je
            do i = is, ie
               vor(i,j,k) = xb%coefy(i,j) * (-u(i,j+1,k) + u(i,j-1,k))  + & 
                            xb%coefx(i,j) * ( v(i+1,j,k) - v(i-1,j,k))
            end do
         end do

         ! if (global) cycle
      else
         um(is-1:ie+1,js-1:je+1) = u(is-1:ie+1,js-1:je+1,k) / xb%map_factor(is-1:ie+1,js-1:je+1)
         vm(is-1:ie+1,js-1:je+1) = v(is-1:ie+1,js-1:je+1,k) / xb%map_factor(is-1:ie+1,js-1:je+1)

         ! [2.1] Compute finite difference vorticity at interior points:

         do j = js, je
            do i = is, ie
               vor(i,j,k) = (-um(i,j+1) + um(i,j-1) + vm(i+1,j) - vm(i-1,j)) * coeff(i,j)
            end do
         end do

         ! [2.2] Impose zero vorticity gradient condition at boundaries:

         ! [2.2.1] Bottom boundaries:

         if (its == ids) then
            i = its 
            do j = jts, jte
               vor(i,j,k) = one_third * (4.0 * vor(i+1,j,k) - vor(i+2,j,k))
            end do
         end if

         ! [2.2.2] Top boundaries:

         if (ite == ide) then
            i = ite
            do j = jts, jte
               vor(i,j,k) = one_third * (4.0 * vor(i-1,j,k) - vor(i-2,j,k))
            end do
         end if

         ! [2.2.3] Left boundaries:

         if (jts == jds) then
            j = jts
            do i = its, ite
               vor(i,j,k) = one_third * (4.0 * vor(i,j+1,k) - vor(i,j+2,k))
            end do
         end if

         ! [2.2.4] right boundaries:

         if (jte == jde) then
            j = jte
            do i = its, ite
               vor(i,j,k) = one_third * (4.0 * vor(i,j-1,k) - vor(i,j-2,k))
            end do
         end if
      end if
   end do

   if (global) then
      call da_set_boundary_3d(vor)
   end if

   if (trace_use) call da_trace_exit("da_uv_to_vorticity")

end subroutine da_uv_to_vorticity


subroutine da_wz_base(xb,WZ_b)

   !--------------------------------------------------------------------
   ! Purpose: TBD
   !--------------------------------------------------------------------

   implicit none

   type (xb_type), intent(in)    :: xb                ! First guess structure.

   real, dimension(ims:ime,jms:jme,kms:kme), intent(inout) :: WZ_b

   integer                       :: is, ie       ! 1st dim. end points.
   integer                       :: js, je       ! 2nd dim. end points.

   integer                       :: I,J,K

   real  ::  TERM3

   real, dimension(ims:ime,jms:jme,kms:kme) :: URHO, VRHO
   real, dimension(ims:ime,jms:jme,kms:kme) :: DIV

   if (trace_use) call da_trace_entry("da_wz_base")


   ! Computation to check for edge of domain:

   is = its
   ie = ite
   js = jts
   je = jte
   if (its == ids) is = ids+1
   if (ite == ide) ie = ide-1
   if (jts == jds) js = jds+1
   if (jte == jde) je = jde-1


   do K=kts,kte
      do J=js,je
         do I=is,ie
            WZ_b(I,J,K)=-xb%u(I,J,K)*(xb%p(I+1,J,K)-xb%p(I-1,J,K))*xb%coefx(I,J)
         end do
      end do

      do J=js,je
         do I=is,ie
            WZ_b(I,J,K)=WZ_b(I,J,K)-xb%v(I,J,K)*(xb%p(I,J+1,K)-xb%p(I,J-1,K))*xb%coefy(I,J)
         end do
      end do
   end do

   if (its == ids) then
      i = its
      do K=kts,kte
         do J=js,je
            WZ_b(I,J,K)=WZ_b(I+1,J,K)
         end do
      end do
   end if

   if (ite == ide) then
      i = ite
      do K=kts,kte
         do J=js,je
            WZ_b(I,J,K)=WZ_b(I-1,J,K)
         end do
      end do
   end if

   if (jts == jds) then
      j = jts
      do K=kts,kte
         do I=its, ite
            WZ_b(I,J,K)=WZ_b(I,J+1,K)
         end do
      end do
   end if

   if (jte == jde) then
      j = jte
      do K=kts,kte
         do I=its, ite
            WZ_b(I,J,K)=WZ_b(I,J-1,K)
         end do
      end do
   end if


   call da_uv_to_divergence(xb, xb%u,xb%v, DIV)

   do K=kts,kte
      do J=jts,jte
         do I=its,ite
            WZ_b(I,J,K)=WZ_b(I,J,K)-GAMMA*xb%p(I,J,K)*DIV(I,J,K)
         end do
      end do
   end do


   ! Computation to check for edge of domain:
   is = its-1; ie = ite+1; js = jts-1; je = jte+1
   if (its == ids) is = ids; if (ite == ide) ie = ide
   if (jts == jds) js = jds; if (jte == jde) je = jde

   do K=kts,kte
      do J=js,je
         do I=is,ie
            URHO(I,J,K)=xb%rho(I,J,K)*xb%u(I,J,K)
            VRHO(I,J,K)=xb%rho(I,J,K)*xb%v(I,J,K)
         end do
      end do
   end do

   call da_uv_to_divergence(xb, URHO, VRHO, DIV)

   do J=jts,jte
      do I=its,ite
         TERM3=0.0

         do K=kte-1,kts,-1
            TERM3=TERM3+GRAVITY*(DIV(I,J,K+1)+DIV(I,J,K))*0.5  &
                       *(xb%h(I,J,K+1)-xb%h(I,J,K))
            WZ_b(I,J,K)=WZ_b(I,J,K)+TERM3
         end do
      end do
   end do


   do K=kts,kte
      do J=jts,jte
         do I=its,ite
            WZ_b(I,J,K)=WZ_b(I,J,K)/(GAMMA*xb%p(I,J,K))
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_wz_base")

end subroutine da_wz_base


subroutine da_wpec_constraint(grid, xbx)

   !---------------------------------------------------------------------------
   !  Purpose: Calculates balance equation G(x)
   !---------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid
   type(xbx_type), intent(in)    :: xbx                            ! Header & non-gridded vars.

   real   :: u(ims:ime,jms:jme,kms:kme)     ! u wind comp.
   real   :: v(ims:ime,jms:jme,kms:kme)     ! v wind comp.
   real   :: phi_b(ims:ime,jms:jme,kms:kme) 

   integer :: i, j, k                 ! Loop counters.
   integer :: is, ie                  ! 1st dim. end points.
   integer :: js, je                  ! 2nd dim. end points.

   real    :: coefx(ims:ime,jms:jme)  ! Multiplicative coefficient.
   real    :: coefy(ims:ime,jms:jme)  ! Multiplicative coefficient.
   real    :: term_x(ims:ime,jms:jme) ! Balance eqn x term
   real    :: term_y(ims:ime,jms:jme) ! Balance eqn y term
   real    :: phi_b_x(ims:ime,jms:jme) ! Balance eqn x term
   real    :: phi_b_y(ims:ime,jms:jme) ! Balance eqn y term
   


   if (trace_use) call da_trace_entry("da_wpec_constraint")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise i and set multiplicative constants
   !---------------------------------------------------------------------------

   ! Computation to check for edge of domain:

   is = its; ie = ite; js = jts; je = jte
   if (.not.global .and. its == ids) is = ids+1
   if (.not.global .and. ite == ide) ie = ide-1
   if (jts == jds ) js = jds+1; if (jte == jde) je = jde-1


   if (fg_format == fg_format_kma_global) then
      coefx = grid%xb%coefx
      coefy = grid%xb%coefy
   else if( fg_format == fg_format_wrf_arw_regional) then
      coefx = grid%xb%coefx
      coefy = grid%xb%coefy
   else if (fg_format == fg_format_wrf_arw_global) then
      write (unit=message(1),fmt='(A,I3)') ' needs work for fg_format_wrf_arw_global  = ',fg_format
      call da_error("da_wpec_constraint.inc",51,message(1:1))
   else if (fg_format == fg_format_wrf_nmm_regional) then
      write (unit=message(1),fmt='(A,I3)') ' needs work for fg_format_wrf_nmm_regional = ',fg_format
      call da_error("da_wpec_constraint.inc",54,message(1:1))
   else
      write(unit=message(1),fmt='(A,I3)') 'Wrong FG_FORMAT = ',fg_format
      call da_error("da_wpec_constraint.inc",57,message(1:1))
   end if

   ! [1.1] 

   phi_b  =   grid%xb%p

   do k = kts, kte

      term_x(ims:ime,jms:jme) = 0.0
      term_y(ims:ime,jms:jme) = 0.0
      phi_b_x(ims:ime,jms:jme) = 0.0
      phi_b_y(ims:ime,jms:jme) = 0.0

      !---------------------------------------------------------------------------
      ! [2.0] Calculate RHS of balance equation in gridpoint space:
      !---------------------------------------------------------------------------

      ! [2.1] Include geostrophic terms in balance eqn if requested:

      if (balance_type == balance_geo .OR. balance_type == balance_geocyc ) then
         call da_wpec_constraint_geoterm (grid%xb % cori, grid%xb % rho(:,:,k), grid%xb %u(:,:,k), grid%xb %v(:,:,k), &
            term_x, term_y)
      end if
      
      ! [2.2] Include cyclostrophic terms in balance eqn if requested:

      if (balance_type == balance_cyc .OR. balance_type == balance_geocyc ) then
         call da_wpec_constraint_cycloterm (grid%xb % rho(:,:,k), grid%xb % u(:,:,k),   &
            grid%xb % v(:,:,k), grid%xb % coefx(:,:), grid%xb % coefy(:,:), &
            term_x, term_y)
      end if

      ! [2.3] Include phi_b terms in balance eqn
      do j = js, je
         do i = is, ie
            phi_b_x(i,j) = coefx(i,j)*(phi_b(i+1,j,k)-phi_b(i-1,j,k)) 
            phi_b_y(i,j) = coefy(i,j)*(phi_b(i,j+1,k)-phi_b(i,j-1,k))                 
         end do
      end do

   !------------------------------------------------------------------------------
   !  [3.0] Solve Grad_p for balance eqn :
   !------------------------------------------------------------------------------
      do j = js, je
         do i = is, ie
            grid%xb%xb_p_x(i,j,k)=phi_b_x(i,j)+term_x(i,j)
            grid%xb%xb_p_y(i,j,k)=phi_b_y(i,j)+term_y(i,j)
         end do
      end do

    end do


   if (trace_use) call da_trace_exit("da_wpec_constraint")

end subroutine da_wpec_constraint


subroutine da_wpec_constraint_adj(grid, xbx)

   !---------------------------------------------------------------------------
   ! Purpose: Calculates ADM of balance equation G(x)
   !---------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)               :: grid

   type (xbx_type),intent(in) :: xbx              ! Header & non-gridded vars.

   real :: p(ims:ime,jms:jme,kms:kme) ! pressure increment.
   real :: geoh(ims:ime,jms:jme,kms:kme) ! geopotential height increment.
   real :: u(ims:ime,jms:jme,kms:kme) ! u wind comp. (dot pts)
   real :: v(ims:ime,jms:jme,kms:kme) ! v wind comp. (dot pts)

   integer             :: i, j, k          ! Loop counters.
   integer             :: is, ie           ! 1st dim. end points.
   integer             :: js, je           ! 2nd dim. end points.

   real, dimension(ims:ime,jms:jme) :: coefx, &   ! Multiplicative coefficient.
                                       coefy, &   ! Multiplicative coefficient.
                                       term_x, &  ! Balance eqn x term
                                       term_y     ! Balance eqn y term
   real    :: phi_b_x(ims:ime,jms:jme) ! Balance eqn x term
   real    :: phi_b_y(ims:ime,jms:jme) ! Balance eqn y term

   real                 :: data(ims:ime,jms:jme)        ! Work array.
   real                 :: var, uar

   if (trace_use) call da_trace_entry("da_wpec_constraint_adj")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   is = its; ie = ite; js = jts; je = jte
   if (.not. global .and. its == ids) is = ids+1
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds ) js = jds+1
   if (jte == jde ) je = jde-1

   if (fg_format == fg_format_kma_global) then
      coefx = grid%xb%coefx
      coefy = grid%xb%coefy
   else if (fg_format == fg_format_wrf_arw_regional) then
      coefx = grid%xb%coefx
      coefy = grid%xb%coefy
   else if (fg_format == fg_format_wrf_arw_global) then
      write (unit=message(1),fmt='(A,I3)') ' needs work for fg_format  = ',fg_format
      call da_error("da_wpec_constraint_adj.inc",52,message(1:1))
   else if (fg_format == fg_format_wrf_nmm_regional) then
      write (unit=message(1),fmt='(A,I3)') ' needs work for fg_format  = ',fg_format
      call da_error("da_wpec_constraint_adj.inc",55,message(1:1))
   else
      write (unit=message(1),fmt='(A,I3)') ' Wrong FG_FORMAT = ',fg_format
      call da_error("da_wpec_constraint_adj.inc",58,message(1:1))
   end if
   
   u       = 0.0
   v       = 0.0
   p       = 0.0
   geoh    = 0.0

   do k = kts, kte

      term_x = 0.0
      term_y = 0.0
      phi_b_x = 0.0
      phi_b_y = 0.0

      !---------------------------------------------------------------------------
      ! [2.0] Solve Grad_p for balance eqn
      !---------------------------------------------------------------------------

      phi_b_x = grid%xa%grad_p_x(:,:,k)
      phi_b_y = grid%xa%grad_p_y(:,:,k)
      term_x = grid%xa%grad_p_x(:,:,k)
      term_y = grid%xa%grad_p_y(:,:,k)

      !---------------------------------------------------------------------------
      ! [3.0] Calculate RHS of balance equation in gridpt space
      !---------------------------------------------------------------------------

      ! [3.1] Include phi_b terms in balance eqn

      do j = je, js, -1
         do i = ie, is, -1

            p(i+1,j,k) = p(i+1,j,k) + coefx(i,j) * phi_b_x(i,j)
            p(i-1,j,k) = p(i-1,j,k) - coefx(i,j) * phi_b_x(i,j)
            p(i,j+1,k) = p(i,j+1,k) + coefy(i,j) * phi_b_y(i,j)
            p(i,j-1,k) = p(i,j-1,k) - coefy(i,j) * phi_b_y(i,j)

            geoh(i+1,j,k) = geoh(i+1,j,k) + coefx(i,j) * phi_b_x(i,j) * grid%xb % rho(i,j,k)
            geoh(i-1,j,k) = geoh(i-1,j,k) - coefx(i,j) * phi_b_x(i,j) * grid%xb % rho(i,j,k)
            geoh(i,j+1,k) = geoh(i,j+1,k) + coefy(i,j) * phi_b_y(i,j) * grid%xb % rho(i,j,k)
            geoh(i,j-1,k) = geoh(i,j-1,k) - coefy(i,j) * phi_b_y(i,j) * grid%xb % rho(i,j,k)

         end do
      end do

      ! [3.2] Include cyclostrophic terms in balance eqn if requested:

      if (balance_type == balance_cyc .OR. balance_type == balance_geocyc ) then

         ! [3.2.1] Calculate adjoint of term_y = rho M ( u'dv/dx + v'dv/dy + udv'/dx + vdv'/dy ):
         data = grid%xb%rho(:,:,k) * term_y

         do j = je, js, -1
            do i = ie, is, -1
               uar = coefx(i,j) * grid%xb%u(i,j,k) * data(i,j)
               var = coefy(i,j) * grid%xb%v(i,j,k) * data(i,j)

                u(i,j,k) = u(i,j,k) + coefx(i,j)*data(i,j)*( grid%xb%v(i+1,j,k) - grid%xb%v(i-1,j,k) )
                v(i,j,k) = v(i,j,k) + coefy(i,j)*data(i,j)*( grid%xb%v(i,j+1,k) - grid%xb%v(i,j-1,k) )
                v(i+1,j,k) = v(i+1,j,k) + uar
                v(i-1,j,k) = v(i-1,j,k) - uar
                v(i,j+1,k) = v(i,j+1,k) + var
                v(i,j-1,k) = v(i,j-1,k) - var
            end do
         end do

         ! [3.2.2] Calculate adjoint of term_x = rho M ( u'du/dx + v'du/dy + udu'/dx + vdu'/dy ):
         data = grid%xb%rho(:,:,k) * term_x

         do j = je, js, -1
            do i = ie, is, -1
               uar = coefx(i,j) * grid%xb%u(i,j,k) * data(i,j)
               var = coefy(i,j) * grid%xb%v(i,j,k) * data(i,j)

                u(i,j,k) = u(i,j,k) + coefx(i,j)*( grid%xb%u(i+1,j,k) - grid%xb%u(i-1,j,k) ) * data(i,j)
                v(i,j,k) = v(i,j,k) + coefy(i,j)*( grid%xb%u(i,j+1,k) - grid%xb%u(i,j-1,k) ) * data(i,j)
                u(i+1,j,k) = u(i+1,j,k) + uar
                u(i-1,j,k) = u(i-1,j,k) - uar
                u(i,j+1,k) = u(i,j+1,k) + var
                u(i,j-1,k) = u(i,j-1,k) - var
            end do
         end do
      end if

      
      ! [3.3] Calculate geostrophic terms in balance eqn:
 
      if (balance_type == balance_geo .OR. balance_type == balance_geocyc ) then
         ! [3.3.1] Calculate term_y = f rho u~:
         u(:,:,k) = u(:,:,k) + grid%xb%rho(:,:,k) * grid%xb%cori * term_y

         ! [3.3.2] Calculate term_x = -f rho v~:
         v(:,:,k) = v(:,:,k) - grid%xb%rho(:,:,k) * grid%xb%cori * term_x

      end if

   end do

   !---------------------------------------------------------------------------
   ! [4.0] Store results in grid%xa structure
   !---------------------------------------------------------------------------


   grid%xa%u=u
   grid%xa%v=v
   grid%xa%p=p
   grid%xa%geoh=geoh


   if (trace_use) call da_trace_exit("da_wpec_constraint_adj")

end subroutine da_wpec_constraint_adj


subroutine da_wpec_constraint_cycloterm (rho, u, v, coefx, coefy, term_x, term_y)

   !---------------------------------------------------------------------------
   !  Purpose: Calculates cyclostrophic term in wpec constraint equation.
   !
   !  Method:  Term is rho (u.grad) u on a single level.
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(in)    :: rho(ims:ime,jms:jme)    ! Density.
   real, intent(in)    :: u(ims:ime,jms:jme)      ! u wind increment
   real, intent(in)    :: v(ims:ime,jms:jme)      ! v wind increment
   real, intent(in)    :: coefx(ims:ime,jms:jme)  ! Multiplicative const.
   real, intent(in)    :: coefy(ims:ime,jms:jme)
   real, intent(inout) :: term_x(ims:ime,jms:jme) ! x component of term.
   real, intent(inout) :: term_y(ims:ime,jms:jme) ! y component of term.

   integer :: i, j                         ! Loop counters.
   integer :: is, ie                       ! 1st dim. end points.
   integer :: js, je                       ! 2nd dim. end points.
   
   real    :: data(ims:ime,jms:jme)        ! Temporary storage.

   real    :: varb

   if (trace_use) call da_trace_entry("da_wpec_constraint_cycloterm")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------
   
   ! Computation to check for edge of domain:
   is = its; ie = ite; js = jts; je = jte
   if (.not. global .and. its==ids) is = ids+1
   if (.not. global .and. ite==ide) ie = ide-1
   if (jts==jds) js = jds+1
   if (jte==jde) je = jde-1
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate term_x = rho M (u.du/dx + v.du/dy):
   !---------------------------------------------------------------------------

   ! [2.1] Interior points:

   do j = js, je
      do i = is, ie
         data(i,j) = u(i,j) * coefx(i,j)*(u(i+1,j) - &
            u(i-1,j)) + v(i,j) * coefy(i,j)*(u(i,j+1) - &
            u(i,j-1))
      end do
   end do
   

   !  [2.7] Multiply by rho  and add to term_x:
      
   term_x(its:ite,jts:jte) = rho(its:ite,jts:jte)*data(its:ite,jts:jte) + term_x(its:ite,jts:jte)

   !---------------------------------------------------------------------------
   ! [3.0] Calculate term_y = rho M (u.dv/dx + v.dv/dy):
   !---------------------------------------------------------------------------

   ! [3.1] Interior points:

   do j = js, je
      do i = is, ie
         data(i,j) = u(i,j) * coefx(i,j)*(v(i+1,j) - v(i-1,j)) + &
                     v(i,j) * coefy(i,j)*(v(i,j+1) - v(i,j-1))
      end do
   end do
   

   ! [3.7] Multiply by rho and add to term_y

   term_y(its:ite,jts:jte)=rho(its:ite,jts:jte)* data(its:ite,jts:jte) + &
      term_y(its:ite,jts:jte)

   if (trace_use) call da_trace_exit("da_wpec_constraint_cycloterm")

end subroutine da_wpec_constraint_cycloterm


subroutine da_wpec_constraint_geoterm (cori,rho, u, v, term_x, term_y)
 
   !---------------------------------------------------------------------------
   ! Purpose: calculates nonlinear geostrophic term in wpec constraint equation.
   !
   ! method:  term is k x rho f u on a single level.
   !---------------------------------------------------------------------------

   implicit none
   real,         intent(in)  :: cori(ims:ime,jms:jme)   ! Coriolis factor.
   real,         intent(in)  :: rho(ims:ime,jms:jme)    ! Density
   real,         intent(in)  :: u(:,:)       ! u wind comp. (dot pts)
   real,         intent(in)  :: v(:,:)       ! v wind comp. (dot pts)
   real,         intent(out) :: term_x(:,:)  ! x component of term.
   real,         intent(out) :: term_y(:,:)  ! y component of term.
   integer :: is, ie                       ! 1st dim. end points.
   integer :: js, je                       ! 2nd dim. end points.

   if (trace_use) call da_trace_entry("da_wpec_constraint_geoterm")

   !---------------------------------------------------------------------------
   !  [1.0] Initialise:
   !---------------------------------------------------------------------------

   ! Computation to check for edge of domain:
   is = its; ie = ite; js = jts; je = jte
   if (.not. global .and. its == ids) is = ids+1
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds) js = jds+1; if (jte == jde) je = jde-1

   !---------------------------------------------------------------------------
   !  [2.0] Calculate term_x = -f rho v~:
   !---------------------------------------------------------------------------

   term_x(is:ie, js:je) = term_x(is:ie, js:je)-rho(is:ie, js:je) &
      * v(is:ie, js:je) * cori(is:ie, js:je)

   !---------------------------------------------------------------------------
   !  [3.0] Calculate term_y = f rho u~:
   !---------------------------------------------------------------------------

   term_y(is:ie, js:je) = term_y(is:ie, js:je)+rho(is:ie, js:je) &
      * u(is:ie, js:je) * cori(is:ie, js:je)

   if (trace_use) call da_trace_exit("da_wpec_constraint_geoterm")

end subroutine da_wpec_constraint_geoterm


subroutine da_wpec_constraint_lin(grid, xbx)

   !---------------------------------------------------------------------------
   !  Purpose: Calculates TLM of balance equation G(x)
   !---------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid
   type(xbx_type), intent(in)    :: xbx                            ! Header & non-gridded vars.

   integer :: i, j, k                 ! Loop counters.
   integer :: is, ie                  ! 1st dim. end points.
   integer :: js, je                  ! 2nd dim. end points.

   real    :: coefx(ims:ime,jms:jme)  ! Multiplicative coefficient.
   real    :: coefy(ims:ime,jms:jme)  ! Multiplicative coefficient.
   real    :: term_x(ims:ime,jms:jme) ! Balance eqn x term
   real    :: term_y(ims:ime,jms:jme) ! Balance eqn y term
   real    :: phi_b_x(ims:ime,jms:jme) ! Balance eqn x term
   real    :: phi_b_y(ims:ime,jms:jme) ! Balance eqn y term
   

   if (trace_use) call da_trace_entry("da_wpec_constraint_lin")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise i and set multiplicative constants
   !---------------------------------------------------------------------------

   is = its; ie = ite; js = jts; je = jte
   if (.not.global .and. its == ids) is = ids+1
   if (.not.global .and. ite == ide) ie = ide-1
   if (jts == jds ) js = jds+1; if (jte == jde) je = jde-1

   if (fg_format == fg_format_kma_global) then
      coefx = grid%xb%coefx
      coefy = grid%xb%coefy
   else if( fg_format == fg_format_wrf_arw_regional) then
      coefx = grid%xb%coefx
      coefy = grid%xb%coefy
   else if (fg_format == fg_format_wrf_arw_global) then
      write (unit=message(1),fmt='(A,I3)') ' needs work for fg_format_wrf_arw_global  = ',fg_format
      call da_error("da_wpec_constraint_lin.inc",43,message(1:1))
   else if (fg_format == fg_format_wrf_nmm_regional) then
      write (unit=message(1),fmt='(A,I3)') ' needs work for fg_format_wrf_nmm_regional = ',fg_format
      call da_error("da_wpec_constraint_lin.inc",46,message(1:1))
   else
      write(unit=message(1),fmt='(A,I3)') 'Wrong FG_FORMAT = ',fg_format
      call da_error("da_wpec_constraint_lin.inc",49,message(1:1))
   end if

   ! [1.1] 

   phi_b_x = 0.0
   phi_b_y = 0.0

 
   do k = kts, kte

      term_x = 0.0
      term_y = 0.0

      !---------------------------------------------------------------------------
      ! [2.0] Calculate RHS of balance equation in gridpoint space:
      !---------------------------------------------------------------------------

      ! [2.1] Include geostrophic terms in balance eqn if requested:

      if (balance_type == balance_geo .OR. balance_type == balance_geocyc ) then

         ! [2.1.1] Calculate term_x = -f rho v~:
         term_x = term_x - grid%xb%rho(:,:,k) * grid%xb%cori * grid%xa%v(:,:,k)

         ! [2.1.2] Calculate term_y = f rho u~:
         term_y = term_y + grid%xb%rho(:,:,k) * grid%xb%cori * grid%xa%u(:,:,k)

      end if
      
      ! [2.2] Include cyclostrophic terms in balance eqn if requested:

      if (balance_type == balance_cyc .OR. balance_type == balance_geocyc ) then

         do j = js, je
            do i = is, ie
               ! [2.2.1] Calculate term_x = rho M ( u'du/dx + v'du/dy + udu'/dx + vdu'/dy )
               term_x(i,j) = term_x(i,j) + grid%xb%rho(i,j,k) * &
                   ( coefx(i,j)*grid%xa%u(i,j,k) * ( grid%xb%u(i+1,j,k) - grid%xb%u(i-1,j,k)) + &
                     coefy(i,j)*grid%xa%v(i,j,k) * ( grid%xb%u(i,j+1,k) - grid%xb%u(i,j-1,k)) + &
                     coefx(i,j)*grid%xb%u(i,j,k) * ( grid%xa%u(i+1,j,k) - grid%xa%u(i-1,j,k)) + &
                     coefy(i,j)*grid%xb%v(i,j,k) * ( grid%xa%u(i,j+1,k) - grid%xa%u(i,j-1,k)))

               ! [2.2.2] Calculate term_y = rho M ( u'dv/dx + v'dv/dy + udv'/dx + vdv'/dy )
               term_y(i,j) = term_y(i,j) + grid%xb%rho(i,j,k) * &
                   ( coefx(i,j)*grid%xa%u(i,j,k) * ( grid%xb%v(i+1,j,k) - grid%xb%v(i-1,j,k)) + &
                     coefy(i,j)*grid%xa%v(i,j,k) * ( grid%xb%v(i,j+1,k) - grid%xb%v(i,j-1,k)) + &
                     coefx(i,j)*grid%xb%u(i,j,k) * ( grid%xa%v(i+1,j,k) - grid%xa%v(i-1,j,k)) + &
                     coefy(i,j)*grid%xb%v(i,j,k) * ( grid%xa%v(i,j+1,k) - grid%xa%v(i,j-1,k)))
            end do
         end do

      end if

      ! [2.3] Include phi_b terms in balance eqn
      do j = js, je
         do i = is, ie
            phi_b_x(i,j) = coefx(i,j)*(grid%xa%p(i+1,j,k)-grid%xa%p(i-1,j,k) + grid%xb % rho(i,j,k)*(grid%xa%geoh(i+1,j,k)-grid%xa%geoh(i-1,j,k)) )
            phi_b_y(i,j) = coefy(i,j)*(grid%xa%p(i,j+1,k)-grid%xa%p(i,j-1,k) + grid%xb % rho(i,j,k)*(grid%xa%geoh(i,j+1,k)-grid%xa%geoh(i,j-1,k)) )                 
         end do
      end do

   !------------------------------------------------------------------------------
   !  [3.0] Solve Grad_p for balance equation :
   !------------------------------------------------------------------------------
      do j = js, je
         do i = is, ie
            grid%xa%grad_p_x(i,j,k)=phi_b_x(i,j)+term_x(i,j)
            grid%xa%grad_p_y(i,j,k)=phi_b_y(i,j)+term_y(i,j)
         end do
      end do

   end do

   if (trace_use) call da_trace_exit("da_wpec_constraint_lin")

end subroutine da_wpec_constraint_lin



end module da_dynamics

