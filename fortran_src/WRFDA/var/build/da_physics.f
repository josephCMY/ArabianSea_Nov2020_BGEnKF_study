












module da_physics

   
   
   

   use module_domain, only : domain, x_type
   use module_dm, only : local_communicator, mytask, ntasks, ntasks_x, &
      ntasks_y, data_order_xyz
   use module_comm_dm, only : halo_bal_eqn_adj_sub, halo_xa_cloud_sub
   use da_control, only : ips,ipe,jps,jpe,kps,kpe
   use da_define_structures, only : synop_type, residual_synop_type, infa_type, iv_type
   use da_control, only : gas_constant, gravity,kts,kte, svp3,svpt0, &
      a_ew,  wdk1, wdk2, zdk1, zdk2, zdk3, radian, &
      gas_constant_v, svp1, to, xls, svp2,its,ite,jts,jte,kts,kte, &
      ims,ime,jms,jme,kms,kme,xlv1,cp,ids,ide,jds,jde,kds,kde, test_transforms, &
      trace_use, missing_r, maximum_rh, minimum_rh,cv_options_hum,coeff,l_over_rv, &
      es_gammakelvin, es_gammabeta, rd_over_rv1,t_kelvin, es_alpha, es_gamma, &
      es_beta, rd_over_rv, trace_use_frequent,gamma, print_detail_xa, stdout, &
      cv_options_hum_specific_humidity, trace_use_dull, pi
   use da_par_util, only : da_transpose_z2y, da_transpose_y2x, &
      da_transpose_x2z, da_transpose_z2x, da_transpose_x2y, da_transpose_y2z
   use da_tracing, only : da_trace_entry, da_trace_exit
   use da_interpolation, only : da_interp_lin_2d, da_interp_lin_2d_adj
   use da_dynamics, only : da_w_adjustment_adj, da_uv_to_divergence_adj, &
      da_w_adjustment_lin, da_uv_to_divergence
   use da_reporting, only : da_error, message
   use da_wrf_interfaces, only : wrf_debug
   use da_grid_definitions, only : da_ffdduv_model

   implicit none

   REAL,   DIMENSION(0:1000 ),SAVE          :: PSIMTB,PSIHTB

   contains

subroutine sfclayinit

! adapted from phys/module_sf_sfclay.F
! to be used in subroutine da_sfc_wtq

   INTEGER                   ::      N
   REAL                      ::      ZOLN,X,Y

   DO N=0,1000
      ZOLN=-FLOAT(N)*0.01
      X=(1-16.*ZOLN)**0.25
      PSIMTB(N)=2*ALOG(0.5*(1+X))+ALOG(0.5*(1+X*X))- &
                2.*ATAN(X)+2.*ATAN(1.)
      Y=(1-16*ZOLN)**0.5
      PSIHTB(N)=2*ALOG(0.5*(1+Y))
   ENDDO
end subroutine sfclayinit
subroutine da_prho_to_t_adj(grid) 

   !---------------------------------------------------------------------------
   !  Purpose: Adjoint of da_prho_to_t.
   !
   !  Method:  Standard adjoint coding.
   !
   !  Assumptions: 1) Model level stored top down.
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain), intent(inout)  :: grid

   integer :: is, ie       ! 1st dim. end points.
   integer :: js, je       ! 2nd dim. end points.
   integer :: k            ! Loop counter.
   real    :: temp(its:ite,jts:jte) ! Temporary array.

   if (trace_use) call da_trace_entry("da_prho_to_t_adj")

   !---------------------------------------------------------------------------
   !  [1.0] initialise:
   !---------------------------------------------------------------------------

   is = its
   ie = ite
   js = jts
   je = jte
   
   if (test_transforms) then
      is = its-1
      js = jts-1

      ie = ite+1
      je = jte+1

      if (is < ids) is = ids
      if (js < jds) js = jds

      if (ie > ide) ie = ide
      if (je > jde) je = jde
   end if

   !---------------------------------------------------------------------------
   ! [2.0] Calculate temperature increments:
   !---------------------------------------------------------------------------

   do k = kts, kte
      temp(is:ie,js:je) = grid%xa % t(is:ie,js:je,k) * grid%xb % t(is:ie,js:je,k)

      grid%xa % p(is:ie,js:je,k) = grid%xa % p(is:ie,js:je,k) &
         + temp(is:ie,js:je) / grid%xb % p(is:ie,js:je,k)
                         
      grid%xa % rho(is:ie,js:je,k) = grid%xa % rho(is:ie,js:je,k) &
         - temp(is:ie,js:je) / grid%xb % rho(is:ie,js:je,k)
   end do  

   if (trace_use) call da_trace_exit("da_prho_to_t_adj")
                             
end subroutine da_prho_to_t_adj


subroutine da_prho_to_t_lin(grid) 

   !---------------------------------------------------------------------------
   ! Purpose: Calculates temperature increments from pressure/density increments
   !
   ! Method:  Linearitsed ideal gas law: T~ = T (p"/p - rho"/rho)
   !
   ! Assumptions: 1) Model level stored top down.
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain), intent(inout)    :: grid

   if (trace_use) call da_trace_entry("da_prho_to_t_lin")

   grid%xa % t(its:ite,jts:jte,kts:kte) = grid%xb % t(its:ite,jts:jte,kts:kte) * &
                            (grid%xa % p(its:ite,jts:jte,kts:kte) / &
                              grid%xb % p(its:ite,jts:jte,kts:kte) - &
                              grid%xa % rho(its:ite,jts:jte,kts:kte) / &
                              grid%xb % rho(its:ite,jts:jte,kts:kte)) 

   if (trace_use) call da_trace_exit("da_prho_to_t_lin")    

end subroutine da_prho_to_t_lin


subroutine da_uvprho_to_w_lin(grid)

   !------------------------------------------------------------------------------
   !  Purpose: Calculates vertical velocity increments from Richardson's Eq.
   !
   !  Method: Richardson's Eq., which
   !          combines continuity Eq., thermodynamic Eq. and hrdrostatic Eq.
   !------------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   integer :: is, ie       ! 1st dim. end points.
   integer :: js, je       ! 2nd dim. end points.

   integer :: I,J,K

   real    :: urho(ims:ime,jms:jme,kms:kme)
   real    :: vrho(ims:ime,jms:jme,kms:kme)
   real    :: div(ims:ime,jms:jme,kms:kme)
   real    :: wz(ims:ime,jms:jme,kms:kme)
   real    :: term3

   if (trace_use) call da_trace_entry("da_uvprho_to_w_lin")
   
   ! Computation to check for edge of domain:
   is = its
   ie = ite
   js = jts
   je = jte
   if (its == ids) is = ids+1
   if (ite == ide) ie = ide-1
   if (jts == jds) js = jds+1
   if (jte == jde) je = jde-1

   WZ(:,:,:) = 0.0
   ! Term 1.1: perturbed pressure advection along the basic wind
   do K=kts,kte
      do J=js,je
         do I=is,ie
            WZ(I,J,K)=WZ(I,J,K)-grid%xb%u(I,J,K)*(grid%xa%p(I+1,J,K)-grid%xa%p(I-1,J,K))* &
               grid%xb%coefx(I,J)
            WZ(I,J,K)=WZ(I,J,K)-grid%xb%v(I,J,K)*(grid%xa%p(I,J+1,K)-grid%xa%p(I,J-1,K))* &
               grid%xb%coefy(I,J)
         end do
      end do
   end do

   ! Term 1.2: Basic pressure advection along the perturbed wind

   do K=kts,kte
      do J=js,je
         do I=is,ie
            WZ(I,J,K)=WZ(I,J,K)-grid%xa%u(I,J,K)*(grid%xb%p(I+1,J,K)-grid%xb%p(I-1,J,K))* &
               grid%xb%coefx(I,J)
            WZ(I,J,K)=WZ(I,J,K)-grid%xa%v(I,J,K)*(grid%xb%p(I,J+1,K)-grid%xb%p(I,J-1,K))* &
               grid%xb%coefy(I,J)
         end do
      end do
   end do

   ! Dealing the laterial boundary because of the advection.
   ! boundary too simple? (It is the same as fill in interpf, fill can be used)

   if (its == ids) then
      do K=kts,kte
         do J=js,je
            WZ(its,J,K)=WZ(its+1,J,K)
         end do
      end do
   end if

   if (ite == ide) then
      do K=kts,kte
         do J=js,je
            WZ(ite,J,K)=WZ(ite-1,J,K)
         end do
      end do
   end if

   if (jts == jds) then
      do K=kts,kte
         do I=its, ite
            WZ(I,jts,K)=WZ(I,jts+1,K)       
         end do
      end do
   end if

   if (jte == jde) then
      do K=kts,kte
         do I=its, ite
            WZ(I,jte,K)=WZ(I,jte-1,K)
         end do
      end do
   end if

   ! Term 2.1: Divergence term from perturbed wind

   call da_uv_to_divergence(grid%xb, grid%xa%u, grid%xa%v, DIV)

   WZ(its:ite,jts:jte,kts:kte)=WZ(its:ite,jts:jte,kts:kte)-GAMMA*grid%xb%p(its:ite,jts:jte,kts:kte)*DIV(its:ite,jts:jte,kts:kte)

   ! Term 2.2: Divergence term from basic wind

   call da_uv_to_divergence(grid%xb, grid%xb%u, grid%xb%v, DIV)

   WZ(its:ite,jts:jte,kts:kte)=WZ(its:ite,jts:jte,kts:kte)-GAMMA*grid%xa%p(its:ite,jts:jte,kts:kte)*DIV(its:ite,jts:jte,kts:kte)

   ! Computation to check for edge of domain:
   is = its-1; ie = ite+1; js = jts-1; je = jte+1
   if (its == ids) is = ids; if (ite == ide) ie = ide
   if (jts == jds) js = jds; if (jte == jde) je = jde

   ! Term 3.1: Vertical integration of the perturbed mass divergence

   URHO(is:ie,js:je,kts:kte)=grid%xb%rho(is:ie,js:je,kts:kte)*grid%xa%u(is:ie,js:je,kts:kte)
   VRHO(is:ie,js:je,kts:kte)=grid%xb%rho(is:ie,js:je,kts:kte)*grid%xa%v(is:ie,js:je,kts:kte)

   call da_uv_to_divergence(grid%xb, URHO, VRHO, DIV)

   do J=jts,jte
      do I=its,ite
         TERM3=0.0

         do K=kte-1,kts,-1
            TERM3=TERM3+GRAVITY*(DIV(I,J,K+1)+DIV(I,J,K))*0.5 *(grid%xb%h(I,J,K+1)-grid%xb%h(I,J,K))
            WZ(I,J,K)=WZ(I,J,K)+TERM3
         end do
      end do
   end do

   ! Term 3.2: Vertical integration of the basic mass divergence

   URHO(is:ie,js:je,kts:kte)=grid%xa%rho(is:ie,js:je,kts:kte)*grid%xb%u(is:ie,js:je,kts:kte)
   VRHO(is:ie,js:je,kts:kte)=grid%xa%rho(is:ie,js:je,kts:kte)*grid%xb%v(is:ie,js:je,kts:kte)

   call da_uv_to_divergence(grid%xb, URHO, VRHO, DIV)

   do J=jts,jte
      do I=its,ite
         TERM3=0.0

         do K=kte-1,kts,-1
            TERM3=TERM3+GRAVITY*(DIV(I,J,K+1)+DIV(I,J,K))*0.5*(grid%xb%h(I,J,K+1)-grid%xb%h(I,J,K))
            WZ(I,J,K)=WZ(I,J,K)+TERM3
         end do
      end do
   end do

   ! Term 4: Derivative of basic vertical velocity with respect to z.

   do J=jts,jte
      do I=its,ite
         do K=kts,kte
            WZ(I,J,K)=WZ(I,J,K)-GAMMA*grid%xa%p(I,J,K)*(grid%xb%w(I,J,K+1)-grid%xb%w(I,J,K))/  &
               (grid%xb%hf(I,J,K+1)-grid%xb%hf(I,J,K))
         end do
      end do
   end do

   ! Divide by constant

   WZ(its:ite,jts:jte,kts:kte)=WZ(its:ite,jts:jte,kts:kte)/(GAMMA*grid%xb%p(its:ite,jts:jte,kts:kte))

   ! integration to calculate the vertical velocity 

   call da_w_adjustment_lin(grid%xb,grid%xa%w,WZ)

   do J=jts,jte
      do I=its,ite
         grid%xa%w(I,J,kte+1)=0.0
         do K=kte,kts,-1
            grid%xa%w(I,J,K)=grid%xa%w(I,J,K+1) + WZ(I,J,K)*(grid%xb%hf(I,J,K)-grid%xb%hf(I,J,K+1))
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_uvprho_to_w_lin")

end subroutine da_uvprho_to_w_lin


subroutine da_uvprho_to_w_adj(grid)

   !---------------------------------------------------------------------------
   ! Purpose: Calculates vertical velocity increments from Richardson's Eq.
   !
   ! Method: Richardson's Eq., which
   !         combines continuity Eq., thermodynamic Eq. and hrdrostatic Eq.
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout)     :: grid

   integer                       :: is, ie       ! 1st dim. end points.
   integer                       :: js, je       ! 2nd dim. end points.

   integer                       :: I,J,K

   real, dimension(ims:ime,jms:jme,kms:kme) :: URHO, VRHO
   real, dimension(ims:ime,jms:jme,kms:kme) :: DIV, WZ
   real                                     :: TERM3

   if (trace_use) call da_trace_entry("da_uvprho_to_w_adj")

   ! initialize to zero for some variables because of the adjoint requirements.

   WZ(:,:,:)   = 0.0
   URHO(:,:,:) = 0.0
   VRHO(:,:,:) = 0.0
   DIV(:,:,:)  = 0.0
   TERM3       = 0.0

   ! integration to calculate the vertical velocity

   do J=jts,jte
      do I=its,ite
         do K=kts,kte
            grid%xa%w(I,J,K+1)=grid%xa%w(I,J,K+1)+grid%xa%w(I,J,K)
            WZ(I,J,K)=grid%xa%w(I,J,K)*(grid%xb%hf(I,J,K)-grid%xb%hf(I,J,K+1))
            grid%xa%w(I,J,K)=0.0
         end do
         grid%xa%w(I,J,kte+1)=0.0
      end do
   end do

   call da_w_adjustment_adj(grid%xb,WZ)

   ! Divide by constant

   WZ(its:ite,jts:jte,kts:kte)=WZ(its:ite,jts:jte,kts:kte)/(GAMMA*grid%xb%p(its:ite,jts:jte,kts:kte))

   ! Term 4: Derivative of basic vertical velocity with respect to z.

   do J=jts,jte
      do I=its,ite
         do K=kte,kts,-1
            grid%xa%p(I,J,K)=grid%xa%p(I,J,K)-WZ(I,J,K)*GAMMA*  &
                     (grid%xb%w(I,J,K+1)-grid%xb%w(I,J,K))/  &
                     (grid%xb%hf(I,J,K+1)-grid%xb%hf(I,J,K))
         end do
      end do
   end do

   ! Term 3.2: Vertical integration of the basic mass divergence

   do J=jts,jte
      do I=its,ite
         do K=kts,kte-1
            TERM3=TERM3+WZ(I,J,K)
            DIV(I,J,K+1)=DIV(I,J,K+1)+  &
                      TERM3*GRAVITY*0.5*(grid%xb%h(I,J,K+1)-grid%xb%h(I,J,K))
            DIV(I,J,K)  =DIV(I,J,K)+  &
                      TERM3*GRAVITY*0.5*(grid%xb%h(I,J,K+1)-grid%xb%h(I,J,K))
         end do
        TERM3=0.0
      end do
   end do

   call da_uv_to_divergence_adj(grid, URHO,VRHO, DIV)

   ! Computation to check for edge of domain:
   if (test_transforms) then
      is = its-1; ie = ite+1; js = jts-1; je = jte+1
      if (its == ids) is = ids; if (ite == ide) ie = ide
      if (jts == jds) js = jds; if (jte == jde) je = jde
   else
      is = its
      ie = ite
      js = jts
      je = jte
   end if

   grid%xa%rho(is:ie,js:je,kts:kte)=grid%xa%rho(is:ie,js:je,kts:kte)+VRHO(is:ie,js:je,kts:kte)*grid%xb%v(is:ie,js:je,kts:kte)
   grid%xa%rho(is:ie,js:je,kts:kte)=grid%xa%rho(is:ie,js:je,kts:kte)+URHO(is:ie,js:je,kts:kte)*grid%xb%u(is:ie,js:je,kts:kte)

   URHO(:,:,:)=0.0
   VRHO(:,:,:)=0.0

   ! Term 3.1: Vertical integration of the perturbed mass divergence

   do J=jts,jte
      do I=its,ite
         do K=kts,kte-1
            TERM3=TERM3+WZ(I,J,K)
            DIV(I,J,K+1)=DIV(I,J,K+1)+  &
                      TERM3*GRAVITY*0.5*(grid%xb%h(I,J,K+1)-grid%xb%h(I,J,K))
            DIV(I,J,K)  =DIV(I,J,K)+  &
                      TERM3*GRAVITY*0.5*(grid%xb%h(I,J,K+1)-grid%xb%h(I,J,K))
         end do
         TERM3=0.0
      end do
   end do

   call da_uv_to_divergence_adj(grid, URHO,VRHO, DIV)

   grid%xa%v(is:ie,js:je,kts:kte)=grid%xa%v(is:ie,js:je,kts:kte)+VRHO(is:ie,js:je,kts:kte)*grid%xb%rho(is:ie,js:je,kts:kte)
   grid%xa%u(is:ie,js:je,kts:kte)=grid%xa%u(is:ie,js:je,kts:kte)+URHO(is:ie,js:je,kts:kte)*grid%xb%rho(is:ie,js:je,kts:kte)

   URHO(:,:,:)=0.0
   VRHO(:,:,:)=0.0

   ! Term 2.2: Divergence term from basic wind

   call da_uv_to_divergence(grid%xb, grid%xb%u,grid%xb%v, DIV)

   grid%xa%p(its:ite,jts:jte,kts:kte)=grid%xa%p(its:ite,jts:jte,kts:kte)-WZ(its:ite,jts:jte,kts:kte)*GAMMA*DIV(its:ite,jts:jte,kts:kte)
   
   ! Term 2.1: Divergence term from perturbed wind

   DIV(its:ite,jts:jte,kts:kte)=-WZ(its:ite,jts:jte,kts:kte)*GAMMA*grid%xb%p(its:ite,jts:jte,kts:kte)  ! DIV redefined

   call da_uv_to_divergence_adj(grid, grid%xa%u,grid%xa%v, DIV)

   ! Computation to check for edge of domain:
   is = its
   ie = ite
   js = jts
   je = jte
   if (its == ids) is = ids+1
   if (ite == ide) ie = ide-1
   if (jts == jds) js = jds+1
   if (jte == jde) je = jde-1

   ! Term 1.2: Basic pressure advection along the perturbed wind

   if (jte == jde) then
      j = jte
      do K=kts,kte
         do I=its, ite
            WZ(I,J-1,K)=WZ(I,J-1,K)+WZ(I,J,K)
         end do
      end do
   end if

   if (jts == jds) then
      j = jts
      do K=kts,kte
         do I=its, ite
            WZ(I,J+1,K)=WZ(I,J+1,K)+WZ(I,J,K)
         end do
      end do
   end if

   if (ite == ide) then
      i = ite
      do K=kts,kte
         do J=js,je
            WZ(I-1,J,K)=WZ(I-1,J,K)+WZ(I,J,K)
         end do
      end do
   end if

   if (its == ids) then
      i = its
      do K=kts,kte
         do J=js,je
            WZ(I+1,J,K)=WZ(I+1,J,K)+WZ(I,J,K)
         end do
      end do
   end if

   do K=kts,kte
      do J=js,je
         do I=is,ie
            grid%xa%v(I,J,K)=grid%xa%v(I,J,K)-WZ(I,J,K)*  &
               (grid%xb%p(I,J+1,K)-grid%xb%p(I,J-1,K))*grid%xb%coefy(I,J)
            grid%xa%u(I,J,K)=grid%xa%u(I,J,K)-WZ(I,J,K)*  &
               (grid%xb%p(I+1,J,K)-grid%xb%p(I-1,J,K))*grid%xb%coefx(I,J)
         end do
      end do
   end do

   !-------------------------------------------------------------------------
   ! Computation to check for edge of domain:
   ! This is only for adjoint, as we have to cross the processor boundary
   ! to get the contribution.

   is = its - 1
   ie = ite + 1
   js = jts - 1
   je = jte + 1

   grid%xp%v1z(its:ite, jts:jte, kts:kte) = wz(its:ite, jts:jte, kts:kte)

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_BAL_EQN_ADJ.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_BAL_EQN_ADJ_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE

   if (its == ids) then
      is = ids+1
   else
      wz(is, js:je, kts:kte) = grid%xp%v1z(is, js:je, kts:kte)
   end if

   if (ite == ide) then
      ie = ide-1
   else
      wz(ie, js:je, kts:kte) = grid%xp%v1z(ie, js:je, kts:kte)
   end if

   if (jts == jds) then
      js = jds+1
   else
      wz(is:ie, js, kts:kte) = grid%xp%v1z(is:ie, js, kts:kte)
   end if

   if (jte == jde) then
      je = jde-1
   else
      wz(is:ie, je, kts:kte) = grid%xp%v1z(is:ie, je, kts:kte)
   end if

   ! Term 1.1: Perturbed pressure advection along the basic wind

   do K=kts,kte
      do J=js,je
         do I=is,ie
            grid%xa%p(I,J+1,K)=grid%xa%p(I,J+1,K)-WZ(I,J,K)*grid%xb%v(I,J,K)*grid%xb%coefy(I,J)
            grid%xa%p(I,J-1,K)=grid%xa%p(I,J-1,K)+WZ(I,J,K)*grid%xb%v(I,J,K)*grid%xb%coefy(I,J)
            grid%xa%p(I+1,J,K)=grid%xa%p(I+1,J,K)-WZ(I,J,K)*grid%xb%u(I,J,K)*grid%xb%coefx(I,J)
            grid%xa%p(I-1,J,K)=grid%xa%p(I-1,J,K)+WZ(I,J,K)*grid%xb%u(I,J,K)*grid%xb%coefx(I,J)
         end do
      end do
   end do

   WZ(:,:,:) = 0.0

   if (trace_use) call da_trace_exit("da_uvprho_to_w_adj")

end subroutine da_uvprho_to_w_adj


subroutine da_pt_to_rho_adj(grid) 

   !---------------------------------------------------------------------------
   !  Purpose: Adjoint of da_pt_to_rho.
   !
   !  Assumptions: 1) Model level stored top down
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain), intent(inout)  :: grid

   integer                       :: i,j,k        ! Loop counter.

   integer                       :: is, ie, js, je

   real                          :: temp

   if (trace_use) call da_trace_entry("da_pt_to_rho_adj")
   
   is = its
   js = jts

   ie = ite
   je = jte

   if (test_transforms) then
      is = its-1
      js = jts-1

      ie = ite+1
      je = jte+1

      if (is < ids) is = ids
      if (js < jds) js = jds

      if (ie > ide) ie = ide
      if (je > jde) je = jde
   end if

   !---------------------------------------------------------------------------
   ! Calculate rho increments:
   !---------------------------------------------------------------------------

   do j=js, je
      do k=kts, kte
         do i=is, ie
            temp = grid%xa%rho(i,j,k) * grid%xb%rho(i,j,k)

            grid%xa%p(i,j,k) = grid%xa%p(i,j,k) + temp/grid%xb%p(i,j,k)
                           
            grid%xa%t(i,j,k) = grid%xa%t(i,j,k) - temp/grid%xb%t(i,j,k)

            grid%xa%rho(i,j,k) = 0.0
         end do        
      end do        
   end do        

   if (trace_use) call da_trace_exit("da_pt_to_rho_adj")

end subroutine da_pt_to_rho_adj


subroutine da_pt_to_rho_lin(grid)

   !---------------------------------------------------------------------------
   ! Purpose: Calculates density increments from pressure/temperature increments
   !
   ! Method:  Linearised ideal gas law: rho~/rho = p'/p - T'/T
   !
   ! Assumptions: 1) Model level stored top down.
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   if (trace_use) call da_trace_entry("da_pt_to_rho_lin")
   
   grid%xa % rho(its:ite,jts:jte,kts:kte) = grid%xb % rho(its:ite,jts:jte,kts:kte) * ( &
      grid%xa % p(its:ite,jts:jte,kts:kte) / grid%xb % p(its:ite,jts:jte,kts:kte) - &
      grid%xa % t(its:ite,jts:jte,kts:kte) / grid%xb % t(its:ite,jts:jte,kts:kte))                       


   if (trace_use) call da_trace_exit("da_pt_to_rho_lin")

end subroutine da_pt_to_rho_lin


subroutine da_tpq_to_rh( t, p, q, es, qs, rh )

   !---------------------------------------------------------------------------
   ! Purpose: Convert T/p/q to relative humidity rh.
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)  :: t, p, q
   real, intent(out) :: es, qs, rh

   if (trace_use_dull) call da_trace_entry("da_tpq_to_rh")

   !---------------------------------------------------------------------------
   ! [1.0] Calculate saturation specific humidity:
   !---------------------------------------------------------------------------

   call da_tp_to_qs( t, p, es, qs )
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate relative humidity:
   !---------------------------------------------------------------------------

   rh = 100.0 * q / qs

   if (trace_use_dull) call da_trace_exit("da_tpq_to_rh")

end subroutine da_tpq_to_rh


subroutine da_tpq_to_rh_lin(grid)

   !---------------------------------------------------------------------------
   !  Purpose: Convert T/pressure/q to relative humidity increments.
   !
   !  Method: r~ = r (q~/q - qs~/qs).
   !
   ! When q approaching to zero, the above formula its undefined. The
   ! general formula below must be used:
   ! 
   !  Method: r~ = 100 * (q~/qs - q*(qs~/qs)/qs))
   !             = 100 * q~/qs - (100*q/qs)*(qs~/qs)
   !             = 100 * q~/qs - rh * (qs~/qs) 
   !
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain),  intent(inout) ::grid
   
   real :: qs(its:ite,jts:jte,kts:kte)
   real :: es(its:ite,jts:jte,kts:kte)
   real :: qs_prime_over_qs(its:ite,jts:jte,kts:kte)

   if (trace_use_dull) call da_trace_entry("da_tpq_to_rh_lin")

   !---------------------------------------------------------------------------
   ! [1.0] Calculate saturation specific humidity ratio qs~/qs:
   !---------------------------------------------------------------------------
 
   call da_tp_to_qs_lin(grid, qs_prime_over_qs )

   !--------------------------------------------------------------------------
   ! [2.0] Culcalete background saturation specific humidity qs:
   !--------------------------------------------------------------------------

   call da_tp_to_qs1(grid, es, qs) 
   
   !---------------------------------------------------------------------------
   ! [3.0] Calculate relative humidity increment:
   !---------------------------------------------------------------------------

   grid%xa % rh(its:ite,jts:jte,kts:kte) = 100.0 * &
                                ( grid%xa % q(its:ite,jts:jte,kts:kte) / &
                                      qs(its:ite,jts:jte,kts:kte) ) - &
                                  grid%xb % rh(its:ite,jts:jte,kts:kte) * &
                                  qs_prime_over_qs(its:ite,jts:jte,kts:kte)

   if (trace_use_dull) call da_trace_exit("da_tpq_to_rh_lin")

end subroutine da_tpq_to_rh_lin


subroutine da_tpq_to_rh_lin1( t, p, es, rh, t_prime, p_prime, q_prime, rh_prime )

   !---------------------------------------------------------------------------
   !  Purpose: Convert T/pressure/q to relative humidity increments.
   !
   !  Method: r~ = r (q~/q - qs~/qs).
   !
   !  When q approaching to zero, the above formula is undefined. The
   !  general formula below must be used:
   ! 
   !  Method: r~ = 100 * (q~/qs) - rh*(qs~/qs)
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)  :: t        ! Temperature.
   real, intent(in)  :: p        ! Pressure.
   real, intent(in)  :: es       ! Saturation vapour pressure.
   real, intent(in)  :: rh       ! Relative Humidity.
   real, intent(in)  :: t_prime  ! Temperature increment.
   real, intent(in)  :: p_prime  ! Pressure increment.
   real, intent(in)  :: q_prime  ! Pressure increment.
   real, intent(out) :: rh_prime ! Pressure increment.
   
   real :: es1, qs  ! Saturation specific humidity.
   real :: qs_prime_over_qs ! qs~/qs.

   if (trace_use_dull) call da_trace_entry("da_tpq_to_rh_lin1")

   !---------------------------------------------------------------------------
   ! [1.0] Calculate saturation specific humidity ratio qs~/qs:
   !---------------------------------------------------------------------------

   call da_tp_to_qs_lin1( t, p, es, t_prime, p_prime, qs_prime_over_qs )
   
   !--------------------------------------------------------------------------
   ! [2.0] Culcalete background saturation specific humidity qs:
   !--------------------------------------------------------------------------

   call da_tp_to_qs( t, p, es1, qs) 
   
   !---------------------------------------------------------------------------
   ! [3.0] Calculate relative humidity increment:
   !---------------------------------------------------------------------------

   rh_prime = 100.0 * (q_prime / qs) - rh * qs_prime_over_qs

   if (trace_use_dull) call da_trace_exit("da_tpq_to_rh_lin1")

end subroutine da_tpq_to_rh_lin1


subroutine da_tprh_to_q_adj (grid)

   !---------------------------------------------------------------------------
   !  Purpose: Adjoint of da_tprh_to_q_adj.
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout)  :: grid

   integer :: is, ie   ! 1st dim. end points.
   integer :: js, je   ! 2nd dim. end points.
   integer :: ks, ke   ! 3rd dim. end points.
   integer :: i, j, k  ! Loop counter.
   real    :: qs_prime_over_qs(its:ite,jts:jte,kts:kte) ! Temp.

   if (trace_use) call da_trace_entry("da_tprh_to_q_adj")

   !---------------------------------------------------------------------------
   ! [1.0] initialise:
   !---------------------------------------------------------------------------

   is = its; ie = ite
   js = jts; je = jte
   ks = kts; ke = kte   

   if (test_transforms) then
      is = its-1
      js = jts-1

      ie = ite+1
      je = jte+1

      if ( is < ids ) is = ids
      if ( js < jds ) js = jds

      if ( ie > ide ) ie = ide
      if ( je > jde ) je = jde
   end if

   !---------------------------------------------------------------------------
   ! [2.0] Calculate relative humidity increment:
   !---------------------------------------------------------------------------

   do k = ks, ke
      do j = js, je
         do i = is, ie
            qs_prime_over_qs(i,j,k) = grid%xb % q(i,j,k) * grid%xa % q(i,j,k)

            grid%xa % rh(i,j,k) = grid%xa % rh(i,j,k) + qs_prime_over_qs(i,j,k) / &
                             grid%xb % rh(i,j,k)
         end do
      end do
   end do

   !---------------------------------------------------------------------------
   ! [2.0] Calculate saturation specific humidity ratio qs~/qs:
   !---------------------------------------------------------------------------

   call da_tp_to_qs_adj (grid, qs_prime_over_qs )

   if (trace_use) call da_trace_exit("da_tprh_to_q_adj")

end subroutine da_tprh_to_q_adj

   !subroutine da_tprh_to_q_adj( t, p, es, q, rh, &
   !                             t_prime, p_prime, rh_prime, q_prime, n )

   !---------------------------------------------------------------------------
   !  Purpose: Adjoint of da_tprh_to_q_adj.
   !---------------------------------------------------------------------------

   !   implicit none

   !   integer        i, n
   !   real           t, es, p, q, rh,t_prime, p_prime, rh_prime, q_prime  
   !   dimension      t       (n) ! Temperature.
   !   dimension      es      (n) ! Saturation vapour pressure.
   !   dimension      p       (n) ! Pressure.
   !   dimension      q       (n) ! Specific humidity.
   !   dimension      rh      (n) ! Relative Humidity.
   !   dimension      t_prime (n) ! Temperature increment.
   !   dimension      p_prime (n) ! Pressure increment.
   !   dimension      rh_prime(n) ! Pressure increment.
   !   dimension      q_prime (n) ! Pressure increment.

   !   real        temp, qs_prime_over_qs  ! Temporary storage.
   !   dimension   qs_prime_over_qs(n)     ! qs~/qs.

   !   do i = 1,n
   !   temp = q(i) * q_prime(i)

   !---------------------------------------------------------------------------
   !  [2.0] Calculate relative humidity increment:
   !---------------------------------------------------------------------------

   !   rh_prime(i) = rh_prime(i) + temp / rh(i)
   !   qs_prime_over_qs(i) = temp
   !   end do

   !---------------------------------------------------------------------------
   !  [1.0] Calculate saturation specific humidity ratio qs~/qs:
   !---------------------------------------------------------------------------

   !   call da_tp_to_qs_adj( t, p, es, t_prime, p_prime, qs_prime_over_qs, n )

   !end subroutine da_tprh_to_q_adj


subroutine da_tprh_to_q_adj1( t, p, es, q, rh, t_prime, p_prime, rh_prime, q_prime )

   !---------------------------------------------------------------------------
   !  Purpose: Adjoint of da_tprh_to_q_adj.
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)    :: t        ! Temperature.
   real, intent(in)    :: es       ! Saturation vapour pressure.
   real, intent(in)    :: p        ! Pressure.
   real, intent(in)    :: q        ! Specific humidity.
   real, intent(in)    :: rh       ! Relative Humidity.
   real, intent(inout) :: t_prime  ! Temperature increment.
   real, intent(inout) :: p_prime  ! Pressure increment.
   real, intent(inout) :: rh_prime ! Pressure increment.
   real, intent(in)    :: q_prime  ! Pressure increment.
   
   real :: temp     ! Temporary storage.
   real :: qs_prime_over_qs ! qs~/qs.

   if (trace_use) call da_trace_entry("da_tprh_to_q_adj1")

   temp = q * q_prime

   !---------------------------------------------------------------------------
   ! [2.0] Calculate relative humidity increment:
   !---------------------------------------------------------------------------

   rh_prime = rh_prime + temp / rh
   qs_prime_over_qs = temp

   !---------------------------------------------------------------------------
   ! [1.0] Calculate saturation specific humidity ratio qs~/qs:
   !---------------------------------------------------------------------------

   call da_tp_to_qs_adj1 (t, p, es, t_prime, p_prime, qs_prime_over_qs)

   if (trace_use) call da_trace_exit("da_tprh_to_q_adj1")

end subroutine da_tprh_to_q_adj1


subroutine da_tprh_to_q_lin(grid)

   !---------------------------------------------------------------------------
   ! Purpose: Convert T/pressure/rh to specific humidity increments.
   !
   ! Method: q~ = q (rh~/rh + qs~/qs)
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain), intent(inout) :: grid
   
   real :: qs_prime_over_qs(its:ite,jts:jte,kts:kte) ! qs~/qs.

   if (trace_use) call da_trace_entry("da_tprh_to_q_lin")

   !---------------------------------------------------------------------------
   ! [1.0] Calculate saturation specific humidity ratio qs~/qs:
   !---------------------------------------------------------------------------

   call da_tp_to_qs_lin( grid, qs_prime_over_qs )
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate specific humidity increment:
   !---------------------------------------------------------------------------

   grid%xa % q(its:ite,jts:jte,kts:kte) = grid%xb % q(its:ite,jts:jte,kts:kte) * &
                               ( grid%xa % rh(its:ite,jts:jte,kts:kte) / &
                                 grid%xb % rh(its:ite,jts:jte,kts:kte) + &
                                 qs_prime_over_qs(its:ite,jts:jte,kts:kte) )

   if (trace_use) call da_trace_exit("da_tprh_to_q_lin")

end subroutine da_tprh_to_q_lin


subroutine da_tprh_to_q_lin1( t, p, es, q, rh, t_prime, p_prime, rh_prime, q_prime )

   !---------------------------------------------------------------------------
   !  Purpose: Convert T/pressure/rh to specific humidity increments.
   !
   !  Method: q~ = q (rh~/rh + qs~/qs)
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)  :: t        ! Temperature.
   real, intent(in)  :: p        ! Pressure.
   real, intent(in)  :: es       ! Saturation vapour pressure.
   real, intent(in)  :: q        ! Specific humidity.
   real, intent(in)  :: rh       ! Relative Humidity.
   real, intent(in)  :: t_prime  ! Temperature increment.
   real, intent(in)  :: p_prime  ! Pressure increment.
   real, intent(in)  :: rh_prime ! Pressure increment.
   real, intent(out) :: q_prime  ! Pressure increment.
   
   real :: qs_prime_over_qs ! qs~/qs.

   if (trace_use_dull) call da_trace_entry("da_tprh_to_q_lin1")

   !---------------------------------------------------------------------------
   ! [1.0] Calculate saturation specific humidity ratio qs~/qs:
   !---------------------------------------------------------------------------

   call da_tp_to_qs_lin1( t, p, es, t_prime, p_prime, qs_prime_over_qs )
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate specific humidity increment:
   !---------------------------------------------------------------------------

   q_prime = q * ( rh_prime / rh + qs_prime_over_qs )

   if (trace_use_dull) call da_trace_exit("da_tprh_to_q_lin1")

end subroutine da_tprh_to_q_lin1


subroutine da_tp_to_qs( t, p, es, qs)

   !---------------------------------------------------------------------------
   ! Purpose: Convert T/p to saturation specific humidity.
   !
   !  Method: qs = es_alpha * es / ( p - ( 1 - rd_over_rv ) * es ).
   !          use Rogers & Yau (1989) formula: es = a exp( bTc / (T_c + c) )
   !--------------------------------------------------------------------------

   implicit none

   real, intent(in)  :: t, p
   real, intent(out) :: es, qs
   
   real              :: t_c              ! T in degreesC.

   if (trace_use_dull) call da_trace_entry("da_tp_to_qs")

   !---------------------------------------------------------------------------
   ! [1.0] initialise:
   !---------------------------------------------------------------------------
   
   t_c = t - t_kelvin
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate saturation vapour pressure:
   !---------------------------------------------------------------------------

   es = es_alpha * exp( es_beta * t_c / ( t_c + es_gamma ) )
    
   !---------------------------------------------------------------------------
   ! [3.0] Calculate saturation specific humidity:
   !---------------------------------------------------------------------------

   qs = rd_over_rv * es / ( p - rd_over_rv1 * es )

   if (trace_use_dull) call da_trace_exit("da_tp_to_qs")

end subroutine da_tp_to_qs


subroutine da_tp_to_qs1 (grid, es, qs)

   !---------------------------------------------------------------------------
   !  Purpose: Convert T/p to saturation specific humidity.
   !
   !  Method: qs = es_alpha * es / ( p - ( 1 - rd_over_rv ) * es ).
   !          use Rogers & Yau (1989) formula: es = a exp( bTc / (T_c + c) ).
   !
   !  This da_tp_to_qs1 was added and called by the corrected subroutine
   !       da_tpq_to_rh_lin.
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(in)  :: grid
   real,          intent(out) :: es(its:ite,jts:jte,kts:kte)
   real,          intent(out) :: qs(its:ite,jts:jte,kts:kte)

   integer :: i, j, k      ! Loop counters.
   real    :: t_c          ! Working variable.

   if (trace_use_dull) call da_trace_entry("da_tp_to_qs1")

   !---------------------------------------------------------------------------
   ! [1.0] initialise:
   !---------------------------------------------------------------------------
    

   do k = kts, kte
      do j = jts, jte
         do i = its, ite

            !------------------------------------------------------------------
            ! [1.0] initialise:
            !------------------------------------------------------------------

            t_c = grid%xb % t(i,j,k) - t_kelvin
   
            !------------------------------------------------------------------
            ! [2.0] Calculate saturation vapour pressure:
            !------------------------------------------------------------------

            es(i,j,k) = es_alpha * exp( es_beta * t_c / ( t_c + es_gamma ) )
   
            !------------------------------------------------------------------
            ! [3.0] Calculate saturation specific humidity:
            !------------------------------------------------------------------

            qs(i,j,k) = rd_over_rv * es(i,j,k) / &
                     (grid%xb % p(i,j,k) - rd_over_rv1 * es(i,j,k))

         end do
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_tp_to_qs1")

end subroutine da_tp_to_qs1


subroutine da_tp_to_qs_adj (grid, qs_prime_over_qs )

   !---------------------------------------------------------------------------
   ! Purpose: Adjoint of da_tp_to_qs_lin
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout)  :: grid
   real,          intent(in)     :: qs_prime_over_qs(its:ite,jts:jte,kts:kte)

   integer :: is, ie       ! 1st dim. end points.
   integer :: js, je       ! 2nd dim. end points.
   integer :: ks, ke       ! 3rd dim. end points.
   integer :: i, j, k      ! Loop counters.
   real    :: temp         ! Temporary array.
   real    :: es_prime_over_es ! Sat Vap pressure ratio.

   if (trace_use) call da_trace_entry("da_tp_to_qs_adj")

   !---------------------------------------------------------------------------
   ! [1.0] initialise:
   !---------------------------------------------------------------------------

   is = its; ie = ite
   js = jts; je = jte
   ks = kts; ke = kte      

   if ( test_transforms ) then
      is = its-1
      js = jts-1

      ie = ite+1
      je = jte+1

      if ( is < ids ) is = ids
      if ( js < jds ) js = jds

      if ( ie > ide ) ie = ide
      if ( je > jde ) je = jde
   end if

   !---------------------------------------------------------------------------
   ! [3.0] Calculate saturation specific humidity increment:
   !---------------------------------------------------------------------------

   do k = ks, ke
      do j = js, je
         do i = is, ie

            temp = qs_prime_over_qs(i,j,k) / &
                   ( grid%xb % p(i,j,k) - rd_over_rv1 * grid%xb % es(i,j,k) )
   
            es_prime_over_es = temp * grid%xb % p(i,j,k)

            grid%xa % p(i,j,k) = grid%xa % p(i,j,k) - temp
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate saturation vapour pressure increment:
   !---------------------------------------------------------------------------

            temp = grid%xb % t(i,j,k) + es_gammakelvin

            grid%xa % t(i,j,k) = grid%xa % t(i,j,k) + es_gammabeta * es_prime_over_es / &
                            ( temp * temp )
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_tp_to_qs_adj")

end subroutine da_tp_to_qs_adj

   !subroutine da_tp_to_qs_adj( t, p, es, t_prime, p_prime, &
   !                            qs_prime_over_qs, n )

   !---------------------------------------------------------------------------
   ! Purpose: Adjoint of da_tp_to_qs_lin
   !---------------------------------------------------------------------------

   ! implicit none

   ! integer      i, n
   ! real         t, p, es, t_prime, p_prime, qs_prime_over_qs
   ! dimension    t               (n) ! Temperature.
   ! dimension    p               (n) ! Pressure.
   ! dimension    es              (n) ! Sat. vapour pressure.
   ! dimension    t_prime         (n) ! Temperature increment.
   ! dimension    p_prime         (n) ! Pressure increment.
   ! dimension    qs_prime_over_qs(n) ! qs~/qs.

   ! real         temp             ! Temporary storage.
   ! real         es_prime_over_es ! es~/es
   !    
   ! do i = 1,n
      !------------------------------------------------------------------------
      !    [3.0] Calculate saturation specific humidity increment:
      !------------------------------------------------------------------------

      ! temp = qs_prime_over_qs(i) / ( p(i) - rd_over_rv1 * es(i) )

      ! es_prime_over_es = temp * p(i)

      ! p_prime(i) = p_prime(i) - temp

      !------------------------------------------------------------------------
      ! [2.0] Calculate saturation vapour pressure increment:
      !------------------------------------------------------------------------

      ! temp = t(i) + es_gammakelvin

      ! t_prime(i) = t_prime(i) + es_gammabeta * es_prime_over_es / ( temp * temp )
   ! end do

   ! end subroutine da_tp_to_qs_adj


subroutine da_tp_to_qs_adj1( t, p, es, t_prime, p_prime, qs_prime_over_qs )

   !---------------------------------------------------------------------------
   !  Purpose: Adjoint of da_tp_to_qs_lin.
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(in)    :: t                ! Temperature.
   real, intent(in)    :: p                ! Pressure.
   real, intent(in)    :: es               ! Sat. vapour pressure.
   real, intent(inout) :: t_prime          ! Temperature increment.
   real, intent(inout) :: p_prime          ! Pressure increment.
   real, intent(in)    :: qs_prime_over_qs ! qs~/qs.
   
   real :: temp             ! Temporary storage.
   real :: es_prime_over_es ! es~/es

   if (trace_use) call da_trace_entry("da_tp_to_qs_adj1")
      
   !---------------------------------------------------------------------------
   ! [3.0] Calculate saturation specific humidity increment:
   !---------------------------------------------------------------------------

   temp = qs_prime_over_qs / ( p - rd_over_rv1 * es )
   
   es_prime_over_es = temp * p

   p_prime = p_prime - temp
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate saturation vapour pressure increment:
   !---------------------------------------------------------------------------

   temp = t + es_gammakelvin

   t_prime = t_prime + es_gammabeta * es_prime_over_es / ( temp * temp )

   if (trace_use) call da_trace_exit("da_tp_to_qs_adj1")

end subroutine da_tp_to_qs_adj1


subroutine da_tp_to_qs_lin(grid, qs_prime_over_qs )

   !---------------------------------------------------------------------------
   !  Purpose: Convert es/p/es_prime to saturation specific humidity increment.
   !
   !  Method: qs~ = qs * ( p es'/es - p' ) / ( p - (1-rd_over_rv) es ).
   !          use Rogers & Yau (1989) formula: es = a exp( bTc / (T_c + c) ).

   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid
   real,          intent(out)   :: qs_prime_over_qs(its:ite,jts:jte,kts:kte)

   integer :: i, j, k      ! Loop counters.
   real    :: temp         ! Temporary array.
   real    :: es_prime_over_es ! Sat Vap pressure ratio.

   if (trace_use_dull) call da_trace_entry("da_tp_to_qs_lin")

   do k = kts, kte
      do j = jts, jte
         do i = its, ite
            temp = grid%xb % t(i,j,k) + es_gammakelvin
            !-----------------------------------------------------------------
            ! [2.0] Calculate saturation vapour pressure increment:
            !-----------------------------------------------------------------

            es_prime_over_es = es_gammabeta * grid%xa % t(i,j,k) / ( temp * temp )

            !-----------------------------------------------------------------
            ! [3.0] Calculate saturation specific humidity increment:
            !-----------------------------------------------------------------

            qs_prime_over_qs(i,j,k) = ( grid%xb % p(i,j,k) * es_prime_over_es - &
                                        grid%xa % p(i,j,k) ) / &
                                      ( grid%xb % p(i,j,k) - rd_over_rv1 * &
                                        grid%xb % es(i,j,k) )
         end do
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_tp_to_qs_lin")

end subroutine da_tp_to_qs_lin


subroutine da_tp_to_qs_lin1( t, p, es, t_prime, p_prime, qs_prime_over_qs )

   !---------------------------------------------------------------------------
   ! Purpose: Convert es/p/es_prime to saturation specific humidity increment.
   !
   !  Method: qs~ = qs * ( p es'/es - p' ) / ( p - (1-rd_over_rv) es ).
   !          use Rogers & Yau (1989) formula: es = a exp( bTc / (T_c + c) ).
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(in)  :: t                ! Temperature.
   real, intent(in)  :: p                ! Pressure.
   real, intent(in)  :: es               ! Sat. vapour pressure.
   real, intent(in)  :: t_prime          ! Temperature increment.
   real, intent(in)  :: p_prime          ! Pressure increment.
   real, intent(out) :: qs_prime_over_qs ! qs~/qs.
   
   real :: temp           ! Temporary value.
   real :: es_prime_over_es ! es~/es

   if (trace_use_dull) call da_trace_entry("da_tp_to_qs_lin1")

   !---------------------------------------------------------------------------
   ! [1.0] initialise:
   !---------------------------------------------------------------------------

   temp = t + es_gammakelvin
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate saturation vapour pressure increment:
   !---------------------------------------------------------------------------

   es_prime_over_es = es_gammabeta * t_prime / ( temp * temp )

   !---------------------------------------------------------------------------
   ! [3.0] Calculate saturation specific humidity increment:
   !---------------------------------------------------------------------------

   qs_prime_over_qs = (p * es_prime_over_es - p_prime) / (p - rd_over_rv1 * es)

   if (trace_use_dull) call da_trace_exit("da_tp_to_qs_lin1")


end subroutine da_tp_to_qs_lin1


subroutine da_trh_to_td (grid)

   !---------------------------------------------------------------------
   !
   !                       function f_td_from_rh
   !                     **************************
   !
   !  purpose:
   !  -------
   !     compute dew point from temperature and relative humidity
   !
   !   method:
   !   ------
   !     invert the relation
   !
   !     rh = 100.0 * exp (l_over_rv * (1.0/t - 1.0/td))
   !
   !   input:
   !   -----
   !      t_k:   temperature       in k
   !      rh:    relative humidity in %
   !
   !   output:
   !   ------
   !      td:    dew point in k
   !
   !   references:
   !   -----------
   !    R. R. Rogers and M. K. Yau, 1989: a short course in cloud physics,
   !                                   3nd edition, pergamon press, page 14-19.
   !
   !   verification set:
   !   -----------------
   !    t_k  = 268.15 k,  
   !    td_k = 262.55 k
   !    rh   = 65 %, 
   !    p_pa = 80000  pa, 
   !    qv   = 2.11e-03 kg/kg,
   !
   !  modifications:
   !   ------------
   !    parallel implementation. -al bourgeoits
   ! 
   !-------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   integer :: i, j, k, ij

   real    :: invdifftd, invtd

   if (trace_use_dull) call da_trace_entry("da_trh_to_td")

   !$OMP PARALLEL DO &
   !$OMP PRIVATE( ij, i, j, k, invdifftd, invtd )
   do ij = 1 , grid%num_tiles

   do k=kts,kte
      do j=grid%j_start(ij), grid%j_end(ij)
         do i=its,ite
            if (grid%xb%rh(i,j,k) < 10.0) then
               grid%xb%rh(i,j,k) = 10.0
            else if (grid%xb%rh(i,j,k) > 105.0) then
               grid%xb%rh(i,j,k) = 105.0
            end if

            invdifftd = log (grid%xb%rh(i,j,k)/100.0) / l_over_rv

            invtd = 1/grid%xb%t(i,j,k)  - invdifftd

            grid%xb%td(i,j,k)  = 1.0 / invtd

            if (grid%xb%td(i,j,k) > grid%xb%t(i,j,k)) &
               grid%xb%td(i,j,k) = grid%xb%t(i,j,k)
         end do
      end do
   end do

   end do
   !$OMP END PARALLEL DO

   if (trace_use_dull) call da_trace_exit("da_trh_to_td")

end subroutine da_trh_to_td


subroutine da_tpq_to_slp ( t, q, p, terr, psfc, slp)

   !-----------------------------------------------------------------------
   ! purpose:  computes sea level pressure from the rule                
   !              t1/t2=(p1/p2)**(gamma*r/g).                              
   !                                                                       
   !     input       t        temperature
   !                 q        mixing ratio
   !                 p        pressure
   !                 terr     terrain
   !                 psfc     surface pressure
   !                                                                       
   !     output      slp      sea level pressure    
   !-----------------------------------------------------------------------        

   implicit none

   real, intent(in)    :: terr, psfc
   real, intent(in)    :: t(kms:kme)
   real, intent(in)    :: q(kms:kme)
   real, intent(in)    :: p(kms:kme)
   real, intent(inout) :: slp

   integer         :: k, klo, khi
   real            :: pl, t0, ts, xterm,tlo, thi, tl
                                          
   real, parameter :: gamma  = 6.5e-3
   real, parameter :: tc     = t_kelvin+17.5
   real, parameter :: pconst = 10000.0
   real, parameter :: eps    = 0.622

   if (trace_use) call da_trace_entry("da_tpq_to_slp")
                                                                       
   ! sea level pressure                                            
                                                                         
   xterm=gamma* gas_constant / gravity                                                   
                                                                       
   ! compute pressure at pconst mb above surface (pl)              
                                                                        
   if (terr <= 0.0) then
      slp = psfc
      if (trace_use) call da_trace_exit("da_tpq_to_slp")
      return
   end if

   pl  = psfc - pconst                                        
   klo = 0

   ! find 2 levels on sigma surfaces surrounding pl at each i,j    

   do k=kts, kte-1
      if ((p(k) >= pl) .and. (p(k+1) < pl)) then
         khi = k+1
         klo = k
         exit
      end if
   end do

   if (klo < 1) then                                      
      write(unit=message(1),fmt='(a,f11.3,a)') &
         'error finding pressure level ',pconst,' mb above the surface'
      write(unit=message(2),fmt='(a,f11.3,2x,a,f11.3)') 'pl=',pl,'  psfc=',psfc
      call da_error("da_tpq_to_slp.inc",63,message(1:2))                                
   end if                                                         

   ! get temperature at pl (tl), extrapolate t at surface (ts)     
   ! and t at sea level (t0) with 6.5 k/km lapse rate              

   tlo=t(klo) * (eps+q(klo))/(eps*(1.0+q(klo)))
   thi=t(khi) * (eps+q(khi))/(eps*(1.0+q(khi)))
   tl=thi-(thi-tlo)*log(pl/p(khi)) &
                      /log(p(klo)/p(khi))               
   ts=tl*(psfc/pl)**xterm                           
   t0=ts +gamma*terr

   ! correct sea level temperature if too hot                      

   if ( t0 >= tc ) then
      if ( ts <= tc ) then
        t0 = tc
      else
        t0 = tc-0.005*(ts-tc)**2
      end if
   end if

   ! compute sea level pressure                                    

   slp=psfc*exp(2.0*gravity*terr/(gas_constant*(ts+t0)))

   if (trace_use) call da_trace_exit("da_tpq_to_slp")

end subroutine da_tpq_to_slp


subroutine da_wrf_tpq_2_slp (grid)

   !---------------------------------------------------------------------
   ! Purpose: computes sea level pressure from the rule                
   !          t1/t2=(p1/p2)**(gamma*r/g).                              
   !                                                                     
   ! input       t        temperature                cross    3d       
   !             q        mixing ratio               cross    2d       
   !             p        pressure perturbation      cross    2d       
   !             terr     terrain                    cross    2d       
   !             psfc     surface pressure           cross    2d       
   !                                                                       
   ! output      slp      sea level pressure         cross    2d    
   !---------------------------------------------------------------------

   implicit none   
              

   type (domain), intent(inout) :: grid

   integer              :: I, J, K, KLO, KHI, ij
   real                 :: PL, T0, TS, XTERM,    &
                          TLO, THI, TL

   real, parameter      :: GAMMA = 6.5E-3,  &
                          TC=t_kelvin+17.5,  &
                          PCONST=10000.0 ,  &
                          EPS   = 0.622

   if (trace_use_dull) call da_trace_entry("da_wrf_tpq_2_slp")

                                                                     
   ! sea level pressure                                            
                                                                         
   xterm=gamma* gas_constant / gravity                                                   
                                                                       
   ! compute pressure at pconst mb above surface (pl)              
                                                                          

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (ij, i, j, k, PL, klo, khi, tlo, thi, tl, ts, t0)
   do ij = 1, grid%num_tiles

   j_loop: do j=grid%j_start(ij), grid%j_end(ij)
      i_loop: do i=its, ite
         if (grid%xb%terr(i,j) <= 0.0) then
            grid%xb%slp(i,j) = grid%xb%psfc(i,j)
            cycle i_loop
         end if

         PL  = grid%xb%psfc(i,j) - PCONST                                        
         klo = 0

         ! Find 2 levels on sigma surfaces surrounding PL at each I,J    

         k_loop: do k=kts, kte-1
            if ((grid%xb%p(i,j,k) >= pl) .and. (grid%xb%p(i,j,k+1) < pl)) then
               khi = k+1
               klo = k
               exit k_loop
            end if
         end do k_loop

         if (klo < 1) then 
            write(unit=message(1),fmt='(A,F11.3,A)') &
               'Cannot find pressure level ',PCONST,' Pa above the surface'
            write(unit=message(2),fmt='(A,F11.3,2X,A,F11.3)') &
               'PL=',PL,'  PSFC=',grid%xb%psfc(i,j)
            call da_error("da_wrf_tpq_2_slp.inc",69,message(1:2))                                               
         end if                                                         

         ! get temperature at pl (tl), extrapolate t at surface (ts)     
         ! and t at sea level (t0) with 6.5 k/km lapse rate              

         tlo=grid%xb%t(i,j,klo) * (eps+grid%xb%q(i,j,klo))/(eps*(1.0+grid%xb%q(i,j,klo)))
         thi=grid%xb%t(i,j,khi) * (eps+grid%xb%q(i,j,khi))/(eps*(1.0+grid%xb%q(i,j,khi)))
         tl=thi-(thi-tlo)*log(pl/grid%xb%p(i,j,khi)) &
                      /log(grid%xb%p(i,j,klo)/grid%xb%p(i,j,khi))               
         ts=tl*(grid%xb%psfc(i,j)/pl)**xterm                           
         t0=ts +gamma*grid%xb%terr(i,j)

         ! correct sea level temperature if too hot                      

         if (t0 >= tc) then
            if (ts <= tc) then
               t0 = tc
            else
               t0 = tc-0.005*(ts-tc)**2
            end if
         end if

         ! compute sea level pressure                                    

         grid%xb%slp(i,j)=grid%xb%psfc(i,j)*EXP(2.0*gravity*grid%xb%terr(i,j)/ &
            (gas_constant*(TS+T0)))
      end do i_loop
   end do j_loop

   end do 
   !$OMP END PARALLEL DO  

   if (trace_use_dull) call da_trace_exit("da_wrf_tpq_2_slp")

end subroutine da_wrf_tpq_2_slp


subroutine da_tpq_to_slp_lin ( T, Q, P, TERR, PSFC, T9, Q9, P9, PSFC9, SLP9)

   !-----------------------------------------------------------------------
   ! Purpose: computes sea level pressure from the rule                
   !              t1/t2=(p1/p2)**(gamma*r/g).                              
   !                                                                       
   !     input       t        temperature
   !                 q        mixing ratio
   !                 p        pressure
   !                 terr     terrain
   !                 psfc     surface pressure
   !                                                                       
   !     output      slp      sea level pressure
   !-----------------------------------------------------------------------

   implicit none
              
   real, intent(in)    :: terr, psfc, psfc9
   real, intent(in)    :: t(kms:kme)
   real, intent(in)    :: q(kms:kme)
   real, intent(in)    :: p(kms:kme)
   real, intent(in)    :: t9(kms:kme)
   real, intent(in)    :: q9(kms:kme)
   real, intent(in)    :: p9(kms:kme)
   ! real                :: slp
   real, intent(out)   :: slp9

   integer :: k, klo, khi
   real    :: pl, t0, ts, xterm,tlo, thi, tl
   real    :: pl9,t09,ts9,tlo9,thi9,tl9,coef1,coef2

   real, parameter      :: gamma  = 6.5e-3
   real, parameter      :: tc     = t_kelvin+17.5
   real, parameter      :: pconst = 10000.0
   real, parameter      :: eps    = 0.622

   if (trace_use) call da_trace_entry("da_tpq_to_slp_lin")
                                                                          
   !     ... sea level pressure                                            
                                                                          
   xterm=gamma* gas_constant / gravity                                                   
                                                                          
   ! compute pressure at pconst mb above surface (pl)                                                                                     

   if (terr <= 0.0) then
      slp9 = psfc9
      ! slp = psfc
      if (trace_use) call da_trace_exit("da_tpq_to_slp_lin")
      return
   end if


   pl9  = psfc9 
   pl  = psfc - pconst                                        
   klo = 0

   ! find 2 levels on sigma surfaces surrounding pl at each i,j    

   k_loop: do k=kts, kte-1
      if ((p(k) >= pl) .and. (p(k+1) < pl)) then
         khi = k+1
         klo = k
         exit k_loop
      end if
   end do k_loop

   if (klo < 1) then                                      
      write(unit=message(1),fmt='(A,F11.3,A)') &
           'error finding pressure level ',pconst,' mb above the surface'
      write(unit=message(2),fmt='(A,F11.3,2X,A,F11.3)') 'PL=',PL,'  PSFC=',psfc
      call da_error("da_tpq_to_slp_lin.inc",71,message(1:2))                                               
   end if                                                         

   ! get temperature at pl (tl), extrapolate t at surface (ts)     
   ! and t at sea level (t0) with 6.5 k/km lapse rate              

   tlo9=t9(klo) * (eps+q(klo))/(eps*(1.0+q(klo))) + &
        q9(klo)*t(klo)*(1.0-eps)/(eps*(1.0+q(klo))**2)
   tlo=t(klo) * (eps+q(klo))/(eps*(1.0+q(klo)))
   thi9=t9(khi) * (eps+q(khi))/(eps*(1.0+q(khi)))+   &
        q9(khi)*t(khi)*(1.0-eps)/(eps*(1.0+q(khi))**2)
   thi=t(khi) * (eps+q(khi))/(eps*(1.0+q(khi)))
   coef1=alog(pl/p(khi))
   coef2=alog(p(klo)/p(khi))
   tl9=(1.0-coef1/coef2)*thi9+coef1/coef2*tlo9       &
       -(thi-tlo)/(coef2*pl)*pl9                 &
       +((thi-tlo)/(coef2*p(khi))*(1.0-coef1/coef2))*p9(khi)   &
       +(thi-tlo)*coef1/(coef2*coef2*p(klo))*p9(klo)
   tl=thi-(thi-tlo)*coef1/coef2
   ts9=tl9*(psfc/pl)**xterm+psfc9*xterm*(tl/pl)*(psfc/pl)**  &
       (xterm-1)-pl9*xterm*(tl*psfc/(pl*pl))*(psfc/pl)**(xterm-1)
   ts=tl*(psfc/pl)**xterm                           
   t09=ts9
   t0=ts +gamma*terr

   ! correct sea level temperature if too hot                      

   if ( t0 >= tc ) then
      if ( ts <= tc ) then
         t09 = 0.0
         t0 = tc
      else
         t09 = -0.01*(ts-tc)*ts9
         t0 = tc-0.005*(ts-tc)**2
      end if
   end if

   ! compute sea level pressure                                    

   slp9=psfc9*exp(2.0*gravity*terr/(gas_constant*(ts+t0)))  &
          -psfc*exp(2.0*gravity*terr/(gas_constant*(ts+t0)))*  &
          2.0*gravity*terr/(gas_constant*(ts+t0)**2)*(ts9+t09)
   ! slp=psfc*exp(2.0*gravity*terr/(gas_constant*(ts+t0)))

   if (trace_use) call da_trace_exit("da_tpq_to_slp_lin")

end subroutine da_tpq_to_slp_lin


subroutine da_tpq_to_slp_adj (T, Q, P, TERR, PSFC, T9, Q9, P9, PSFC9, SLP9)

   !-----------------------------------------------------------------------
   ! Purpose: computes sea level pressure from the rule                
   !              t1/t2=(p1/p2)**(gamma*r/g).                              
   !                                                                       
   !     input       t        temperature
   !                 q        mixing ratio
   !                 p        pressure
   !                 terr     terrain
   !                 psfc     surface pressure
   !                                                                       
   !     output      slp      sea level pressure
   !              
   !-----------------------------------------------------------------------

   implicit none

   real,              intent(in)  :: TERR, PSFC
   real,              intent(in)  :: T(kms:kme), Q(kms:kme), P(kms:kme)
   real,              intent(out) :: T9(kms:kme), Q9(kms:kme), P9(kms:kme)
   real,              intent(in)  :: SLP9
   real,              intent(out) :: PSFC9

   ! real :: SLP

   integer              :: K, KLO, KHI
   real                 :: PL, T0, TS, XTERM,TLO, THI, TL
   real                 :: PL9,T09,TS9,TLO9,THI9,TL9,COEF1,COEF2
   real                 :: T08,TS8

   real, parameter      :: GAMMA  = 6.5E-3
   real, parameter      :: TC     = t_kelvin+17.5
   real, parameter      :: PCONST = 10000.0
   real, parameter      :: EPS    = 0.622

   if (trace_use) call da_trace_entry("da_tpq_to_slp_adj")
                                                                        
   ! sea level pressure                                            
                                                                          
   xterm=gamma* gas_constant / gravity
   
   ! compute pressure at pconst mb above surface (pl)              
                                                                          
   T9(:) = 0.0
   Q9(:) = 0.0
   P9(:) = 0.0

   if (terr <= 0.0) then
      psfc9=slp9
      if (trace_use) call da_trace_exit("da_tpq_to_slp_adj")
      return
   end if

   PL  = psfc - PCONST                                        
   klo = 0

   ! find 2 levels on sigma surfaces surrounding pl at each i,j    

   do k=kts, kte-1
      if ((p(k) >= pl) .and. (p(k+1) < pl)) then
         khi = k+1
         klo = k
         exit
      end if
   end do

   if (klo < 1) then                                      
      write(unit=message(1),fmt='(A,F11.3,A)') &
         'Error finding pressure level ',PCONST,' Mb above the surface'
      write(unit=message(2),fmt='(A,F11.3,2X,A,F11.3)') &
         'PL=',PL,'  PSFC=',psfc
      call da_error("da_tpq_to_slp_adj.inc",73,message(1:2))
   end if                                                         

   ! get temperature at pl (tl), extrapolate t at surface (ts)     
   ! and t at sea level (t0) with 6.5 k/km lapse rate              

   TLO   = t(KLO) * (EPS+q(KLO))/(EPS*(1.0+q(KLO)))
   THI   = t(KHI) * (EPS+q(KHI))/(EPS*(1.0+q(KHI)))
   COEF1 = ALOG(PL/p(KHI))
   COEF2 = ALOG(p(KLO)/p(KHI))
   TL    = THI-(THI-TLO)*COEF1/COEF2
   TS    = TL*(psfc/PL)**XTERM                           
   TS8   = TS
   T0    = TS +GAMMA*terr
   T08   = T0

   ! Correct sea level temperature if too hot                      

   if (t0 >= tc) then
      if (ts <= tc) then
         t0 = tc
      else
         t0 = tc-0.005*(ts-tc)**2
      end if
   end if
   
   psfc9 = slp9*EXP(2.0*gravity*TERR/(gas_constant*(TS+T0)))
   TS9 = -slp9*psfc*EXP(2.0*gravity*TERR/(gas_constant*(TS+T0)))*  &
       2.0*gravity*TERR/(gas_constant*(TS+T0)**2)
   T09 = -slp9*psfc*EXP(2.0*gravity*TERR/(gas_constant*(TS+T0)))*  &
       2.0*gravity*TERR/(gas_constant*(TS+T0)**2)

   if ( t08 >= tc ) then
      if ( ts8 > tc ) then
         ts9=ts9-0.01*(ts-tc)*t09
      end if
      t09=0.0
   end if

   TS9     = TS9+T09
   TL9     = TS9*(psfc/PL)**XTERM
   psfc9   = psfc9+TS9*XTERM*(TL/PL)*(psfc/PL)**(XTERM-1)      
   PL9     = -TS9*XTERM*(TL*psfc/(PL*PL))*(psfc/PL)**(XTERM-1)
   THI9    = (1.0-COEF1/COEF2)*TL9
   TLO9    = COEF1/COEF2*TL9
   PL9     = PL9-(THI-TLO)/(COEF2*PL)*TL9
   p9(KHI) = p9(KHI)+((THI-TLO)/(COEF2*p(KHI))*(1.0-COEF1/COEF2))*TL9
   p9(KLO) = p9(KLO)+(THI-TLO)*COEF1/(COEF2*COEF2*p(KLO))*TL9
   t9(KHI) = t9(KHI)+THI9* (EPS+q(KHI))/(EPS*(1.0+q(KHI)))
   q9(KHI) = q9(KHI)+THI9*t(KHI)*(1.0-EPS)/(EPS*(1.0+q(KHI))**2)
   t9(KLO) = t9(KLO)+TLO9* (EPS+q(KLO))/(EPS*(1.0+q(KLO)))
   q9(KLO) = q9(KLO)+TLO9*t(KLO)*(1.0-EPS)/(EPS*(1.0+q(KLO))**2)

   psfc9=psfc9+PL9

   if (trace_use) call da_trace_exit("da_tpq_to_slp_adj")

end subroutine da_tpq_to_slp_adj


subroutine da_tv_profile(grid, info, n, pre_ma, tv_ma)

   !--------------------------------------------------------------------------
   ! Purpose: Calculates virtual temperature (tv_ma) on each level
   ! (pre_ma, pressure at the level) at the observed location (i,j). 
   ! dx, dxm, dy, dym are horizontal interpolation weighting.
   !--------------------------------------------------------------------------

   implicit none

   type (domain),         intent(in)  :: grid
   type (infa_type),      intent(in)  :: info
   integer,               intent(in)  :: n
   real,                  intent(out) :: pre_ma(kts-1:kte+1)
   real,                  intent(out) :: tv_ma(kts-1:kte+1)
                          
   integer :: ii,jj    ! index dimension.
   real    :: tv_m(2,2,kts:kte)     ! Virtual temperatures

   integer :: i, j      ! OBS location
   real    :: dx, dxm   ! interpolation weights.
   real    :: dy, dym   ! interpolation weights.

   if (trace_use_dull) call da_trace_entry("da_tv_profile")

   i   = info%i(1,n)
   j   = info%j(1,n)
   dx  = info%dx(1,n)
   dy  = info%dy(1,n)
   dxm = info%dxm(1,n)
   dym = info%dym(1,n)

   ! Virtual temperature

   do ii=i,i+1
      do jj=j,j+1
         tv_m(ii-i+1,jj-j+1,kts:kte) = grid%xb%t(ii,jj,kts:kte) * &
            (1.0 + 0.61*grid%xb%q(ii,jj,kts:kte))
      end do
   end do

   ! Horizontal interpolation to the obs. pt.

   pre_ma(kts:kte) = dym* ( dxm * grid%xb%p(i,j,kts:kte) + dx * grid%xb%p(i+1,j,kts:kte) ) + &
      dy * ( dxm * grid%xb%p(i,j+1,kts:kte) + dx * grid%xb%p(i+1,j+1,kts:kte) )

   tv_ma (kts:kte) = dym* ( dxm * tv_m (1,1,kts:kte) + dx * tv_m (2,1,kts:kte) ) + &
      dy * ( dxm * tv_m (1,2,kts:kte) + dx * tv_m (2,2,kts:kte) )

   if (trace_use_dull) call da_trace_exit("da_tv_profile")

end subroutine da_tv_profile


subroutine da_find_layer_tl(layer,tv,pre,pre_ma,tv_ma,ks,ke,TGL_tv,TGL_pre_ma,TGL_tv_ma)

   !-----------------------------------------------------------------------
   ! Purpose: tangent-linear routine for da_find_layer
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: ks, ke
   integer, intent(out)   :: layer
   real,    intent(in)    :: tv_ma(ks-1:ke+1)
   real,    intent(inout) :: pre_ma(ks-1:ke+1)
   real,    intent(in)    :: TGL_tv_ma(ks-1:ke+1)
   real,    intent(inout) :: TGL_pre_ma(ks-1:ke+1)
   real,    intent(in)    :: pre
   real,    intent(out)   :: tv
   real,    intent(out)   :: TGL_tv

   real    :: TGL_alpha
   integer :: k
   real    :: alpha, coef1, coef2

   if (trace_use_frequent) call da_trace_entry("da_find_layer_tl")

   ! coef1, coef2 are temporarily used in this routine

   if (pre >= pre_ma(ks)) then
      ! Below model bottom
      layer = ks
      coef1=log(pre/pre_ma(ks+1))/(pre_ma(ks)*     &
            (log(pre_ma(ks)/pre_ma(ks+1)))**2)
      coef2=log(pre_ma(ks)/pre)/(pre_ma(ks+1)*     &
            (log(pre_ma(ks)/pre_ma(ks+1)))**2)
      TGL_alpha = coef1 * TGL_pre_ma(ks) + coef2 * TGL_pre_ma(ks+1)
      alpha = log(pre_ma(ks)/pre)/log(pre_ma(ks)/pre_ma(ks+1))

      TGL_tv = (1.0-alpha)*TGL_tv_ma(ks) +               &
               (tv_ma(ks+1)-tv_ma(ks))*TGL_alpha +     &
               alpha*TGL_tv_ma(ks+1)
      TGL_pre_ma(ks-1) = 0.0
      tv = tv_ma(ks) * (1.0-alpha) + tv_ma(ks+1) * alpha
      pre_ma(ks-1) = pre
   else if (pre <= pre_ma(ke)) then
      ! Above model top
      layer = ke+1
      coef1=log(pre/pre_ma(ke))/(pre_ma(ke-1)*           &
            (log(pre_ma(ke-1)/pre_ma(ke)))**2)
      coef2=log(pre_ma(ke-1)/pre)/(pre_ma(ke)*           &
            (log(pre_ma(ke-1)/pre_ma(ke)))**2)
      TGL_alpha = coef1 * TGL_pre_ma(ke-1) + coef2 * TGL_pre_ma(ke)
      alpha = log(pre_ma(ke-1)/pre)/log(pre_ma(ke-1)/pre_ma(ke))

      TGL_tv = (1.0-alpha)*TGL_tv_ma(ke-1) +                 &
               (tv_ma(ke)-tv_ma(ke-1))*TGL_alpha +           &
               alpha*TGL_tv_ma(ke)
      TGL_pre_ma(ke+1) = 0.0
      tv = tv_ma(ke-1) * (1.0-alpha) + tv_ma(ke) * alpha
      pre_ma(ke+1) = pre
   else
      ! Between model layers
      do k=ks,ke-1
         if (pre>=pre_ma(k+1) .and. pre<pre_ma(k)) then
            layer = k+1
            coef1=log(pre/pre_ma(k+1))/(pre_ma(k)*   &
                  (log(pre_ma(k)/pre_ma(k+1)))**2)
            coef2=log(pre_ma(k)/pre)/(pre_ma(k+1)*   &
                  (log(pre_ma(k)/pre_ma(k+1)))**2)
            TGL_alpha = coef1 * TGL_pre_ma(k) + coef2 * TGL_pre_ma(k+1)
            alpha = log(pre_ma(k)/pre)/log(pre_ma(k)/pre_ma(k+1))
            TGL_tv = (1.0-alpha)*TGL_tv_ma(k) +                 &
                     (tv_ma(k+1)-tv_ma(k))*TGL_alpha +         &
                      alpha*TGL_tv_ma(k+1)
            tv = tv_ma(k) * (1.0-alpha) + tv_ma(k+1) * alpha
            exit
         end if
      end do
   end if

   if (trace_use_frequent) call da_trace_exit("da_find_layer_tl")
 
end subroutine da_find_layer_tl


subroutine da_thickness_adj(pre_ma,tv_ma,ks,ke,tv1,tv2,layer1,layer2,pre1,pre2,   &
   ADJ_pre_ma,ADJ_tv_ma,ADJ_tv1,ADJ_tv2,ADJ_thk)

   !-----------------------------------------------------------------------
   ! Purpose: adjoint routine for thickness
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: layer1,layer2
   integer, intent(in)    :: ks,ke
   real,    intent(in)    :: pre_ma(ks-1:ke+1)
   real,    intent(in)    :: tv_ma(ks-1:ke+1)
   real,    intent(inout) :: ADJ_pre_ma(ks-1:ke+1)
   real,    intent(inout) :: ADJ_tv_ma(ks-1:ke+1)
   real,    intent(in)    :: tv1,tv2
   real,    intent(inout) :: ADJ_tv1,ADJ_tv2
   real,    intent(in)    :: pre1,pre2
   real,    intent(inout) :: ADJ_thk

   integer :: k
   real    :: p_tmp(ks-1:ke+1)
   real    :: ADJ_p_tmp(ks-1:ke+1)

   if (trace_use_dull) call da_trace_entry("da_thickness_adj")

   ! p_tmp and ADJ_p_tmp are temporary (local) variables

   ADJ_p_tmp(:)=0.0

   p_tmp(layer1) = pre1
   p_tmp(layer2-1) = pre2
   do k=layer2,layer1-1
      p_tmp(k) = pre_ma(k)
   end do

   do k=layer2,layer1-1
      ADJ_p_tmp(k+1)  = ADJ_p_tmp(k+1) - 0.5*gas_constant/gravity *     &
                        ADJ_thk*tv_ma(k)/p_tmp(k+1)
      ADJ_p_tmp(k-1)  = ADJ_p_tmp(k-1) + 0.5*gas_constant/gravity *     &
                        ADJ_thk*tv_ma(k)/p_tmp(k-1)
      ADJ_tv_ma(k)    = ADJ_tv_ma(k)   + 0.5*gas_constant/gravity *     &
                        ADJ_thk*log(p_tmp(k-1)/p_tmp(k+1))
   end do

   do k=layer2,layer1-1
      ADJ_pre_ma(k) = ADJ_pre_ma(k) + ADJ_p_tmp(k)
   end do

   ADJ_thk = 0.5 * gas_constant/gravity * ADJ_thk
   ADJ_pre_ma(layer2) = ADJ_pre_ma(layer2) - ADJ_thk*tv2/pre_ma(layer2)
   ADJ_tv2 = ADJ_tv2 + ADJ_thk*log(pre2/pre_ma(layer2))
   ADJ_pre_ma(layer1-1) = ADJ_pre_ma(layer1-1) +              &
                          ADJ_thk*tv1/pre_ma(layer1-1)
   ADJ_tv1 = ADJ_tv1 + ADJ_thk*log(pre_ma(layer1-1)/pre1)

   if (trace_use_dull) call da_trace_exit("da_thickness_adj")

end subroutine da_thickness_adj


subroutine da_find_layer(layer,tv,pre,pre_ma,tv_ma,ks,ke)

   !-----------------------------------------------------------------------
   ! Purpose: find the vertical location in the Tv profile given
   ! a specific pressure and vertically interpolate Tv to that height.
   ! pre_ma,tv_ma give vertical profile of virtual temperature
   ! pre is a given pressure, alpha is the percentage of pre in the layer.
   ! layer,tv are calculated vertical layer and interpolated virtual temp.
   !-----------------------------------------------------------------------
 
   implicit none

   integer, intent(in)    :: ks, ke
   integer, intent(inout) :: layer
   real,    intent(inout) :: pre_ma(ks-1:ke+1)
   real,    intent(in)    :: tv_ma(ks-1:ke+1)
   real,    intent(in)    :: pre
   real,    intent(inout) :: tv

   integer :: k
   real    :: alpha

   if (trace_use_frequent) call da_trace_entry("da_find_layer")

   if (pre >= pre_ma(ks)) then
      ! Below model bottom
      layer = ks
      alpha = log(pre_ma(ks)/pre)/log(pre_ma(ks)/pre_ma(ks+1))
      tv = tv_ma(ks) * (1.0-alpha) + tv_ma(ks+1) * alpha
      pre_ma(ks-1)=pre
   else if (pre <= pre_ma(ke)) then
      ! Above model top
      layer = ke+1
      alpha = log(pre_ma(ke-1)/pre)/log(pre_ma(ke-1)/pre_ma(ke))
      tv = tv_ma(ke-1) * (1.0-alpha) + tv_ma(ke) * alpha
      pre_ma(ke+1) = pre
   else
      ! Between model layers 
      do k=ks,ke-1
         if (pre>=pre_ma(k+1) .and. pre<pre_ma(k)) then
            layer = k+1
            alpha = log(pre_ma(k)/pre)/log(pre_ma(k)/pre_ma(k+1))
            tv = tv_ma(k) * (1.0-alpha) + tv_ma(k+1) * alpha
            exit
         end if
      end do
   end if

   if (trace_use_frequent) call da_trace_exit("da_find_layer")
 
end subroutine da_find_layer


subroutine da_tv_profile_adj(grid,jo_grad_x,info, n,pre_ma,tv_ma,ADJ_pre_ma, ADJ_tv_ma)

   !-----------------------------------------------------------------------
   ! Purpose: adjoint routine for da_tv_profile
   !-----------------------------------------------------------------------

   implicit none

   type (x_type),         intent(inout) :: jo_grad_x ! grad_x(jo)
   type (infa_type),      intent(in)    :: info
   integer,               intent(in)    :: n
   type (domain),         intent(in)    :: grid
   real,                  intent(in)    :: pre_ma(kts-1:kte+1)
   real,                  intent(in)    :: tv_ma(kts-1:kte+1)
   real,                  intent(inout) :: ADJ_pre_ma(kts-1:kte+1)
   real,                  intent(inout) :: ADJ_tv_ma(kts-1:kte+1)

   integer :: ii,jj
   real    :: ADJ_tv_m(2,2,kts:kte)
   integer :: i, j      ! OBS location
   real    :: dx, dxm   ! interpolation weights.
   real    :: dy, dym   ! interpolation weights.

   if (trace_use_dull) call da_trace_entry("da_tv_profile_adj")

   i   = info%i(1,n)
   j   = info%j(1,n)
   dx  = info%dx(1,n)
   dy  = info%dy(1,n)
   dxm = info%dxm(1,n)
   dym = info%dym(1,n)

   ADJ_tv_m(1,1,kts:kte) = dym*dxm * ADJ_tv_ma (kts:kte)
   ADJ_tv_m(2,1,kts:kte) = dym*dx *  ADJ_tv_ma (kts:kte)
   ADJ_tv_m(1,2,kts:kte) = dy*dxm*   ADJ_tv_ma (kts:kte)
   ADJ_tv_m(2,2,kts:kte) = dy*dx*    ADJ_tv_ma (kts:kte)

   jo_grad_x%p(i,j,kts:kte)     = jo_grad_x%p(i,j,kts:kte) + dym*dxm  * ADJ_pre_ma(kts:kte)
   jo_grad_x%p(i+1,j,kts:kte)   = jo_grad_x%p(i+1,j,kts:kte) + dym*dx * ADJ_pre_ma(kts:kte)
   jo_grad_x%p(i,j+1,kts:kte)   = jo_grad_x%p(i,j+1,kts:kte) + dy*dxm * ADJ_pre_ma(kts:kte)
   jo_grad_x%p(i+1,j+1,kts:kte) = jo_grad_x%p(i+1,j+1,kts:kte) + dy*dx* ADJ_pre_ma(kts:kte)

   ADJ_tv_ma (kts:kte) = 0.0
   ADJ_pre_ma(kts:kte) = 0.0

   do ii=i,i+1
      do jj=j,j+1
         jo_grad_x%t(ii,jj,kts:kte) = jo_grad_x%t(ii,jj,kts:kte) + &
            ADJ_tv_m(ii-i+1,jj-j+1,kts:kte)*(1.0+0.61*grid%xb%q(ii,jj,kts:kte))
         jo_grad_x%q(ii,jj,kts:kte) = jo_grad_x%q(ii,jj,kts:kte) + &
            0.61*grid%xb%t(ii,jj,kts:kte)*ADJ_tv_m(ii-i+1,jj-j+1,kts:kte)
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_tv_profile_adj")

end subroutine da_tv_profile_adj


subroutine da_thickness(pre_ma,tv_ma,ks,ke,tv1,tv2,layer1,layer2,pre1,pre2,thk)

   !-----------------------------------------------------------------------
   ! Purpose: calculates the thickness between two layers 
   ! using vertical integration of virtual temperatures.
   ! pre1 and pre2 are two pressures for the two layers
   !-----------------------------------------------------------------------

   implicit none

   integer,intent(in)  :: layer1,layer2         ! two layers
   real,   intent(in)  :: tv1,tv2               ! virtual temp.
   real,   intent(in)  :: pre1,pre2             ! pressure

   integer,intent(in)  :: ks,ke
   real,   intent(in)  :: pre_ma(ks-1:ke+1)
   real,   intent(in)  :: tv_ma(ks-1:ke+1)
   real,   intent(out) :: thk                   ! thickness

   integer :: k
   real    :: p_tmp(ks-1:ke+1)

   if (trace_use_dull) call da_trace_entry("da_thickness")

   ! Thickness at the top and bottom parts of the layer.

   thk = 0.5 * gas_constant/gravity * (tv1*log(pre_ma(layer1-1)/pre1) +  &
      tv2*log(pre2/pre_ma(layer2)) )

   ! Temporary pressure

   p_tmp(layer1) = pre1
   p_tmp(layer2-1) = pre2
   do k = layer2, layer1-1
      p_tmp(k) = pre_ma(k)
   end do

   ! Vertical integration of the virtual temperature

   do k=layer2,layer1-1
      thk = thk + 0.5 * gas_constant/gravity * tv_ma(k) * log(p_tmp(k-1)/p_tmp(k+1))
   end do

   if (trace_use_dull) call da_trace_exit("da_thickness")

end subroutine da_thickness


subroutine da_find_layer_adj (layer,tv,pre,pre_ma,tv_ma,ks,ke,ADJ_tv,ADJ_pre_ma,ADJ_tv_ma)

   !-----------------------------------------------------------------------
   ! Purpose: adjoint routine for da_find_layer
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: ks, ke
   integer, intent(out)   :: layer
   real,    intent(in)    :: pre_ma(ks-1:ke+1)
   real,    intent(in)    :: tv_ma(ks-1:ke+1)
   real,    intent(inout) :: ADJ_pre_ma(ks-1:ke+1)
   real,    intent(inout) :: ADJ_tv_ma(ks-1:ke+1)
   real,    intent(in)    :: pre, tv
   real,    intent(inout) :: ADJ_tv

   integer :: k
   real    :: alpha, coef1, coef2
   real    :: ADJ_alpha

   if (trace_use_frequent) call da_trace_entry("da_find_layer_adj")

   if (pre >= pre_ma(ks)) then
      layer = ks
      coef1=log(pre/pre_ma(ks+1))/(pre_ma(ks)*     &
            (log(pre_ma(ks)/pre_ma(ks+1)))**2)
      coef2=log(pre_ma(ks)/pre)/(pre_ma(ks+1)*     &
            (log(pre_ma(ks)/pre_ma(ks+1)))**2)
      alpha = log(pre_ma(ks)/pre)/log(pre_ma(ks)/pre_ma(ks+1))

      ADJ_pre_ma(ks-1)= 0.0
      ADJ_tv_ma(ks)   = ADJ_tv_ma(ks) + (1.0-alpha)*ADJ_tv
      ADJ_alpha        = (tv_ma(ks+1)-tv_ma(ks))*ADJ_tv
      ADJ_tv_ma(ks+1) = ADJ_tv_ma(ks+1) + alpha*ADJ_tv

      ADJ_pre_ma(ks)    = ADJ_pre_ma(ks) + coef1 * ADJ_alpha
      ADJ_pre_ma(ks+1)  = ADJ_pre_ma(ks+1) + coef2 * ADJ_alpha
   else if (pre <= pre_ma(ke)) then
      layer = ke+1
      coef1=log(pre/pre_ma(ke))/(pre_ma(ke-1)*           &
            (log(pre_ma(ke-1)/pre_ma(ke)))**2)
      coef2=log(pre_ma(ke-1)/pre)/(pre_ma(ke)*           &
            (log(pre_ma(ke-1)/pre_ma(ke)))**2)
      alpha = log(pre_ma(ke-1)/pre)/log(pre_ma(ke-1)/pre_ma(ke))

      ADJ_pre_ma(ke+1)    = 0.0
      ADJ_tv_ma(ke-1)     = ADJ_tv_ma(ke-1) + (1.0-alpha)*ADJ_tv
      ADJ_alpha        = (tv_ma(ke)-tv_ma(ke-1))*ADJ_tv
      ADJ_tv_ma(ke)     = ADJ_tv_ma(ke) + alpha*ADJ_tv

      ADJ_pre_ma(ke-1) = ADJ_pre_ma(ke-1) + coef1 * ADJ_alpha
      ADJ_pre_ma(ke) = ADJ_pre_ma(ke) + coef2 * ADJ_alpha
   else
      do k=ks,ke-1
         if (pre>=pre_ma(k+1) .and. pre<pre_ma(k)) then
            layer = k+1
            coef1=log(pre/pre_ma(k+1))/(pre_ma(k)*   &
                  (log(pre_ma(k)/pre_ma(k+1)))**2)
            coef2=log(pre_ma(k)/pre)/(pre_ma(k+1)*   &
                  (log(pre_ma(k)/pre_ma(k+1)))**2)
            alpha = log(pre_ma(k)/pre)/log(pre_ma(k)/pre_ma(k+1))

            ADJ_tv_ma(k)     = ADJ_tv_ma(k) + (1.0-alpha)*ADJ_tv
            ADJ_alpha        = (tv_ma(k+1)-tv_ma(k))*ADJ_tv
            ADJ_tv_ma(k+1)   = ADJ_tv_ma(k+1) + alpha * ADJ_tv

            ADJ_pre_ma(k)   = ADJ_pre_ma(k) + coef1 * ADJ_alpha
            ADJ_pre_ma(k+1) = ADJ_pre_ma(k+1) + coef2 * ADJ_alpha
            exit
         end if
      end do
   end if

   ADJ_tv           = 0.0

   if (trace_use_frequent) call da_trace_exit("da_find_layer_adj")
 
end subroutine da_find_layer_adj


subroutine da_thickness_tl(pre_ma,tv_ma,ks,ke,tv1,tv2,layer1,layer2,pre1,pre2,   &
   TGL_pre_ma,TGL_tv_ma,TGL_tv1,TGL_tv2,TGL_thk)

   !-----------------------------------------------------------------------
   ! Purpose: tangent-linear routine for thickness
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)  :: layer1,layer2
   integer, intent(in)  :: ks,ke
   real,    intent(in)  :: pre_ma(ks-1:ke+1)
   real,    intent(in)  :: tv_ma(ks-1:ke+1)
   real,    intent(in)  :: TGL_pre_ma(ks-1:ke+1)
   real,    intent(in)  :: TGL_tv_ma(ks-1:ke+1)
   real,    intent(in)  :: tv1,tv2
   real,    intent(in)  :: TGL_tv1,TGL_tv2
   real,    intent(in)  :: pre1,pre2
   real,    intent(out) :: TGL_thk

   integer :: k
   real    :: p_tmp(ks-1:ke+1)
   real    :: TGL_p_tmp(ks-1:ke+1)

   if (trace_use_dull) call da_trace_entry("da_thickness_tl")

   TGL_thk = TGL_tv1*log(pre_ma(layer1-1)/pre1) +          &
             TGL_pre_ma(layer1-1)*tv1/pre_ma(layer1-1) +    &
             TGL_tv2*log(pre2/pre_ma(layer2)) -            &
             TGL_pre_ma(layer2)*tv2/pre_ma(layer2)
   TGL_thk = 0.5 * gas_constant/gravity * TGL_thk

   TGL_p_tmp(layer1) = 0.0
   p_tmp(layer1) = pre1
   TGL_p_tmp(layer2-1) = 0.0
   p_tmp(layer2-1) = pre2

   do k=layer2,layer1-1
      TGL_p_tmp(k) = TGL_pre_ma(k)
      p_tmp(k) = pre_ma(k)
   end do

   ! Vertical integration of the virtual temperature

   do k=layer2,layer1-1
      TGL_thk = TGL_thk + 0.5 * gas_constant/gravity * (   &
                TGL_tv_ma(k)*log(p_tmp(k-1)/p_tmp(k+1)) + &
                TGL_p_tmp(k-1)*tv_ma(k)/p_tmp(k-1) -       &
                TGL_p_tmp(k+1)*tv_ma(k)/p_tmp(k+1)     )
   end do

   if (trace_use_dull) call da_trace_exit("da_thickness_tl")

end subroutine da_thickness_tl


subroutine da_tv_profile_tl(grid,i,j,dx,dxm,dy,dym,pre_ma,tv_ma,TGL_pre_ma,TGL_tv_ma)

   !--------------------------------------------------------------------------
   ! Purpose: tangent-linear routine for da_tv_profile
   !--------------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)  :: grid
   integer,        intent(in)  :: i, j     ! OBS location
   real,           intent(in)  :: dx, dxm  ! interpolation weights.
   real,           intent(in)  :: dy, dym  ! interpolation weights.
   real,           intent(out) :: TGL_pre_ma(kts-1:kte+1)
   real,           intent(out) :: TGL_tv_ma(kts-1:kte+1)
   real,           intent(out) :: pre_ma(kts-1:kte+1)
   real,           intent(out) :: tv_ma(kts-1:kte+1)

   integer :: ii,jj
   real    :: tv_m(2,2,kts:kte)
   real    :: TGL_tv_m(2,2,kts:kte)

   if (trace_use_dull) call da_trace_entry("da_tv_profile_tl")

   do ii=i,i+1
      do jj=j,j+1
         TGL_tv_m(ii-i+1,jj-j+1,kts:kte) = grid%xa%t(ii,jj,kts:kte)*(1.0+0.61*grid%xb%q(ii,jj,kts:kte)) +  &
                                      0.61*grid%xb%t(ii,jj,kts:kte)*grid%xa%q(ii,jj,kts:kte)
         tv_m(ii-i+1,jj-j+1,kts:kte)     = grid%xb%t(ii,jj,kts:kte)*(1.0+0.61*grid%xb%q(ii,jj,kts:kte))
      end do
   end do
 
   TGL_pre_ma(kts:kte) = dym* ( dxm * grid%xa%p(i,j,kts:kte) + dx * grid%xa%p(i+1,j,kts:kte) ) + &
                       dy * ( dxm * grid%xa%p(i,j+1,kts:kte) + dx * grid%xa%p(i+1,j+1,kts:kte) )
   TGL_tv_ma (kts:kte) = dym* ( dxm * TGL_tv_m(1,1,kts:kte) + dx * TGL_tv_m(2,1,kts:kte) ) + &
                       dy * ( dxm * TGL_tv_m(1,2,kts:kte) + dx * TGL_tv_m(2,2,kts:kte) )
   pre_ma(kts:kte) = dym* ( dxm * grid%xb%p(i,j,kts:kte) + dx * grid%xb%p(i+1,j,kts:kte) ) + &
                   dy * ( dxm * grid%xb%p(i,j+1,kts:kte) + dx * grid%xb%p(i+1,j+1,kts:kte) )
   tv_ma (kts:kte) = dym* ( dxm * tv_m (1,1,kts:kte) + dx * tv_m (2,1,kts:kte) ) + &
                   dy * ( dxm * tv_m (1,2,kts:kte) + dx * tv_m (2,2,kts:kte) )

   if (trace_use_dull) call da_trace_exit("da_tv_profile_tl")
 
end subroutine da_tv_profile_tl



subroutine da_transform_xtoztd(grid)

!------------------------------------------------------------------------
!  Purpose: to compute the Zenith Total Delay, and save it to xb%ztd.
!
!  Both of the wet and dry delay are computed based on Vedel and Huang,
!          J. Meteor. Soc., 82, 459-472, 2004.
!
!       ** Equation (3) in Vedel and Huang is wrong.
! 
!                 ported by Yong-Run Guo  05/12/2008 from wrf3dvar.
!------------------------------------------------------------------------

   implicit none
   
   type (domain), intent(inout) :: grid

   integer :: i, j, k

   real    :: const, part, term1, term2, wzd, hzd, zf

   if (trace_use) call da_trace_entry("da_transform_xtoztd")

!--WEIGHTED SUM OF VERTICAL COLUMN
   do j=jts, jte
      do i=its, ite

! Wet delay:
      wzd = 0.0
      do k=kts, kte
        const  = (grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k)) / a_ew
        part   = grid%xb%p(i,j,k)*grid%xb%q(i,j,k) / grid%xb%t(i,j,k)
        term1  = part * const * wdk1
        term2  = part * const * wdk2 / grid%xb%t(i,j,k)
        wzd    = wzd + term1 + term2
      enddo

! Hydrostatic delay (Saastamoinen 1972):
       zf = (1.0 - zdk2*cos(2.0*grid%xb%lat(i,j)*radian) - zdk3*grid%xb%terr(i,j))
      hzd = zdk1 * grid%xb%psfc(i,j) / zf

!-----To save the ZTD in cm to ztd:
      grid%xb%ztd(i,j) = (wzd + hzd) * 1.e2
    enddo
   enddo

end subroutine da_transform_xtoztd

SUBROUTINE DA_Transform_XToZTD_Lin( grid )

!------------------------------------------------------------------------
!  Purpose: to compute the Zenith Total Delay, and save it to grid%xa%ztd.
!
!                                Yong-Run Guo  05/20/2008
!------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   integer :: i, j, K

   real    :: const, part, parta, term1, term2, wzd, hzd, zf

!--WEIGHTED SUM OF VERTICAL COLUMN
   do j=jts, jte
   do i=its, ite

! Wet delay:
      wzd = 0.0
      do k=kts, kte
        const  = (grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k)) / a_ew
        part   = grid%xb%p(i,j,k)*grid%xb%q(i,j,k) / grid%xb%t(i,j,k)

        parta  = (grid%xb%q(i,j,k)*grid%xa%p(i,j,k) + grid%xb%p(i,j,k)*grid%xa%q(i,j,k) &
                  - grid%xb%p(i,j,k)*grid%xb%q(i,j,k)*grid%xa%t(i,j,k) / grid%xb%t(i,j,k)) &
                 / grid%xb%t(i,j,k)  
        term1 = parta * const * wdk1
        term2 = (parta * const * wdk2                                &
                 - part * const * wdk2 * grid%xa%t(i,j,k) / grid%xb%t(i,j,k))  &
                / grid%xb%t(i,j,k)
        wzd   = wzd + term1 + term2
      enddo

! Hydrostatic delay (Saastamoinen 1972):
       zf = (1.0 - zdk2*cos(2.0*grid%xb%lat(i,j)*radian) - zdk3*grid%xb%terr(i,j))
      hzd = zdk1 * grid%xa%psfc(i,j) / zf

!-----To save the ZTD in cm to ztd:
      grid%xa%ztd(i,j) = (wzd + hzd) * 1.e2
    enddo
   enddo
 
END SUBROUTINE DA_Transform_XToZTD_Lin

SUBROUTINE DA_Transform_XToZTD_Adj( grid )

!------------------------------------------------------------------------
!  Purpose: to compute the adjoint of the Zenith Total Delay, 
!           and save it to grid%xa.
!
!                                Yong-Run Guo  05/20/2008
!------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   integer :: i, j, K

   real    :: const, part, parta, term1, term2, wzd, hzd, zf

!--WEIGHTED SUM OF VERTICAL COLUMN
   do j=jts, jte
   do i=its, ite
     wzd = grid%xa%ztd(i,j) * 1.e2
     hzd = grid%xa%ztd(i,j) * 1.e2
!
      zf = (1.0 - zdk2*cos(2.0*grid%xb%lat(i,j)*radian) - zdk3*grid%xb%terr(i,j))
      grid%xa%psfc(i,j) = grid%xa%psfc(i,j) + zdk1 * hzd / zf
!
      do k=kts, kte
        const  = (grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k)) / a_ew
        part   = grid%xb%p(i,j,k)*grid%xb%q(i,j,k) / grid%xb%t(i,j,k)

        term1  = wzd
        term2  = wzd

        parta  = term2 * const * wdk2 / grid%xb%t(i,j,k)
        grid%xa%t(i,j,k) = grid%xa%t(i,j,k) - part * const * wdk2 * term2 &
                                    / (grid%xb%t(i,j,k)*grid%xb%t(i,j,k))

        parta  = parta + term1 * const * wdk1

        grid%xa%p(i,j,k) = grid%xa%p(i,j,k) + grid%xb%q(i,j,k)*parta/grid%xb%t(i,j,k)
        grid%xa%q(i,j,k) = grid%xa%q(i,j,k) + grid%xb%p(i,j,k)*parta/grid%xb%t(i,j,k)
        grid%xa%t(i,j,k) = grid%xa%t(i,j,k) - grid%xb%p(i,j,k)*grid%xb%q(i,j,k)*parta &
                                    / (grid%xb%t(i,j,k)*grid%xb%t(i,j,k))
      enddo

    enddo
   enddo
 
END SUBROUTINE DA_Transform_XToZTD_Adj

subroutine da_transform_xtotpw(grid)

   !---------------------------------------------------------------------
   ! Purpose: weighted sum of vertical column
   !---------------------------------------------------------------------

   implicit none
   
   type (domain),  intent(inout) :: grid

   integer :: i, j, k

   real    :: pw

   if (trace_use) call da_trace_entry("da_transform_xtotpw")

   do j=jts, jte
      do i=its, ite
         pw = 0.0
         do k=kts, kte
            pw = pw + (grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k)) &
                    * (grid%xa%q(i,j,k)*grid%xb%rho(i,j,k) &
                    +  grid%xb%q(i,j,k)*grid%xa%rho(i,j,k))
         end do
 
         ! To convert the unit of PW to cm:
         grid%xa%tpw(i,j) = 0.1*pw
      end do
   end do

   if (trace_use) call da_trace_exit("da_transform_xtotpw")
 
end subroutine da_transform_xtotpw


subroutine da_transform_xtotpw_adj(grid)

   !---------------------------------------------------------------------
   ! Purpose: TBD
   !---------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid

   integer :: i, j, k, is, js, ie, je

   real    :: pw, dzpw

   if (trace_use) call da_trace_entry("da_transform_xtotpw_adj")

   is = its
   js = jts

   ie = ite
   je = jte

   if (test_transforms) then
      is = its-1
      js = jts-1

      ie = ite+1
      je = jte+1

      if (is < ids) is = ids
      if (js < jds) js = jds

      if (ie > ide) ie = ide
      if (je > jde) je = jde
   end if

   do j=js, je
      do i=is, ie
         pw = 0.1*grid%xa%tpw(i,j)

         do k=kts, kte
            dzpw = (grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k))*pw

            grid%xa%  q(i,j,k)=grid%xa%  q(i,j,k)+dzpw*grid%xb%rho(i,j,k)
            grid%xa%rho(i,j,k)=grid%xa%rho(i,j,k)+dzpw*grid%xb%  q(i,j,k)
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_transform_xtotpw_adj")
 
end subroutine da_transform_xtotpw_adj


subroutine da_transform_xtogpsref(grid)

   !-------------------------------------------------------------------
   ! Purpose: TBD
   !-------------------------------------------------------------------

   implicit none

   ! input : grid%xb%q, grid%xb%p, grid%xb%t, and xp
   ! output: grid%xb%ref

   type (domain), intent(inout) :: grid
   
   integer :: i, j, k, ij
   real    :: partone, parttwo, dividnd

   if (trace_use_dull) call da_trace_entry("da_transform_xtogpsref")

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, i, j, k, partone, dividnd, parttwo)
   do ij = 1 , grid%num_tiles

   do k = kts, kte
      do j = grid%j_start(ij), grid%j_end(ij)
         do i = its, ite
            ! calculate refractivity

            !  1, Hydrostatic part of refractivity:
            !     Note: p in Pascal.

            partone  = 0.776*grid%xb%p(i,j,k)/grid%xb%t(i,j,k)

            !  2, (Wet part) / (hydrostatic part):
            !     Note: grid%xb%q its the specific humidity --- an analysis variable

            dividnd  = grid%xb%t(i,j,k)*(0.622+0.378*grid%xb%q(i,j,k))
            parttwo  = 1.0+coeff*grid%xb%q(i,j,k)/dividnd

            !  3, Refractivity:
            grid%xb%ref(i,j,k)= partone * parttwo
            grid%xb%reflog(i,j,k)= log(grid%xb%ref(i,j,k))
         end do
      end do
   end do

   end do
   !$OMP END PARALLEL DO

   if (trace_use_dull) call da_trace_exit("da_transform_xtogpsref")
   
end subroutine da_transform_xtogpsref


subroutine da_transform_xtogpsref_adj(grid)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! input: grid%xb%q, grid%xb%p, grid%xb%t, grid%xa%ref, and xp
   ! output: grid%xa%q, grid%xa%p, grid%xa%t
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid

   integer :: i, j, k, is, ie, js, je, ks, ke
   real    :: partone, parttwo, dividnd
   real    :: partone9,parttwo9,dividnd9

   if (trace_use_dull) call da_trace_entry("da_transform_xtogpsref_adj")

   ! 1.0 Get the range for a tile

   is = its; ie = ite
   js = jts; je = jte
   ks = kts; ke = kte      

   if (test_transforms) then
      is = its-1
      js = jts-1

      ie = ite+1
      je = jte+1

      if (is < ids) is = ids
      if (js < jds) js = jds

      if (ie > ide) ie = ide
      if (je > jde) je = jde
   end if

   !------------------------------------------------------------------------------
   !  [2.0] Calculate the adjoint for GPS refractivity:
   !------------------------------------------------------------------------------

   do k = ks, ke
      do j = js, je
         do i = is, ie
            !  Note: p in Pascal, q is the specific humidity

            !   2.1 basic state

            partone9 = 0.776*grid%xb%p(i,j,k)/grid%xb%t(i,j,k)
            dividnd9 = grid%xb%t(i,j,k)*(0.622+0.378*grid%xb%q(i,j,k))
            parttwo9 = 1.0+coeff*grid%xb%q(i,j,k)/dividnd9

            !   2.2 Adjoint of partone and parttwo

            partone = grid%xa%ref(i,j,k) * parttwo9
            parttwo = grid%xa%ref(i,j,k) * partone9

            !   2.3 Adjoint of q and dividnd

            grid%xa%q(i,j,k) = grid%xa%q(i,j,k) + coeff*parttwo/dividnd9
            dividnd=-coeff*grid%xb%q(i,j,k)*parttwo/(dividnd9*dividnd9)

            !   2.4 Adjoint of t, q, and p

            grid%xa%t(i,j,k) = grid%xa%t(i,j,k) + dividnd*(0.622+0.378*grid%xb%q(i,j,k))
            grid%xa%q(i,j,k) = grid%xa%q(i,j,k) + grid%xb%t(i,j,k)*0.378*dividnd

! Zero out the REF impact on xa%p (YRG, 10/08/2009):
!            grid%xa%p(i,j,k) = grid%xa%p(i,j,k) + 0.776*partone/grid%xb%t(i,j,k)
            
            grid%xa%t(i,j,k) = grid%xa%t(i,j,k) - &
               0.776*grid%xb%p(i,j,k)*partone/(grid%xb%t(i,j,k)*grid%xb%t(i,j,k))

            !  [3.0] Clean the grid%xa%ref    

            grid%xa%ref(i,j,k) = 0.0
         end do
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_transform_xtogpsref_adj")

end subroutine da_transform_xtogpsref_adj


subroutine da_transform_xtogpsref_lin(grid)

   !---------------------------------------------------------------------
   ! Purpose: TBD                            
   !---------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid
   
   integer :: i, j, k
   real    :: partone, parttwo, dividnd
   real    :: partone9,parttwo9,dividnd9

   if (trace_use_dull) call da_trace_entry("da_transform_xtogpsref_lin")

   ! 2.0 Loop for a tile

   do k = kts, kte
      do j = jts, jte
         do i = its, ite
            !-----calculate refractivity
            !     Note: p in Pascal, q its the specific humidity

            ! 2.1 Part one: hydrostatic part

! Do not consider the pressure perturbation (YRG, 10/08/2009):
!            partone  = 0.776*grid%xa%p(i,j,k)/grid%xb%t(i,j,k) &

            partone  = &
               - 0.776*grid%xb%p(i,j,k)*grid%xa%t (i,j,k)/ &
               (grid%xb%t(i,j,k)*grid%xb%t(i,j,k))

            partone9 = 0.776*grid%xb%p(i,j,k)/grid%xb%t(i,j,k)

            ! 2.2 Part two considering the moitsture:

            dividnd  = grid%xa%t (i,j,k)*(0.622+0.378*grid%xb%q(i,j,k)) &
                     + grid%xb%t(i,j,k)*       0.378*grid%xa%q(i,j,k)
            dividnd9 = grid%xb%t(i,j,k)*(0.622+0.378*grid%xb%q(i,j,k))
            parttwo  =     coeff*grid%xa%q(i,j,k)/dividnd9 &
                     -     coeff*grid%xb%q(i,j,k)*dividnd /(dividnd9*dividnd9)
            parttwo9 = 1.0+coeff*grid%xb%q(i,j,k)/dividnd9

            ! 2.3 Total refractivity

            grid%xa%ref(i,j,k)= partone9 * parttwo + partone  * parttwo9
         end do
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_transform_xtogpsref_lin")
   
end subroutine da_transform_xtogpsref_lin


subroutine da_check_rh(grid)

   !---------------------------------------------------------------------------
   ! Purpose: Remove RH over 100% or under 10%
   !          Made the modification to those levels, which RH are less than 95%
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   integer   :: imod(kts:kte)
   real      :: rhtol(kts:kte)
   real      :: x_qs(kts:kte)
   real      :: dz(kts:kte)

   integer :: i, j, k, ij
   real    :: tol_adjust_moist, tol_moist, oldrha, each_moist, es, weight
   real    :: upper_modify_rh, lower_modify_rh
   real    :: t_tropp(ims:ime,jms:jme)
   integer :: k_tropp(ims:ime,jms:jme)

   if (trace_use) call da_trace_entry("da_check_rh")


   upper_modify_rh = 95.0
   lower_modify_rh = 11.0

   ! To get the grid%xa%rh for the moisture adjustments
   call da_tpq_to_rh_lin (grid)

   ! find the k index of tropopause
   t_tropp = 999.0  ! initialize
   k_tropp = kte    ! initialize
   do k = kte, kts, -1
      do j = jts, jte
         do i = its, ite
            if ( grid%xb%t(i,j,k) < t_tropp(i,j) .and.  &
                 grid%xb%p(i,j,k) > 5000.0 ) then  ! tropopause should not
                                                   ! be higher than 50 hPa
               t_tropp(i,j) = grid%xb%t(i,j,k)
               k_tropp(i,j) = k
            end if
         end do
      end do
   end do

   !$OMP PARALLEL DO  SCHEDULE (DYNAMIC, 1) &
   !$OMP PRIVATE ( i, j, k, tol_adjust_moist, tol_moist) &
   !$OMP PRIVATE ( weight, oldrha, each_moist, imod, dz, x_qs, rhtol, es)
!  do ij = 1 , grid%num_tiles

   do i=its,ite
      !do j=grid%j_start(ij), grid%j_end(ij)
      do j=jts, jte

         tol_adjust_moist = 0.0
         tol_moist        = 0.0

         do k=kts,kte
            dz(k)=grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k)

            imod(k)           = 0
            x_qs(k)           = 0.0
            rhtol(k)          = grid%xb%rh(i,j,k)+grid%xa%rh(i,j,k)
         end do

         do k=kts, k_tropp(i,j)
            if (rhtol(k) .gt. maximum_rh) then
               oldrha       = grid%xa%rh(i,j,k)
               ! modify grid%xa%rh
               grid%xa%rh(i,j,k) = maximum_rh - grid%xb%rh(i,j,k)

               call da_tp_to_qs(grid%xb%t(i,j,k)+grid%xa%t(i,j,k), &
                  grid%xb%p(i,j,k)+grid%xa%p(i,j,k), es, x_qs(k))

               ! calculate grid%xa%q
               call da_tprh_to_q_lin1(grid%xb%t(i,j,k), grid%xb%p(i,j,k), &
                  grid%xb%es(i,j,k), grid%xb%q(i,j,k), grid%xb%rh(i,j,k),  grid%xa%t(i,j,k), &
                  grid%xa%p(i,j,k), grid%xa%rh(i,j,k), grid%xa%q(i,j,k))

               tol_adjust_moist = tol_adjust_moist + x_qs(k)*(oldrha - &
                  grid%xa%rh(i,j,k))* dz(k)*(grid%xb%rho(i,j,k)+grid%xa%rho(i,j,k))
               imod(k)=-1
            end if

            if (rhtol(k) .lt. minimum_rh) then
               oldrha           = grid%xa%rh(i,j,k)
               grid%xa%rh(i,j,k)     = minimum_rh - grid%xb%rh(i,j,k)
               call da_tp_to_qs(grid%xb%t(i,j,k)+grid%xa%t(i,j,k), &
                  grid%xb%p(i,j,k)+grid%xa%p(i,j,k), es, x_qs(k))

               call da_tprh_to_q_lin1(grid%xb%t(i,j,k), grid%xb%p(i,j,k), &
                  grid%xb%es(i,j,k), grid%xb%q(i,j,k), grid%xb%rh(i,j,k),  grid%xa%t(i,j,k), &
                  grid%xa%p(i,j,k), grid%xa%rh(i,j,k), grid%xa%q(i,j,k))


               tol_adjust_moist = tol_adjust_moist + x_qs(k)*(oldrha - &
                  grid%xa%rh(i,j,k))* dz(k)*(grid%xb%rho(i,j,k)+grid%xa%rho(i,j,k))
               imod(k)=-1
            end if
         end do

         if (tol_adjust_moist .gt. 0.0) then
            do k=kts, k_tropp(i,j)
               if (rhtol(k) .lt. upper_modify_rh .and. imod(k) .eq. 0) then
                  call da_tp_to_qs(grid%xb%t(i,j,k)+grid%xa%t(i,j,k), &
                                    grid%xb%p(i,j,k)+grid%xa%p(i,j,k),es,x_qs(k))

                  each_moist   = rhtol(k)*x_qs(k)* &
                                 dz(k)*(grid%xb%rho(i,j,k)+grid%xa%rho(i,j,k))
                  tol_moist    = tol_moist + each_moist
                  imod(k)      = 1
               end if
            end do
         end if

         if (tol_adjust_moist .lt. 0.0) then
            do k=kts, k_tropp(i,j)
               if (rhtol(k) .gt. lower_modify_rh .and. imod(k) .eq. 0) then
                  call da_tp_to_qs(grid%xb%t(i,j,k)+grid%xa%t(i,j,k), &
                                    grid%xb%p(i,j,k)+grid%xa%p(i,j,k), es, x_qs(k))

                  each_moist   = rhtol(k)*x_qs(k)* &
                                 dz(k)*(grid%xb%rho(i,j,k)+grid%xa%rho(i,j,k))
                  tol_moist    = tol_moist + each_moist
                  imod(k)      = 1
               end if
            end do
         end if

         if (tol_moist > 0) then
           weight       = tol_adjust_moist/tol_moist
           do k=kts, k_tropp(i,j)
             if (imod(k) .eq. 1) then
               grid%xa%rh(i,j,k) = grid%xa%rh(i,j,k)+(grid%xb%rh(i,j,k)+grid%xa%rh(i,j,k))*weight

! To guarantee the adjusted relative humidity will not be out of the range (YRG, 10/21/2008):
               oldrha = grid%xa%rh(i,j,k)+grid%xb%rh(i,j,k)
               if ( (oldrha > maximum_rh) ) then
                  grid%xa%rh(i,j,k) = maximum_rh - grid%xb%rh(i,j,k)
                  if (print_detail_xa ) &
                     write(unit=stdout, fmt='(3I5," Warning=== Adjusted RH > maximum_rh=",f10.2,&
                            & " total_rh, xb%rh, xa%rh:",3f10.2)') &
                          i, j, k, maximum_rh, oldrha, grid%xb%rh(i,j,k), grid%xa%rh(i,j,k)
               else if ( oldrha < minimum_rh ) then
                  grid%xa%rh(i,j,k) = minimum_rh - grid%xb%rh(i,j,k)
                  if (print_detail_xa ) &
                     write(unit=stdout, fmt='(3I5," Warning=== Adjusted RH < minimum_rh=",f10.2,&
                            & " total_rh, xb%rh, xa%rh:",3f10.2)') &
                          i, j, k, minimum_rh, oldrha, grid%xb%rh(i,j,k), grid%xa%rh(i,j,k)
               endif
! ...........................................................................................

               call da_tprh_to_q_lin1(grid%xb%t(i,j,k), grid%xb%p(i,j,k), grid%xb%es(i,j,k), &
                                      grid%xb%q(i,j,k), grid%xb%rh(i,j,k),  grid%xa%t(i,j,k), &
                                      grid%xa%p(i,j,k), grid%xa%rh(i,j,k), grid%xa%q(i,j,k))

             end if
           end do
         end if
      end do
   end do

!  end do
   !$OMP END PARALLEL DO

   if (trace_use) call da_trace_exit("da_check_rh")

end subroutine da_check_rh


subroutine da_check_rh_simple (grid)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout)  :: grid

   integer :: i, j, k ! Loop counters.
   real    :: q_new   ! Modified rh.
   real    :: es      ! Dummy output.
   real    :: t_new, p_new
   real    :: ramax,ramin,dqmax,dqmin    
   real    :: t_tropp(ims:ime,jms:jme)
   integer :: k_tropp(ims:ime,jms:jme)

   if (trace_use_dull) call da_trace_entry("da_check_rh_simple")        

   ! ramax=maximum_rh/100.0
   ! ramin=minimum_rh/100.0
   ramax = 1.0
   ramin = 1.0e-6

   ! find the k index of tropopause
   t_tropp = 999.0  ! initialize
   k_tropp = kte    ! initialize
   do k = kte, kts, -1
      do j = jts, jte
         do i = its, ite
            if ( grid%xb%t(i,j,k) < t_tropp(i,j) .and.  &
                 grid%xb%p(i,j,k) > 5000.0 ) then  ! tropopause should not
                                                   ! be higher than 50 hPa
               t_tropp(i,j) = grid%xb%t(i,j,k)
               k_tropp(i,j) = k
            end if
         end do
      end do
   end do

   do k = kts, kte
      do j = jts, jte
         do i = its, ite
            p_new  = grid%xb % p(i,j,k) + grid%xa %  p(i,j,k)
            t_new  = grid%xb % t(i,j,k) + grid%xa %  t(i,j,k)

            if ( k > k_tropp(i,j) ) then  ! limit q incement above tropopause
               grid%xa % q(i,j,k) = 0.0
            else
               call da_tp_to_qs(t_new, p_new, es, q_new)
               dqmax=q_new*ramax - grid%xb % q(i,j,k)        
               dqmin=q_new*ramin - grid%xb % q(i,j,k)        
               grid%xa % q(i,j,k)=min(max(grid%xa % q(i,j,k),dqmin),dqmax)
            end if
         end do
      end do
   end do  

   if (trace_use_dull) call da_trace_exit("da_check_rh_simple") 

end subroutine da_check_rh_simple


subroutine da_get_q_error( p, t, q, t_error, rh_error, q_error)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)  :: p        ! Observed pressure.
   real, intent(in)  :: t        ! Observed temperature.
   real, intent(in)  :: q        ! Observed specific humidity.
   real, intent(in)  :: t_error  ! Temperature observation error.
   real, intent(in)  :: rh_error ! RH observation error.
   real, intent(out) :: q_error  ! q observation error.

   real :: es    ! Saturation vapor pressure.
   real :: qs    ! Saturation specific humidity.
   real :: rh    ! Relative humidity.

   real :: f_qv_from_rh
   external f_qv_from_rh

   if (trace_use_frequent) call da_trace_entry("da_get_q_error")

   if (ABS(p - missing_r) > 1.0 .AND. ABS(t - missing_r) > 1.0 .AND. ABS(q - missing_r) > 1) then
      ! Calculate rh:
      call da_tp_to_qs( t, p, es, qs)

      rh = 100.0 * q / qs
      if (rh > 100.0) rh = 100.0

      ! Get observation error for qv. Is this right?
      q_error = f_qv_from_rh( rh_error, t_error, rh, t, p)
   else
      q_error = missing_r
   end if

   if (q_error == 0.0) q_error = missing_r

   if (trace_use_frequent) call da_trace_exit("da_get_q_error")

end subroutine da_get_q_error


subroutine da_roughness_from_lanu(ltbl, mminlu, date, lanu, rough)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer             ,   intent(in)    :: ltbl
   character (len=*)   ,   intent(in)    :: mminlu
   character (len=19)  ,   intent(in)    :: date
   real, dimension(ims:ime,jms:jme),   intent(in)    :: lanu 
   real, dimension(ims:ime,jms:jme),   intent(out)   :: rough 

   integer                               :: LS, LC, LI, LUCATS, LuseAS, &
                                           LUMATCH, year, month, day,  &
                                           julday, Isn, io_error, &
                                           m1, m2, n1, n2 
   real                                  :: albd, slmo, sfem
   real(kind=4), dimension(50,2)         :: sfz0
   character (len=256)                   :: LUtype
   logical                               :: exist

   if (trace_use) call da_trace_entry("da_roughness_from_lanu")

   read(unit=date,fmt='(I4,1x,I2,1X,I2)') year, month, day
   call da_julian_day(year,month,day,Julday, 1)
   Isn = 1
   if (JULDAY < 105 .OR. JULDAY > 288) Isn=2

   inquire (file = 'LANDUSE.TBL', exist = exist)

   if (exist) then
      open (unit = ltbl, file = 'LANDUSE.TBL', form='formatted', &
                     action = 'read', iostat = io_error)
   else
      call da_error("da_roughness_from_lanu.inc",37,&
         (/"Cannot open file LANDUSE.TBL for conversion of roughness"/))
   end if

   lumatch=0  

   do
      read (unit=ltbl,fmt='(A)', iostat=io_error) lutype
      if (io_error /= 0) exit
      read (unit=ltbl,fmt=*, iostat=io_error) lucats,luseas

      if (trim(lutype) == trim(mminlu)) lumatch=1 

      do LS=1,LuseAS 
         read (unit=ltbl,fmt=*)  
         do lc=1,lucats 
            if (trim(lutype) == trim(mminlu)) then 
               read (unit=ltbl,fmt=*) li, albd, slmo, sfem, sfz0(lc,ls)
               ! prevent compiler whinge
               if (albd == 0.0 .or. sfem == 0.0 .or. slmo == 0.0) then
               end if
               if (LC /= LI) then
                 call da_error("da_roughness_from_lanu.inc",59, &
                   (/"Missing landuse: lc"/))
               end if
            else 
               read (unit=ltbl,fmt=*) 
            end if 
         end do 
      end do
   end do

   close (unit=ltbl)

   if (lumatch == 0)then
    call da_error("da_roughness_from_lanu.inc",72,&
      (/"landuse in input file does not match lutable"/))
   end if   

   m1 = its
   m2 = ite
   n1 = jts
   n2 = jte

   do lc = m1,m2
      do ls = n1,n2
         Li = int(lanu(lc,ls)+0.001)
         rough(lc,ls) =  sfz0(Li,Isn)/100.0
      end do
   end do

   if (trace_use) call da_trace_exit("da_roughness_from_lanu")

end subroutine da_roughness_from_lanu


subroutine da_julian_day(NY,NM,ND,JD,METHOD) 

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! method = 1: input ---- ny, nm, nd.  output ---- jd                         
   ! method = 2: input ---- ny, jd.      output ---- nm, nd            
   !----------------------------------------------------------------------- 
    
   implicit none

   integer,  intent(in)    :: METHOD, NY
   integer,  intent(inout) :: NM, ND, JD

   integer                 :: LOOP
   integer, dimension(12)  :: MDAY = (/31,28,31,30,31,30,31,31,30,31,30,31/)

   if (trace_use_dull) call da_trace_entry("da_julian_day")

   if (mod(ny,4) == 0) then
      mday(2)=29      
      if (mod(ny,100) == 0) then
         mday(2)=28
         if (mod(ny,400) == 0) then
            mday(2)=29
         end if
      end if
   end if

   if (method == 1) then                  
      jd=nd
      juday: do loop=1,nm-1                  
         jd=jd+mday(loop)                
      end do juday                           
   else if (method == 2) then             
      nm=1                               
      nd=jd
      do loop=1,11                    
         if (nd <= mday(loop)) exit

         nd=nd-mday(loop)     
         nm=nm+1                      
      end do
   end if   

   if (trace_use_dull) call da_trace_exit("da_julian_day")                             

end subroutine da_julian_day


subroutine da_sfc_wtq (psfc, tg, ps, ts, qs, us, vs, &
   hs, roughness, xland, dx, u10, v10, t2, q2, regime, &
   has_lsm, regime_wrf, qsfc_wrf, znt_wrf, ust_wrf, mol_wrf, hfx, qfx, pblh)

   !---------------------------------------------------------------------------
   ! Purpose: Calculate the  10m wind, 2m temperature and moisture based on the
   ! similarity theory/
   !
   !  The unit for pressure   : psfc, ps          is Pa.
   !  The unit for temperature: tg, ts, t2        is K.
   !  The unit for moisture   : qs, q2            is kg/kg.
   !  The unit for wind       : us, vs, u10, v10  is m/s.
   !  The unit for height     : hs, roughness     is m.
   !  xland and regime are dimensionless.
   !
   ! History: Nov 2010 - improve calculation consistency with WRF model (Eric Chiang)
   !          Jul 2015 - further improvement on consistency
   !
   ! Reference:
   ! ---------
   !
   !  input Variables:
   ! 
   !   psfc, tg               : surface pressure and ground temperature
   !   ps, ts, qs, us, vs, hs : model variable at lowlest half sigma level
   !   dx  (m)                : horizontal resolution
   !
   !
   !  Constants:
   !
   !   hs                     : height at the lowest half sigma level
   !   roughness              : roughness
   !   xland                  : land-water-mask (=2 water, =1 land)
   !
   !  output Variables:
   !
   !   regime                 : PBL regime
   !   u10, v10               : 10-m high observed wind components
   !   t2 , q2                : 2-m high observed temperature and mixing ratio
   !
   !---------------------------------------------------------------------------
   !  
   !                      psim  : mechanical psi at lowlest sigma level
   !                      psim2 : mechanical psi at 2m 
   !                      psimz : mechanical psi at 10m 
   !
   !---------------------------------------------------------------------------

   implicit none

   real, intent (in)  :: ps , ts , qs , us, vs
   real, intent (in)  :: psfc, tg
   real, intent (in)  :: hs, roughness, xland
   real, intent (out) :: regime
   real, intent (out) :: u10, v10, t2, q2
   logical, intent(in), optional :: has_lsm
   real,    intent(in), optional :: regime_wrf, qsfc_wrf, znt_wrf, ust_wrf, mol_wrf
   real,    intent(in), optional :: hfx, qfx, pblh

   logical :: use_table = .true.
   logical :: use_ust_wrf = .false.
   logical :: vconv_wrf
   integer :: nn, nz, n2
   real    :: rr, rz, r2
   real    :: cqs2, chs2, rho, rhox, fluxc, visc, restar, z0t, z0q

   ! h10 is the height of 10m where the wind observed
   ! h2  is the height of 2m where the temperature and 
   !        moisture observed.

   real, parameter :: h10 = 10.0, h2 = 2.0
   
   ! Default roughness over the land

   real, parameter :: zint0 = 0.01 
   
   ! Von Karman constant

   real, parameter :: k_kar = 0.4
   
   ! Working variables

   real :: Vc2, Va2, V2, vc, wspd
   real :: rib, rcp, xx, yy, cc
   real :: psiw, psiz, mol, ust, hol, holz, hol2
   real :: psim, psimz, psim2, psih, psihz, psih2
   real :: psit, psit2, psiq, psiq2
   real :: gzsoz0, gz10oz0, gz2oz0
   real :: eg, qg, tvg, tvs, tvs2
   real :: ths, thg, thvs, thvg, thvs2, vsgd, vsgd2, dx
   real :: zq0, z0

   real, parameter :: ka = 2.4E-5

   if (trace_use_dull) call da_trace_entry("da_sfc_wtq")

   rcp = gas_constant/cp

   ! 1 Compute the roughness length based upon season and land use 

   ! 1.1 Define the roughness length

   z0 = roughness

   if (z0 < 0.0001) z0 = 0.0001

   if ( present(znt_wrf) ) then
      if ( znt_wrf > 0.0 ) then
         z0 = znt_wrf
      end if
   end if

   ! 1.2 Define the rouhgness length for moisture

   if (xland .ge. 1.5) then
      zq0 = z0
   else
      zq0 =  zint0
   end if

   ! 1.3 Define the some constant variable for psi

   gzsoz0 = log(hs/z0)

   gz10oz0 = log(h10/z0)

   gz2oz0 = log(h2/z0)


   ! 2. Calculate the virtual temperature

   ! 2.1 Compute Virtual temperature on the lowest half sigma level

   tvs  = ts * (1.0 + 0.608 * qs)

   ! 2.2 Convert ground virtual temperature assuming it's saturated

   call da_tp_to_qs(tg, psfc, eg, qg) !output qg is specific humidity
   qg = qg*(1.0-qg) !hcl convert to mixing ratio
   if ( present(qsfc_wrf) ) then
      if ( qsfc_wrf > 0.0 ) then
         qg = qsfc_wrf
      end if
   endif

   tvg  = tg * (1.0 + 0.608 * qg)

   ! 3.  Compute the potential temperature

   ! 3.1 Potential temperature on the lowest half sigma level

   ths  = ts * (1000.0 / (ps/100.0)) ** rcp

   ! 3.2 Potential temperature at the ground

   thg  = tg * (1000.0 / (psfc/100.0)) ** rcp

   ! 4. Virtual potential temperature

   ! 4.1 Virtual potential temperature on the lowest half sigma level

   thvs = tvs * (1000.0 / (ps/100.0)) ** rcp

   ! 4.2 Virtual potential temperature at ground

   thvg = tvg * (1000.0 / (psfc/100.0)) ** rcp


   ! 5.  BULK RICHARDSON NUMBER AND MONI-OBUKOV LENGTH

   ! 5.1 Velocity
   
   !     Wind speed:

   Va2 =   us*us + vs*vs
   !  
   !     Convective velocity:

   vconv_wrf = .false.
   if ( present(hfx) .and. present(qfx) .and. present(pblh) ) then
      ! calculating vconv over land following wrf method
      if ( pblh > 0.0 ) then
         vconv_wrf = .true.
      end if
   end if

   if (thvg >= thvs) then
      ! prior to V3.7, Vc2 = 4.0 * (thvg - thvs)
      Vc2 = thvg - thvs
   else
      Vc2 = 0.0
   end if
   if ( xland < 1.5 ) then !land
      if ( vconv_wrf ) then
         ! following the calculation as in module_sf_sfclay.F
         rhox = psfc/(gas_constant*tvg)
         fluxc = max(hfx/rhox/cp+0.608*tvg*qfx/rhox, 0.0)
         vc = (gravity/tg*pblh*fluxc)**0.33
         vc2 = vc*vc
      end if
   end if

   ! Calculate Mahrt and Sun low-res correction         ! Add by Eric Chiang ( July 2010 )
   vsgd = 0.32 * (max(dx/5000.-1.,0.))**0.33            ! Add by Eric Chiang ( July 2010 )
   vsgd2 = vsgd * vsgd                                  ! Add by Eric Chiang ( July 2010 )
   
   V2  = Va2 + Vc2 + vsgd2                              ! Add by Eric Chiang ( July 2010 )  
   wspd = sqrt(v2)
   wspd = max(wspd,0.1)
   v2 = wspd*wspd

   ! 5.2 Bulk richardson number

   rib = (gravity * hs / ths) * (thvs - thvg) / V2

   ! if previously unstable, do not let into regime 1 and 2
   if ( present(mol_wrf) ) then
      if ( mol_wrf < 0.0 ) rib = min(rib, 0.0)
   end if

   !  Calculate   ust, m/L (mol), h/L (hol)

   psim = 0.0
   psih = 0.0

   ! Friction speed

   if ( present(ust_wrf) ) then
      if ( ust_wrf > 0.0 ) then
         use_ust_wrf = .true.
         ust = ust_wrf
      end if
   end if
   if ( .not. use_ust_wrf ) then
      !ust = 0.0001  !init value as in phys/module_physics_init.F
      ust = k_kar * sqrt(v2) /(gzsoz0 - psim)
   end if

   ! Heat flux factor

   if ( present(mol_wrf) ) then
      mol = mol_wrf
   else
      mol = k_kar * (ths - thg)/(gzsoz0 - psih)
      !mol = 0.0
   end if

   ! set regimes based on rib
   if (rib .GE. 0.2) then
      ! Stable conditions (REGIME 1)
      regime = 1.1
   else if ((rib .LT. 0.2) .AND. (rib .GT. 0.0)) then
      ! Mechanically driven turbulence (REGIME 2)
      regime = 2.1
   else if (rib .EQ. 0.0) then
      ! Unstable Forced convection (REGIME 3)
      regime = 3.1
   else 
      ! Free convection (REGIME 4)
      regime = 4.1
   end if

   if ( present(regime_wrf) ) then
      if ( regime_wrf > 0.0 ) then
         regime = regime_wrf
      end if
   end if

   ! 6.  CALCULATE PSI BASED UPON REGIME

   !if (rib .GE. 0.2) then
   if ( nint(regime) == 1 ) then
      ! 6.1 Stable conditions (REGIME 1)
      !     ---------------------------
      regime = 1.1
      psim = -10.0*gzsoz0
      psim = max(psim,-10.0)
      psimz = h10/hs * psim
      psimz = max(psimz,-10.0)
      psim2 = h2 /hs * psim
      psim2 = max(psim2,-10.0)
      psih = psim
      psihz = psimz
      psih2 = psim2

   !else if ((rib .LT. 0.2) .AND. (rib .GT. 0.0)) then
   else if ( nint(regime) == 2 ) then

      ! 6.2 Mechanically driven turbulence (REGIME 2)

      regime = 2.1
      psim = (-5.0 * rib) * gzsoz0 / (1.1 - 5.0*rib)
      psim = max(psim,-10.0)
      psimz = h10/hs * psim
      psimz = max(psimz,-10.0)
      psim2 = h2 /hs * psim
      psim2 = max(psim2,-10.0)
      psih = psim
      psihz = psimz
      psih2 = psim2

   !else if (rib .EQ. 0.0) then
   else if ( nint(regime) == 3 ) then
      ! 6.3 Unstable Forced convection (REGIME 3)

      regime = 3.1
      psim = 0.0
      psimz = 0.0
      psim2 = 0.0
      psih = psim
      psihz = psimz
      psih2 = psim2

   else 
      ! 6.4 Free convection (REGIME 4)
      regime = 4.1

      cc = 2.0 * atan(1.0)

      ! Ratio of PBL height to Monin-Obukhov length

      if (ust .LT. 0.01) then
         hol = rib * gzsoz0
      else
         hol = k_kar * gravity * hs * mol / (ths * ust * ust)
      end if

      ! 6.4.2  Calculate n, nz, R, Rz

      holz = (h10 / hs) * hol
      hol2 = (h2 / hs) * hol

      hol = min(hol,0.0)
      hol = max(hol,-9.9999)

      holz = min(holz,0.0)
      holz = max(holz,-9.9999)

      hol2 = min(hol2,0.0)
      hol2 = max(hol2,-9.9999)

      ! 6.4.3 Calculate Psim & psih

      if ( use_table ) then
         ! Using the look-up table:
         nn = int(-100.0 * hol)
         rr = (-100.0 * hol) - nn
         psim = psimtb(nn) + rr * (psimtb(nn+1) - psimtb(nn))
         psih = psihtb(nn) + rr * (psihtb(nn+1) - psihtb(nn))
      else
         ! Using the continuous function:
         xx = (1.0 - 16.0 * hol) ** 0.25
         yy = log((1.0+xx*xx)/2.0)
         psim = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
         psih = 2.0 * yy
      end if

      if ( use_table ) then
         ! Using the look-up table:
         nz = int(-100.0 * holz)
         rz = (-100.0 * holz) - nz
         psimz = psimtb(nz) + rz * (psimtb(nz+1) - psimtb(nz))
         psihz = psihtb(nz) + rz * (psihtb(nz+1) - psihtb(nz))
      else
         ! Using the continuous function:
         xx = (1.0 - 16.0 * holz) ** 0.25
         yy = log((1.0+xx*xx)/2.0)
         psimz = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
         psihz = 2.0 * yy
      end if

      if ( use_table ) then
         ! Using the look-up table:
         n2 = int(-100.0 * hol2)
         r2 = (-100.0 * hol2) - n2
         psim2 = psimtb(n2) + r2 * (psimtb(n2+1) - psimtb(n2))
         psih2 = psihtb(n2) + r2 * (psihtb(n2+1) - psihtb(n2))
      else
         ! Using the continuous function:
         xx = (1.0 - 16.0 * hol2) ** 0.25
         yy = log((1.0+xx*xx)/2.0)
         psim2 = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
         psih2 = 2.0 * yy
      end if

      ! 6.4.4 Define the limit value for psim & psih

      psim = min(psim,0.9*gzsoz0)
      psimz = min(psimz,0.9*gz10oz0)
      psim2 = min(psim2,0.9*gz2oz0)
      psih = min(psih,0.9*gzsoz0)
      psihz = min(psihz,0.9*gz10oz0)
      psih2 = min(psih2,0.9*gz2oz0)
   end if  ! Regime

   ! 7.  Calculate psi for wind, temperature and moisture

   psiw = gzsoz0 - psim
   psiz = gz10oz0 - psimz
   psit = max(gzsoz0-psih, 2.0)
   psit2 = gz2oz0 - psih2

   if ( .not. use_ust_wrf ) then
      ! re-calculate ust since psim is now available
      ust = k_kar * sqrt(v2) /(gzsoz0 - psim)
   end if

   psiq  = log(k_kar*ust*hs/ka + hs / zq0) - psih
   psiq2 = log(k_kar*ust*h2/ka + h2 / zq0) - psih2

   !V3.7, as in module_sf_sfclay.F
   if ( xland >= 1.5 ) then !water
      visc = (1.32+0.009*(ts-273.15))*1.e-5
      restar = ust*z0/visc
      z0t = (5.5e-5)*(restar**(-0.60))
      z0t = min(z0t,1.0e-4)
      z0t = max(z0t,2.0e-9)
      z0q = z0t
      psiq  = max(log((hs+z0q)/z0q)-psih,  2.)
      psit  = max(log((hs+z0t)/z0t)-psih,  2.)
      psiq2 = max(log((2.+z0q)/z0q)-psih2, 2.)
      psit2 = max(log((2.+z0t)/z0t)-psih2, 2.)
   end if

   ! 8.  Calculate 10m wind, 2m temperature and moisture

   u10 = us * psiz / psiw
   v10 = vs * psiz / psiw
   t2 = (thg + (ths - thg)*psit2/psit)*((psfc/100.0)/1000.0)**rcp
   q2 = qg + (qs - qg)*psiq2/psiq 

   if ( present(has_lsm) ) then
      if ( has_lsm ) then
         !cqs2: 2m surface exchange coefficient for moisture
         !chs2: 2m surface exchange coefficient for heat
         cqs2 = ust*k_kar/psiq2
         if (xland .ge. 1.5) then
            !water
            chs2 = ust*k_kar/psit2
         else
            !land
            chs2 = cqs2 !as in subroutine lsm in phys/module_sf_noahdrv.F
         end if

         !re-calculate T2/Q2 as in module_sf_sfcdiags.F
         rho  = psfc/(gas_constant*tg)
         if ( cqs2 < 1.e-5 ) then
            q2 = qg
         else
            if ( present(qfx) ) then
               q2 = qg - qfx/(rho*cqs2)
            end if
         end if
         if ( chs2 < 1.e-5 ) then
            t2 = tg
         else
            if ( present(hfx) ) then
               t2 = tg - hfx/(rho*cp*chs2)
            end if
         end if
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_sfc_wtq")

end subroutine da_sfc_wtq


subroutine da_sfc_wtq_lin(psfc, tg, ps, ts, qs, us, vs, regime,           &
   psfc_prime,tg_prime,ps_prime,ts_prime,qs_prime, &
   us_prime,vs_prime, hs,roughness,xland,dx,         & ! Modified by Eric Chiang (JULY 2010)
   u10_prime,v10_prime,t2_prime,q2_prime) 

   !---------------------------------------------------------------------------
   ! Purpose: Calculate the  10m wind, 2m temperature and moisture based on the
   ! similarity theory/
   !
   ! Reference:
   ! ---------
   !
   !  input Variables(basic state):
   ! 
   !   psfc, tg               : surface pressure and ground temperature
   !   ps, ts, qs, us, vs, hs : model variable at lowlest half sigma level
   !   regime                 : PBL regime
   !
   !  input Variables(pertubation):
   ! 
   !   psfc_prime, tg_prime   : Surface pressure and ground temperature
   !   ps_prime, ts_prime,    : Model variables at the lowest half sigma
   !   qs_prime, us_prime,    : level 
   !   vs_prime               : 
   !
   !  Constants:
   !
   !   hs                     : height at the lowest half sigma level
   !   roughness              : roughness
   !   xland                  : land-water-mask (=2 water, =1 land)
   !
   !  output Variables(pertubation):
   !
   !   u10_prime, v10_prime   : 10-m high observed wind components
   !   t2_prime , q2_prime    : 2-m high observed temperature and mixing ratio
   !
   !---------------------------------------------------------------------------
   !  
   !                      psim  : mechanical psi at lowlest sigma level
   !                      psim2 : mechanical psi at 2m 
   !                      psimz : mechanical psi at 10m 
   !
   !---------------------------------------------------------------------------

   implicit none

   real, intent (in)  :: regime
   real, intent (in)  :: ps , ts , qs , us, vs, psfc, tg
   real, intent (in)  :: ps_prime, ts_prime, qs_prime, us_prime, vs_prime, psfc_prime, tg_prime
   real, intent (in)  :: hs, roughness, xland
   real, intent (out) :: u10_prime, v10_prime, t2_prime, q2_prime

   ! Maximum number of iterations in computing psim, psih

   integer, parameter :: k_iteration = 10 
   !      integer, parameter :: k_iteration = 1

   ! h10 is the height of 10m where the wind observed
   ! h2  is the height of 2m where the temperature and 
   !        moisture observed.

   real, parameter :: h10 = 10.0, h2 = 2.0
   !
   ! Default roughness over the land

   real, parameter :: zint0 = 0.01 
   !
   ! Von Karman constant

   real, parameter :: k_kar = 0.4
   !
   ! Working variables

   real :: Vc2, Va2, V2 
   real :: rib, rcp, xx, yy, cc, Pi
   real :: psiw, psiz, mol, ust, hol, holz, hol2
   real :: psim, psimz, psim2, psih, psihz, psih2
   real :: psit, psit2, psiq, psiq2
   real :: gzsoz0, gz10oz0, gz2oz0
   real :: eg, qg, tvg, tvs
   real :: ths, thg, thvs, thvg, vsgd, vsgd2, dx
   real :: zq0, z0

   real :: Vc2_prime, Va2_prime, V2_prime
   real :: rib_prime, xx_prime, yy_prime
   real :: psiw_prime, psiz_prime, mol_prime, ust_prime, &
           hol_prime, holz_prime, hol2_prime
   real :: psim_prime, psimz_prime, psim2_prime, &
           psih_prime, psihz_prime, psih2_prime
   real :: psit_prime, psit2_prime, &
           psiq_prime, psiq2_prime
   real :: qg_prime, tvg_prime, tvs_prime
   real :: ths_prime, thg_prime, thvs_prime, thvg_prime 

   real, parameter :: ka = 2.4E-5

   integer :: iregime

   if (trace_use) call da_trace_entry("da_sfc_wtq_lin")

   rcp = gas_constant/cp

   ! 1 Compute the roughness length based upon season and land use 
   ! =====================================

   ! 1.1 Define the rouhness length
   !     -----------------

   z0 = roughness

   if (z0 < 0.0001) z0 = 0.0001

   ! 1.2 Define the rouhgness length for moisture
   !     -----------------

   if (xland .ge. 1.5) then
      zq0 = z0
    else
      zq0 =  zint0
   end if

   ! 1.3 Define the some constant variable for psi
   !     -----------------

   gzsoz0 = log(hs/z0)

   gz10oz0 = log(h10/z0)

   gz2oz0 = log(h2/z0)


   ! 2. Calculate the virtual temperature
   ! =====================================

   ! 2.1 Compute Virtual temperature on the lowest half sigma level
   !     ---------------------------------------------------------

   tvs_prime  = ts_prime * (1.0 + 0.608 * qs) + 0.608 * ts * qs_prime
   tvs  = ts * (1.0 + 0.608 * qs)

   ! 2.2 Compute the ground saturated mixing ratio and the ground virtual 
   !     temperature
   !     ----------------------------------------------------------------

   call da_tp_to_qs(tg, psfc, eg, qg)
   call da_tp_to_qs_lin1(tg, psfc, eg, tg_prime, psfc_prime, qg_prime)

   qg_prime = qg_prime * qg

   tvg_prime  = tg_prime * (1.0 + 0.608 * qg) + 0.608 * tg * qg_prime
   tvg  = tg * (1.0 + 0.608 * qg)

   ! 3.  Compute the potential temperature and virtual potential temperature
   ! =======================================================================

   ! 3.1 Potential temperature on the lowest half sigma level
   !     ----------------------------------------------------

   Pi = (100000.0 / ps) ** rcp
   ths_prime  = (ts_prime - ps_prime * rcp * ts/ps) * Pi 
   ths  = ts * Pi

   ! 3.2 Virtual potential temperature on the lowest half sigma level
   !     ------------------------------------------------------------

   thvs_prime  = (tvs_prime - ps_prime * rcp * tvs/ps) * Pi 
   thvs = tvs * Pi

   ! 3.3 Potential temperature at the ground
   !     -----------------------------------

   Pi = (100000.0 / psfc) ** rcp
   thg_prime  = (tg_prime - psfc_prime * rcp * tg/psfc) * Pi 
   thg  = tg * Pi

   ! 3.4 Virtual potential temperature at ground
   !     ---------------------------------------

   thvg_prime  = (tvg_prime - psfc_prime * rcp * tvg/psfc) * Pi
   thvg = tvg * Pi

   ! 4.  BULK RICHARDSON NUMBER AND MONI-OBUKOV LENGTH
   ! =================================================

   ! 4.1 Velocity
   !     --------
   
   ! Wind speed:

   Va2_prime =  2.0*us*us_prime + 2.0*vs*vs_prime
   Va2 =   us*us + vs*vs
    
   ! Convective velocity:

   if (thvg >= thvs) then
      Vc2_prime = 4.0 * (thvg_prime - thvs_prime)
      Vc2 = 4.0 * (thvg - thvs)
   else
      Vc2_prime = 0.0
      Vc2 = 0.0
   end if
 
   ! Calculate Mahrt and Sun low-res correction                    ! Add by Eric Chiang ( July 2010 )
   ! dx is a constant, so vsgd is also a constant, 
   ! the perturnations of vsgd and vsgd2 are zero. (YRG, 09/15/2011)
                                     
   vsgd = 0.32 * (max(dx/5000.-1.,0.))**0.33                       ! Add by Eric Chiang ( July 2010 )
   vsgd2 = vsgd * vsgd                                             ! Add by Eric Chiang ( July 2010 )

   ! V2_prime should be computed before used below. (YRG, 09/15/2011)    
   V2_prime = Va2_prime + Vc2_prime
   V2  = Va2 + Vc2 + vsgd2                                         ! Modified by Eric Chiang ( July 2010 )

   ! 4.2 Bulk richardson number
   !     ----------------------

   Pi = gravity * hs / (ths*V2)
   rib_prime = (thvs_prime - thvg_prime   &
              - (thvs-thvg)/V2  * V2_prime &
              - (thvs-thvg)/ths * ths_prime) * Pi 
   rib = (thvs - thvg) * Pi

   ! 5.  CALCULATE PSI BASED UPON REGIME
   ! =======================================

   iregime = int(regime)

   select case (iregime) 

   ! 5.1 Stable conditions (REGIME 1)
   !     ---------------------------

   case (1);

      psim_prime  = 0.0
      psimz_prime = 0.0
      psim2_prime = 0.0
      psim  = -10.0*gzsoz0
      psimz = -10.0*gz10oz0
      psim2 = -10.0*gz2oz0
      psim  = max(psim,-10.0)
      psimz = max(psimz,-10.0)
      psim2 = max(psim2,-10.0)

      psih_prime  = psim_prime
      psihz_prime = psimz_prime
      psih2_prime = psim2_prime
      psih  = psim
      psihz = psimz
      psih2 = psim2

   ! 5.2 Mechanically driven turbulence (REGIME 2)
   !     ------------------------------------------

   case (2);

      Pi =  - 1.0 / ((1.1 - 5.0*rib)*(1.1 - 5.0*rib))
      psim_prime  = 5.5 * gzsoz0  * rib_prime * Pi 
      psimz_prime = 5.5 * gz10oz0 * rib_prime * Pi
      psim2_prime = 5.5 * gz2oz0  * rib_prime * Pi

      Pi =  (-5.0 * rib) / (1.1 - 5.0*rib)
      psim  = gzsoz0  * Pi
      psimz = gz10oz0 * Pi
      psim2 = gz2oz0  * Pi

      if (psim >= -10.0) then
         psim = psim
         psim_prime = psim_prime
      else
         psim = -10.0
         psim_prime = 0.0
      end if
      if (psimz >= -10.0) then
         psimz = psimz
         psimz_prime = psimz_prime
      else
         psimz = -10.0
         psimz_prime = 0.0
      end if
      if (psim2 >= -10.0) then
         psim2 = psim2
         psim2_prime = psim2_prime
      else
         psim2 = -10.0
         psim2_prime = 0.0
      end if

      psih_prime  = psim_prime
      psihz_prime = psimz_prime
      psih2_prime = psim2_prime
      psih = psim
      psihz = psimz
      psih2 = psim2

   ! 5.3 Unstable Forced convection (REGIME 3)
   !     -------------------------------------

   case (3);

      psim_prime  = 0.0
      psimz_prime = 0.0
      psim2_prime = 0.0

      psim  = 0.0
      psimz = 0.0
      psim2 = 0.0

      psih_prime  = psim_prime
      psihz_prime = psimz_prime
      psih2_prime = psim2_prime

      psih  = psim
      psihz = psimz
      psih2 = psim2


      ! 5.4 Free convection (REGIME 4)
      !     --------------------------

   case (4);

      ! Calculate psi m and pshi h using iteration method

      psim_prime = 0.0
      psih_prime = 0.0
      psim = 0.0
      psih = 0.0
      cc = 2.0 * atan(1.0)

      ! do k = 1 , k_iteration

      ! 5.4.1  Calculate   ust, m/L (mol), h/L (hol)
      !        --------------------------

      ! Friction speed

      ust = k_kar * sqrt(v2) /(gzsoz0 - psim)
      ust_prime = (0.5/V2 * v2_prime + psim_prime /(gzsoz0 - psim)) * ust

      ! Heat fux factor

      mol = k_kar * (ths - thg)/(gzsoz0 - psih)
      mol_prime = ((ths_prime - thg_prime) /(ths - thg) + &
                     psih_prime /(gzsoz0 - psih)) * mol

      ! Ratio of PBL height to Monin-Obukhov length

      if (ust .LT. 0.01) then
         hol_prime = rib_prime * gzsoz0
         hol = rib * gzsoz0
      else
         hol = k_kar * gravity * hs * mol / (ths * ust * ust)
         hol_prime = (mol_prime / mol - ths_prime / ths &
                        - 2.0* ust_prime / ust) * hol
      end if

      ! 5.4.2  Calculate n, nz, R, Rz
      !        --------------------------

      if (hol >= 0.0) then
         hol_prime = 0.0
         hol = 0.0
      else
         hol_prime = hol_prime
         hol = hol
      end if
      if (hol >= -10.0) then
         hol_prime = hol_prime
         hol = hol
      else
         hol_prime = 0.0
         hol = -10.0
      end if

      holz_prime = (h10 / hs) * hol_prime
      holz = (h10 / hs) * hol
      if (holz >= 0.0) then
         holz_prime = 0.0
         holz = 0.0
      else
         holz_prime = holz_prime
         holz = holz
      end if
      if (holz >= -10.0) then
         holz_prime = holz_prime
         holz = holz
      else
         holz_prime = 0.0
         holz = -10.0
      end if

      hol2_prime = (h2 / hs) * hol_prime
      hol2 = (h2 / hs) * hol
      if (hol2 >= 0.0) then
         hol2_prime = 0.0
         hol2 = 0.0
      else
         hol2_prime = hol2_prime
         hol2 = hol2
      end if
      if (hol2 >= -10.0) then
         hol2_prime = hol2_prime
         hol2 = hol2
      else
         hol2_prime = 0.0
         hol2 = -10.0
      end if

      ! 5.4.3 Calculate Psim & psih
      !        --------------------------

      ! Using the continuous function:
      xx_prime = -4.0* hol_prime /((1.0 - 16.0 * hol) ** 0.75)
      xx = (1.0 - 16.0 * hol) ** 0.25
      yy_prime = 2.0* xx * xx_prime /(1.0+xx*xx)
      yy = log((1.0+xx*xx)/2.0)
      psim_prime = 2 * xx_prime *(1.0/(1.0+xx)-1.0/(1+xx*xx)) + yy_prime 
      psim = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
      psih_prime = 2.0 * yy_prime
      psih = 2.0 * yy

      ! Using the continuous function:
      xx_prime = -4.0* holz_prime /((1.0 - 16.0 * holz) ** 0.75)
      xx = (1.0 - 16.0 * holz) ** 0.25
      yy_prime = 2.0* xx * xx_prime /(1.0+xx*xx)
      yy = log((1.0+xx*xx)/2.0)
      psimz_prime = 2.0* xx_prime *(1.0/(1.0+xx)-1.0/(1+xx*xx)) + yy_prime
      psimz = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
      psihz_prime = 2.0 * yy_prime
      psihz = 2.0 * yy

      ! Using the continuous function:
      xx_prime = -4.0* hol2_prime /((1.0 - 16.0 * hol2) ** 0.75)
      xx = (1.0 - 16.0 * hol2) ** 0.25
      yy_prime = 2.0* xx * xx_prime /(1.0+xx*xx)
      yy = log((1.0+xx*xx)/2.0)
      psim2_prime = 2.0* xx_prime *(1.0/(1.0+xx)-1.0/(1+xx*xx)) + yy_prime
      psim2 = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
      psih2_prime = 2.0 * yy_prime
      psih2 = 2.0 * yy

      ! end do 

      ! 5.4.4 Define the limit value for psim & psih
      !        --------------------------

      if (psim <= 0.9*gzsoz0) then
         psim_prime = psim_prime
         psim = psim
      else
         psim = 0.9*gzsoz0
         psim_prime = 0.0
      end if
      if (psimz <= 0.9*gz10oz0) then
         psimz_prime = psimz_prime
         psimz = psimz
      else
         psimz_prime = 0.0
         psimz = 0.9*gz10oz0
      end if
      if (psim2 <= 0.9*gz2oz0) then
         psim2_prime = psim2_prime
         psim2 = psim2
      else
         psim2_prime = 0.0
         psim2 = 0.9*gz2oz0
      end if
      if (psih <= 0.9*gzsoz0) then
         psih_prime = psih_prime
         psih = psih
      else
         psih_prime = 0.0
         psih = 0.9*gzsoz0
      end if
      if (psihz <= 0.9*gz10oz0) then
         psihz_prime = psihz_prime
         psihz = psihz
      else
         psihz_prime = 0.0
         psihz = 0.9*gz10oz0
      end if
      if (psih2 <= 0.9*gz2oz0) then
         psih2_prime = psih2_prime
         psih2 = psih2
      else
         psih2_prime = 0.0
         psih2 = 0.9*gz2oz0
       end if

    case default;
       write(unit=message(1),fmt='(A,I2,A)') &
          "Regime=",iregime," is invalid."
       call da_error("da_sfc_wtq_lin.inc",494,message(1:1))

   end select

   ! 6.  CALCULATE PSI FOR WinD, TEMPERATURE AND MOISTURE
   ! =======================================

   psiw_prime = - psim_prime
   psiw = gzsoz0 - psim
   psiz_prime = - psimz_prime
   psiz = gz10oz0 - psimz
   psit_prime = - psih_prime
   psit = gzsoz0 - psih
   psit2_prime = - psih2_prime
   psit2 = gz2oz0 - psih2

   ust = k_kar * sqrt(v2) /(gzsoz0 - psim)
   ust_prime = (0.5/V2 * v2_prime + psim_prime /(gzsoz0 - psim)) * ust

   psiq2_prime = k_kar*hs/(ka*(k_kar*ust*hs/ka + hs / zq0))*ust_prime
   psiq_prime  = psiq2_prime - psih_prime
   psiq2_prime = psiq2_prime - psih2_prime

   psiq  = log(k_kar*ust*hs/ka + hs / zq0) - psih
   psiq2 = log(k_kar*ust*h2/ka + h2 / zq0) - psih2

   ! 7.  CALCULATE THE PERTURBATIONS for 10M WinD, 2M TEMPERATURE AND MOISTURE
   ! =======================================

   Pi = psiz / psiw
   u10_prime= (us_prime + us/psiz * psiz_prime - us/psiw * psiw_prime) * Pi
   v10_prime= (vs_prime + vs/psiz * psiz_prime - vs/psiw * psiw_prime) * Pi

   t2_prime = ((1.0-psit2/psit) * thg_prime + (ths_prime + &
                           (ths - thg)/psit2 * psit2_prime - &
                           (ths - thg)/psit  * psit_prime) * psit2/psit &
             + rcp*(thg + (ths - thg)*psit2/psit)/psfc * psfc_prime) &
             * (psfc/100000.0)**rcp

   q2_prime = (1.0-psiq2/psiq) * qg_prime + psiq2/psiq * qs_prime + &
              (qs -qg)*(psiq2/psiq) * (psiq2_prime/psiq2 - psiq_prime/psiq)

   if (trace_use) call da_trace_exit("da_sfc_wtq_lin")

end subroutine da_sfc_wtq_lin



subroutine da_sfc_wtq_adj (psfc, tg, ps, ts, qs, us, vs, regime,& 
   psfc_prime, tg_prime, ps_prime, ts_prime, qs_prime, &
   us_prime, vs_prime,  hs, roughness, xland, dx,      & 
   u10_prime, v10_prime, t2_prime, q2_prime) 

   !---------------------------------------------------------------------------
   ! Purpose: Calculate the  10m wind, 2m temperature and moisture based on the
   ! similarity theory
   !
   ! Reference:
   ! ---------
   !
   !  input Variables(basic state):
   ! 
   !   psfc, tg               : surface pressure and ground temperature
   !   ps, ts, qs, us, vs, hs : model variable at lowlest half sigma level
   !   regime                 : PBL regime
   !
   !  input Variables(pertubation):
   ! 
   !   psfc_prime, tg_prime   : Surface pressure and ground temperature
   !   ps_prime, ts_prime,    : Model variables at the lowest half sigma
   !   qs_prime, us_prime,    : level 
   !   vs_prime               : 
   !
   !  Constants:
   !
   !   hs                     : height at the lowest half sigma level
   !   roughness              : roughness
   !   xland                  : land-water-mask (=2 water, =1 land)
   !
   !  output Variables(pertubation):
   !
   !   u10_prime, v10_prime   : 10-m high observed wind components
   !   t2_prime , q2_prime    : 2-m high observed temperature and mixing ratio
   !
   !---------------------------------------------------------------------------
   !  
   !                      psim  : mechanical psi at lowlest sigma level
   !                      psim2 : mechanical psi at 2m 
   !                      psimz : mechanical psi at 10m 
   !
   !---------------------------------------------------------------------------

   implicit none

   real, intent (in)          :: regime
   real, intent (in)          :: ps , ts , qs , us, vs, psfc, tg
   real, intent (in)          :: hs, roughness, xland
   real, intent (in)          :: u10_prime, v10_prime, t2_prime, q2_prime
   real, intent (inout)       :: ps_prime, ts_prime, qs_prime, us_prime, vs_prime, psfc_prime, tg_prime

   ! Maximum number of iterations in computing psim, psih

   integer, parameter :: k_iteration = 10 
   !      integer, parameter :: k_iteration = 1

   ! h10 is the height of 10m where the wind observed
   ! h2  is the height of 2m where the temperature and 
   !        moisture observed.

   real, parameter :: h10 = 10.0, h2 = 2.0
   !
   ! Default roughness over the land

   real, parameter :: zint0 = 0.01 
   !
   ! Von Karman constant

   real, parameter :: k_kar = 0.4
   !
   ! Working variables

   real :: Vc2, Va2, V2 
   real :: rib, rcp, xx, yy, cc, Pi
   real :: psiw, psiz, mol, ust, hol, holz, hol2
   real :: psim, psimz, psim2, psih, psihz, psih2
   real :: psit,  psit2, psiq, psiq2
   real :: gzsoz0, gz10oz0, gz2oz0
   real :: eg, qg, tvg, tvs
   real :: ths, thg, thvs, thvg
   real :: vsgd, vsgd2, dx                 ! Add by Eric Chiang (AUG 2010)
   real :: zq0, z0
   real :: ust_s, hol_s, psim_s, psim2_s, psimz_s, psih_s, psih2_s, psihz_s

   real :: Vc2_prime, Va2_prime, V2_prime
   real :: rib_prime, xx_prime, yy_prime
   real :: psiw_prime, psiz_prime, mol_prime, ust_prime, &
           hol_prime, holz_prime, hol2_prime
   real :: psim_prime, psimz_prime, psim2_prime, &
           psih_prime, psihz_prime, psih2_prime
   real :: psit_prime, psit2_prime, &
           psiq_prime, psiq2_prime
   real :: qg_prime, tvg_prime, tvs_prime
   real :: ths_prime, thg_prime, thvs_prime, thvg_prime

   real, parameter :: ka = 2.4E-5

   integer :: iregime

   if (trace_use) call da_trace_entry("da_sfc_wtq_adj")

   !-----------------------------------------------------------------------------!
   !  initialize perturbation value
   !------- ----------------------------------------------------------------------!
   !        tg_prime = 0.0
   !        us_prime = 0.0
   !        vs_prime = 0.0
   !        ts_prime = 0.0
   !        ps_prime = 0.0
   !        qs_prime = 0.0
   !      psfc_prime = 0.0

   psim_prime = 0.0
   psimz_prime = 0.0
   psim2_prime = 0.0
   psih2_prime = 0.0
   psihz_prime = 0.0
   psih_prime = 0.0

   psiw_prime = 0.0
   psiz_prime = 0.0
   psit_prime = 0.0
   psit2_prime = 0.0
   psiq_prime = 0.0
   psiq2_prime = 0.0

   qg_prime = 0.0
   ths_prime = 0.0
   thg_prime = 0.0

   thvs_prime = 0.0
   thvg_prime = 0.0
   V2_prime = 0.0
   rib_prime = 0.0
   ust_prime = 0.0
   tvg_prime = 0.0
   tvs_prime = 0.0
   va2_prime = 0.0
   vc2_prime = 0.0

   !  +++++++++++++++++++++++++++++++++ 
   ! A.0  Calculate Basic state
   !  +++++++++++++++++++++++++++++++++ 
   rcp = gas_constant/cp

   ! 1 Compute the roughness length based upon season and land use 
   ! =====================================

   ! 1.1 Define the rouhness length
   !     -----------------

   z0 = roughness

   if (z0 < 0.0001) z0 = 0.0001

   ! 1.2 Define the rouhgness length for moisture
   !     -----------------

   if (xland .ge. 1.5) then
      zq0 = z0
   else
      zq0 =  zint0
   end if

   ! 1.3 Define the some constant variable for psi
   !     -----------------

   gzsoz0 = log(hs/z0)

   gz10oz0 = log(h10/z0)

   gz2oz0 = log(h2/z0)


   ! 2.0 Calculate the virtual temperature
   ! =====================================

   ! 2.1 Compute Virtual temperature on the lowest half sigma level
   !     ---------------------------------------------------------

   tvs  = ts * (1.0 + 0.608 * qs)

   ! 2.2 Convert ground virtual temperature assuming it's saturated
   !     ----------------------------------------------------------
   call da_tp_to_qs(tg, psfc, eg, qg)
   tvg  = tg * (1.0 + 0.608 * qg)

   ! 3.0  Compute the potential temperature and virtual potential temperature
   ! =======================================================================

   ! 3.1 Potential temperature on the lowest half sigma level
   !     ----------------------------------------------------

   Pi = (100000.0 / ps) ** rcp
   ths  = ts * Pi

   ! 3.2 Virtual potential temperature on the lowest half sigma level
   !     ------------------------------------------------------------

   thvs = tvs * Pi

   ! 3.3 Potential temperature at the ground
   !     -----------------------------------

   Pi =  (100000.0 / psfc) ** rcp
   thg  = tg * Pi

   ! 3.4 Virtual potential temperature at ground
   !     ---------------------------------------

   thvg = tvg * Pi


   ! 4.0  BULK RICHARDSON NUMBER AND MONI-OBUKOV LENGTH
   ! =================================================

   ! 4.1 Velocity
   !     --------
   
   !     Wind speed:

   Va2 =   us*us + vs*vs
     
   !     Convective velocity:

   if (thvg >= thvs) then
      Vc2 = 4.0 * (thvg - thvs)
   else
      Vc2 = 0.0
   end if

   ! Calculate Mahrt and Sun low-res correction         ! Add by Eric Chiang ( AUG 2010 )
   vsgd = 0.32 * (max(dx/5000.-1.,0.))**0.33            ! Add by Eric Chiang ( AUG 2010 )
   vsgd2 = vsgd * vsgd                                  ! Add by Eric Chiang ( AUG 2010 )
   
   V2  = Va2 + Vc2 + vsgd2                              ! Modified by Eric Chiang ( AUG 2010 ) 

   ! 4.2 Bulk richardson number
   !     ----------------------

   rib = (gravity * hs / ths) * (thvs - thvg) / V2
   ! 5.0  CALCULATE PSI BASED UPON REGIME
   ! =======================================

   iregime = int(regime)
   select case (iregime) 

   ! 5.1 Stable conditions (REGIME 1)
   !     ---------------------------

   case (1);

      psim = -10.0*gzsoz0
      psimz = -10.0*gz10oz0
      psim2 = -10.0*gz2oz0

      psim_s  = psim
      psimz_s = psimz
      psim2_s = psim2

      psim  = max(psim,-10.0)
      psimz = max(psimz,-10.0)
      psim2 = max(psim2,-10.0)

      psih = psim
      psihz = psimz
      psih2 = psim2

   ! 5.2 Mechanically driven turbulence (REGIME 2)
   !     ------------------------------------------

   case (2);
      Pi =  (-5.0 * rib) / (1.1 - 5.0*rib)
      psim  = gzsoz0  * Pi
      psimz = gz10oz0 * Pi
      psim2 = gz2oz0  * Pi
      psim_s  = psim
      psimz_s = psimz
      psim2_s = psim2
      if (psim >= -10.0) then
          psim = psim
      else
         psim = -10.0
      end if
      if (psimz >= -10.0) then
         psimz = psimz
      else
         psimz = -10.0
      end if
      if (psim2 >= -10.0) then
         psim2 = psim2
      else
         psim2 = -10.0
      end if

      psih = psim
      psihz = psimz
      psih2 = psim2

   ! 5.3 Unstable Forced convection (REGIME 3)
   !     -------------------------------------

   case (3);

      psim = 0.0
      psimz = 0.0
      psim2 = 0.0
      psih = psim
      psihz = psimz
      psih2 = psim2


   ! 5.4 Free convection (REGIME 4)
   !     --------------------------

   case (4);

      ! Calculate psi m and pshi h using iteration method

      psim = 0.0
      psih = 0.0
      cc = 2.0 * atan(1.0)
      !
      !        do k = 1 , k_iteration

      ! 5.4.1  Calculate   ust, m/L (mol), h/L (hol)
      !        --------------------------

      ! Friction speed

      ust = k_kar * sqrt(v2) /(gzsoz0 - psim)

      ! save ust for adjoint:
      ust_s = ust

      ! Heat flux factor

      mol = k_kar * (ths - thg)/(gzsoz0 - psih)

      ! Ratio of PBL height to Monin-Obukhov length

      if (ust .LT. 0.01) then
         hol = rib * gzsoz0
      else
         hol = k_kar * gravity * hs * mol / (ths * ust * ust)
      end if

      ! 5.4.2  Calculate n, nz, R, Rz
      !        --------------------------

      hol_s = hol

      if (hol >= 0.0) then
         hol = 0.0
      else
         hol = hol
      end if
      if (hol >= -10.0) then
         hol = hol
      else
         hol = -10.0
      end if

      holz = (h10 / hs) * hol
      if (holz >= 0.0) then
         holz = 0.0
      else
         holz = holz
      end if
      if (holz >= -10.0) then
         holz = holz
      else
         holz = -10.0
      end if

      hol2 = (h2 / hs) * hol
      if (hol2 >= 0.0) then
         hol2 = 0.0
      else
         hol2 = hol2
      end if
      if (hol2 >= -10.0) then
         hol2 = hol2
      else
         hol2 = -10.0
      end if

      ! 5.4.3 Calculate Psim & psih
      !        --------------------------

      ! Using the continuous function:
      xx = (1.0 - 16.0 * hol) ** 0.25
      yy = log((1.0+xx*xx)/2.0)
      psim = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
      psih = 2.0 * yy

      ! Using the continuous function:
      xx = (1.0 - 16.0 * holz) ** 0.25
      yy = log((1.0+xx*xx)/2.0)
      psimz = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
      psihz = 2.0 * yy

      ! Using the continuous function:
      xx = (1.0 - 16.0 * hol2) ** 0.25
      yy = log((1.0+xx*xx)/2.0)
      psim2 = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
      psih2 = 2.0 * yy

      ! end do 

      ! 5.4.4 Define the limit value for psim & psih
      !        --------------------------
 
      !  Save the values for adjoint:

      psim_s  = psim
      psimz_s = psimz
      psim2_s = psim2
      psih_s  = psih
      psihz_s = psihz
      psih2_s = psih2

      if (psim <= 0.9*gzsoz0) then
         psim = psim
      else
         psim = 0.9*gzsoz0
      end if
      if (psimz <= 0.9*gz10oz0) then
         psimz = psimz
      else
         psimz = 0.9*gz10oz0
      end if
      if (psim2 <= 0.9*gz2oz0) then
         psim2 = psim2
      else
         psim2 = 0.9*gz2oz0
      end if
      if (psih <= 0.9*gzsoz0) then
         psih = psih
      else
         psih = 0.9*gzsoz0
      end if
      if (psihz <= 0.9*gz10oz0) then
         psihz = psihz
      else
         psihz = 0.9*gz10oz0
      end if
      if (psih2 <= 0.9*gz2oz0) then
         psih2 = psih2
      else
         psih2 = 0.9*gz2oz0
      end if

   case default;
      write(unit=message(1),fmt='(A,I2,A)') &
         "Regime=",iregime," is invalid."
      call da_error("da_sfc_wtq_adj.inc",459,message(1:1))

   end select
   
   ! 6.0  CALCULATE PSI FOR WinD, TEMPERATURE AND MOISTURE
   ! =======================================

   psiw = gzsoz0 - psim
   psiz = gz10oz0 - psimz
   psit = gzsoz0 - psih
   psit2 = gz2oz0 - psih2
   ust = k_kar * sqrt(v2) /(gzsoz0 - psim)
   psiq  = log(k_kar*ust*hs/ka + hs / zq0) - psih
   psiq2 = log(k_kar*ust*h2/ka + h2 / zq0) - psih2
   !  +++++++++++++++++++++++++++++++++ 
   !  B.0  Calculate adjoint solution
   !  +++++++++++++++++++++++++++++++++ 

   ! 7.0  CALCULATE 10M WinD, 2M TEMPERATURE AND MOISTURE
   ! =======================================

   qg_prime    = (1.0 - psiq2/psiq) * q2_prime
   qs_prime    = qs_prime + psiq2/psiq * q2_prime
   psiq2_prime =  (qs -qg)/psiq * q2_prime
   psiq_prime  = -(qs -qg)*psiq2/(psiq*psiq) * q2_prime
   ! q2_prime = 0.0
   
   Pi = (psfc/100000.0)**rcp
   thg_prime   = (1.0 - psit2/psit) * Pi * t2_prime
   ths_prime   = (psit2/psit) * Pi * t2_prime
   psit2_prime =   (ths - thg)/psit * Pi * t2_prime
   psit_prime  = - (ths - thg)*psit2/(psit*psit) * Pi * t2_prime
   psfc_prime  = psfc_prime + Pi * rcp*(thg + (ths - thg)*psit2/psit)/psfc * t2_prime 
   ! t2_prime = 0.0
   
   Pi = psiz / psiw
   vs_prime   = vs_prime + Pi * v10_prime
   psiz_prime =   vs / psiz * Pi * v10_prime
   psiw_prime = - vs / psiw * Pi * v10_prime
   ! v10_prime = 0.0 
   
   us_prime   = us_prime + Pi * u10_prime
   psiz_prime =  psiz_prime + us / psiz * Pi * u10_prime
   psiw_prime =  psiw_prime - us / psiw * Pi * u10_prime
   ! u10_prime = 0.0
   
   ! 6.0  CALCULATE PSI FOR WinD, TEMPERATURE AND MOISTURE
   ! =======================================
   
   ! moisture part:
   psih2_prime = - psiq2_prime
   psih_prime  = - psiq_prime
   psiq2_prime = psiq2_prime + psiq_prime
   ust_prime   = k_kar*hs/(ka*(k_kar*ust*hs/ka + hs / zq0)) * psiq2_prime 
   
   v2_prime   = 0.5/V2 * ust * ust_prime
   psim_prime = ust / (gzsoz0 - psim) * ust_prime
   ust_prime = 0.0

   ! temperature part:
   psih2_prime = psih2_prime - psit2_prime
   psih_prime  = psih_prime  - psit_prime

   ! wind part:
   psimz_prime = psimz_prime - psiz_prime
   psim_prime  = psim_prime  - psiw_prime

   ! 5.0  CALCULATE PSI BASED UPON REGIME
   ! =======================================
   select case (iregime) 

   ! 5.1 Stable conditions (REGIME 1)
   !     ---------------------------

   case (1);

      psim2_prime = psim2_prime + psih2_prime
      psimz_prime = psimz_prime + psihz_prime
      psim_prime  = psim_prime  + psih_prime
      psim_prime  = 0.0
      psimz_prime = 0.0
      psim2_prime = 0.0

   ! 5.2 Mechanically driven turbulence (REGIME 2)
   !      ------------------------------------------

   case (2);

     psim2_prime = psim2_prime + psih2_prime
     psimz_prime = psimz_prime + psihz_prime
     psim_prime  = psim_prime  + psih_prime

     psim  = psim_s
     psimz = psimz_s
     psim2 = psim2_S

     if (psim2 >= -10.0) then
        psim2_prime = psim2_prime
     else
        psim2_prime = 0.0
     end if
     if (psimz >= -10.0) then
        psimz_prime = psimz_prime
     else
        psimz_prime = 0.0
     end if
     if (psim >= -10.0) then
        psim_prime = psim_prime
     else
        psim_prime = 0.0
     end if

     Pi = -1.0 / ((1.1 - 5.*rib)*(1.1 - 5.0*rib))
     rib_prime =             5.5 * gz2oz0  * psim2_prime * Pi
     rib_prime = rib_prime + 5.5 * gz10oz0 * psimz_prime * Pi
     rib_prime = rib_prime + 5.5 * gzsoz0  * psim_prime  * Pi

     ! 5.3 Unstable Forced convection (REGIME 3)
     !     -------------------------------------

   case (3);

      psim2_prime = psih2_prime
      psimz_prime = psihz_prime
      psim_prime  = psih_prime

      psim2_prime = 0.0
      psimz_prime = 0.0
      psim_prime  = 0.0

   ! 5.4 Free convection (REGIME 4)
   !     --------------------------

   case (4);

      ! 5.4.4 Define the limit value for psim & psih
      !        -------------------------------------

      ! Recover the values:

      psim = psim_s
      psim2 = psim2_s
      psimz = psimz_s
      psihz = psihz_s
      psih  = psih_s
      psih2 = psih2_s

      if (psih2 <= 0.9*gz2oz0) then
         psih2_prime = psih2_prime
      else
         psih2_prime = 0.0
      end if
      if (psihz <= 0.9*gz10oz0) then
         psihz_prime = psihz_prime
      else
         psihz_prime = 0.0
      end if
      if (psih <= 0.9*gzsoz0) then
         psih_prime = psih_prime
      else
         psih_prime = 0.0
      end if
      if (psim2 <= 0.9*gz2oz0) then
         psim2_prime = psim2_prime
      else
         psim2_prime = 0.0
      end if
      if (psimz <= 0.9*gz10oz0) then
         psimz_prime = psimz_prime
      else
         psimz_prime = 0.0
      end if
      if (psim <= 0.9*gzsoz0) then
         psim_prime = psim_prime
      else
         psim_prime = 0.0
      end if

      ! 5.4.3 Calculate Psim & psih
      !        --------------------------

      ! Recover ust:
      ust = ust_s
      psim = 0.0
      psih = 0.0

      xx = (1.0 - 16.0 * hol2) ** 0.25
      yy = log((1.0+xx*xx)/2.0)
      yy_prime = 2.0 * psih2_prime
      psih2_prime = 0.0
   
      xx_prime = 2.0*(1.0/(1.0+xx)-1.0/(1+xx*xx)) * psim2_prime
      yy_prime = yy_prime + psim2_prime
      psim2_prime = 0.0
  
      xx_prime = xx_prime + 2.0* xx/(1.0+xx*xx) * yy_prime
      yy_prime = 0.0
      hol2_prime = -4.0/((1.0 - 16.0 * hol2) ** 0.75) * xx_prime
      xx_prime = 0.0

      xx = (1.0 - 16.0 * holz) ** 0.25
      yy = log((1.0+xx*xx)/2.0)
      yy_prime = 2.0 * psihz_prime
      psihz_prime = 0.0
      
      xx_prime = 2.0*(1.0/(1.0+xx)-1.0/(1+xx*xx)) * psimz_prime
      yy_prime = yy_prime + psimz_prime
      psimz_prime = 0.0
      
      xx_prime = xx_prime + 2.0* xx/(1.0+xx*xx) * yy_prime
      yy_prime = 0.0
      holz_prime = -4.0/((1.0 - 16.0 * holz) ** 0.75) * xx_prime
      xx_prime = 0.0

      xx = (1.0 - 16.0 * hol) ** 0.25
      yy = log((1.0+xx*xx)/2.0)
      yy_prime = 2.0 * psih_prime
      psih_prime = 0.0
      
      xx_prime = 2.0*(1.0/(1.0+xx)-1.0/(1+xx*xx)) * psim_prime
      yy_prime = yy_prime + psim_prime
      psim_prime = 0.0
      
      xx_prime = xx_prime + 2.0* xx/(1.0+xx*xx)*yy_prime
      yy_prime = 0.0
      hol_prime = -4.0/((1.0 - 16.0 * hol) ** 0.75)*xx_prime
      xx_prime = 0.0

      ! 5.4.2  Calculate n, nz, R, Rz
      !        --------------------------

      !    Recover the values:

      hol = hol_s

      hol2 = (h2 / hs) * hol
      if (hol2 >= -10.0) then
         hol2_prime = hol2_prime
      else
         hol2_prime = 0.0
      end if
      if (hol2 >= 0.0) then       
         hol2_prime = 0.0
      else
         hol2_prime = hol2_prime
      end if
      
      hol_prime = hol_prime + (h2 / hs) * hol2_prime
      hol2_prime = 0.0
      
      holz = (h10 / hs) * hol
      if (holz >= -10.0) then
         holz_prime = holz_prime
      else
         holz_prime = 0.0
      end if
      if (holz >= 0.0) then
         holz_prime = 0.0
      else
         holz_prime = holz_prime
      end if
      
      hol_prime = hol_prime + (h10 / hs) * holz_prime
      holz_prime = 0.0
      
      if (hol >= -10.0) then
         hol_prime = hol_prime
      else
         hol_prime = 0.0
      end if
      if (hol >= 0.0) then
         hol_prime = 0.0
      else
         hol_prime = hol_prime
      end if

      ! 5.4.1  Calculate   ust, m/L (mol), h/L (hol)
      !        --------------------------

      !       Ratio of PBL height to Monin-Obukhov length
      if (ust .LT. 0.01) then
         rib_prime = hol_prime * gzsoz0
         hol_prime = 0.0
      else
         mol_prime =        hol / mol * hol_prime
         ths_prime = ths_prime - hol / ths * hol_prime
         ust_prime = - 2.0 * hol / ust * hol_prime
         hol_prime = 0.0
      end if

      ! Heat flux factor
      ths_prime  = ths_prime + mol/(ths - thg) * mol_prime
      thg_prime  = thg_prime - mol/(ths - thg) * mol_prime
      psih_prime = psih_prime + mol/(gzsoz0 - psih) * mol_prime
      mol_prime = 0.0

      ! Friction speed

      v2_prime   = V2_prime + 0.5 * ust/V2 * ust_prime 
      psim_prime = psim_prime + ust/(gzsoz0 - psim) * ust_prime 
      ust_prime = 0.0

      ! Calculate psi m and pshi h using iteration method

      psim_prime = 0.0
      psih_prime = 0.0

   case default;
      write(unit=message(1),fmt='(A,I2,A)') &
         "Regime=",iregime," is invalid."
      call da_error("da_sfc_wtq_adj.inc",769,message(1:1))

   end select

   ! 4.0  BULK RICHARDSON NUMBER AND MONI-OBUKOV LENGTH
   ! =================================================

   ! 4.2 Bulk richardson number
   !     ----------------------

   Pi = gravity * hs / (ths*V2)
   ths_prime = ths_prime - Pi * (thvs-thvg)/ths * rib_prime
   V2_prime  = V2_prime  - Pi * (thvs-thvg)/V2 * rib_prime
   thvs_prime = thvs_prime + Pi * rib_prime
   thvg_prime = thvg_prime - Pi * rib_prime
   rib_prime = 0.0
   
   ! 4.1 Velocity
   !     --------
   
   ! Convective velocity:

   Va2_prime = V2_prime
   Vc2_prime = V2_prime

   if (thvg >= thvs) then
     thvg_prime = thvg_prime + 4.0 * Vc2_prime
     thvs_prime = thvs_prime - 4.0 * Vc2_prime
     Vc2_prime = 0.0
    else
     Vc2_prime = 0.0
   end if

   ! Wind speed:

   us_prime = us_prime + 2.0 *us * Va2_prime
   vs_prime = vs_prime + 2.0 *vs * Va2_prime
   Va2_prime = 0.0

   ! 3.0 Virtual potential temperature
   ! ==================================

   ! 3.4 Virtual potential temperature at ground
   !     ---------------------------------------

   Pi = (100000.0 / psfc) ** rcp
   tvg_prime = tvg_prime +  Pi * thvg_prime
   psfc_prime = psfc_prime - thvg_prime * rcp * tvg/psfc * Pi
   thvg_prime = 0.0

   ! 3.3 Potential temperature at the ground
   !     -----------------------------------

   tg_prime = tg_prime + Pi * thg_prime 
   psfc_prime = psfc_prime - thg_prime *rcp * tg/psfc * Pi
   thg_prime = 0.0

   ! 3.2 Virtual potential temperature on the lowest half sigma level
   !     ------------------------------------------------------------

   Pi = (100000.0 / ps) ** rcp
   tvs_prime = tvs_prime + PI * thvs_prime
   ps_prime = ps_prime - thvs_prime *rcp * tvs/ps * Pi
   thvs_prime = 0.0

   ! 3.1 Potential temperature on the lowest half sigma level
   !     ----------------------------------------------------
   ts_prime = ts_prime + Pi * ths_prime
   ps_prime = ps_prime - ths_prime *  rcp * ts/ps * Pi
   ths_prime = 0.0

   ! 2.0 Calculate the virtual temperature
   ! =====================================

   ! 2.2 Compute the ground saturated mixing ratio and the ground virtual 
   !     temperature
   !     ----------------------------------------------------------------

   tg_prime = tg_prime + tvg_prime * (1.0 + 0.608 * qg)
   qg_prime = qg_prime + tvg_prime * 0.608 * tg
   tvg_prime = 0.0

   qg_prime = qg_prime * qg
   call da_tp_to_qs_adj1(tg, psfc, eg, tg_prime, psfc_prime, qg_prime)

   ! 2.1 Compute Virtual temperature on the lowest half sigma level
   !     ---------------------------------------------------------

   ts_prime = ts_prime + tvs_prime  * (1.0 + 0.608 * qs)
   qs_prime = qs_prime + tvs_prime *  0.608 * ts
   tvs_prime = 0.0

   if (trace_use) call da_trace_exit("da_sfc_wtq_adj")

end subroutine da_sfc_wtq_adj


subroutine da_transform_xtopsfc(grid, iv, obs_index, synop, y_synop)

   !---------------------------------------------------------------------
   ! Purpose: TBD
   !---------------------------------------------------------------------

   implicit none

   type (domain),              intent(in)  :: grid
   type (iv_type),             intent(in)  :: iv
   integer,                    intent(in)  :: obs_index
   type (synop_type),          intent(in)  :: synop(:)
   type (residual_synop_type), intent(out) :: y_synop(:)

   integer :: n
   real :: to, qo
   real, allocatable :: hsm(:,:)
   real, allocatable :: tsm(:,:)
   real, allocatable :: qsm(:,:)
   real, allocatable :: psm(:,:)
   real, allocatable :: psm_prime(:,:)
   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: t(:,:)
   real, allocatable :: q(:,:)

   if (trace_use) call da_trace_entry("da_transform_xtopsfc")

   allocate (hsm(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (tsm(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (qsm(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (psm(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (psm_prime(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (u(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (v(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (t(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (q(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))

   ! 2.0 Surface assmiilation approach 2 (2-m t and q, 10-m u and v)

   call da_interp_lin_2d(grid%xa%u10,  iv%info(obs_index), 1, u)
   call da_interp_lin_2d(grid%xa%v10,  iv%info(obs_index), 1, v)
   call da_interp_lin_2d(grid%xa%psfc, iv%info(obs_index), 1, psm_prime)
   call da_interp_lin_2d(grid%xa%t2,   iv%info(obs_index), 1, t)
   call da_interp_lin_2d(grid%xa%q2,   iv%info(obs_index), 1, q)
   do n=iv%info(obs_index)%n1,iv%info(obs_index)%n2
      y_synop(n)%u=u(1,n)
      y_synop(n)%v=v(1,n)
      y_synop(n)%t=t(1,n)
      y_synop(n)%q=q(1,n)
   end do

   ! 3.2 model background surface p, t, q, h at observed site:

   call da_interp_lin_2d(grid%xb%terr, iv%info(obs_index), 1, hsm)
   call da_interp_lin_2d(grid%xb%t2,   iv%info(obs_index), 1, tsm)
   call da_interp_lin_2d(grid%xb%q2,   iv%info(obs_index), 1, qsm)
   call da_interp_lin_2d(grid%xb%psfc, iv%info(obs_index), 1, psm)

   do n=iv%info(obs_index)%n1,iv%info(obs_index)%n2
      if (synop(n)%p%qc < 0) then
         y_synop(n)%p = 0.0
      else

         ! 3.0 The pressure at the observed height: 

         ! 3.1 Constants:

          to = -888888.0
          qo = -888888.0
             
         ! Terrain height at the observed site:

         ! 3.3 perturbations t, q, p at the model surface:

         ! 4.0 Compute the perturbation of the surface pressure perturbation 
         !     at the observed height

         if (synop(n)%t%qc >= 0 .and. synop(n)%q%qc >= 0) then
            ! 4.1 Observed value = background + innovation: both t, q available:
            !     ----------------------------------------

            to = tsm(1,n) + synop(n)%t%inv
            qo = qsm(1,n) + synop(n)%q%inv
            call da_sfc_pre_lin(y_synop(n)%p, psm_prime(1,n), y_synop(n)%t, y_synop(n)%q, &
               psm(1,n), tsm(1,n), qsm(1,n), hsm(1,n), synop(n)%h, to, qo)
         else if (synop(n)%t%qc >= 0 .and. synop(n)%q%qc < 0) then

            ! 4.2 Observed value = background + innovation: only t available
            !     ----------------------------------------

            to = tsm(1,n) + synop(n)%t%inv
            call da_sfc_pre_lin(y_synop(n)%p, psm_prime(1,n), y_synop(n)%t, y_synop(n)%q, &
               psm(1,n), tsm(1,n), qsm(1,n), hsm(1,n), synop(n)%h, to)
         else
            ! 4.3 No observed t and q available:
            !     -----------------------------

            call da_sfc_pre_lin(y_synop(n)%p, psm_prime(1,n), y_synop(n)%t, y_synop(n)%q, &
               psm(1,n), tsm(1,n), qsm(1,n), hsm(1,n), synop(n)%h)
         end if
      end if
   end do

   deallocate (hsm)
   deallocate (tsm)
   deallocate (qsm)
   deallocate (psm)
   deallocate (psm_prime)
   deallocate (u)
   deallocate (v)
   deallocate (t)
   deallocate (q)

   if (trace_use) call da_trace_exit("da_transform_xtopsfc")

end subroutine da_transform_xtopsfc


subroutine da_transform_xtopsfc_adj(grid, iv, obs_index, synop, j_synop_y, jo_grad_x) 

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),              intent(in)    :: grid
   type (iv_type),             intent(in)    :: iv
   integer,                    intent(in)    :: obs_index
   type (synop_type),          intent(in)    :: synop(:)
   type (residual_synop_type), intent(inout) :: j_synop_y(:)   ! grad_y(jo)
   type (x_type),              intent(inout) :: jo_grad_x   ! grad_x(jo)
   

   integer :: n

   real :: to, qo
   real, allocatable :: hsm(:,:)
   real, allocatable :: tsm(:,:)
   real, allocatable :: qsm(:,:)
   real, allocatable :: psm(:,:)
   real, allocatable :: psm_prime(:,:)
   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: t(:,:)
   real, allocatable :: q(:,:)

   if (trace_use) call da_trace_entry("da_transform_xtopsfc_adj")

   allocate (hsm(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (tsm(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (qsm(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (psm(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (psm_prime(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (u(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (v(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (t(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))
   allocate (q(1,iv%info(obs_index)%n1:iv%info(obs_index)%n2))

   psm_prime = 0.0

   ! 2.1 Terrain height at the observed site (xj, yi):

   call da_interp_lin_2d(grid%xb%terr, iv%info(obs_index), 1, hsm)
   call da_interp_lin_2d(grid%xb%t2,   iv%info(obs_index), 1, tsm)
   call da_interp_lin_2d(grid%xb%q2,   iv%info(obs_index), 1, qsm)
   call da_interp_lin_2d(grid%xb%psfc, iv%info(obs_index), 1, psm)

   do n=iv%info(obs_index)%n1,iv%info(obs_index)%n2
      if (synop(n)%p%qc >= 0) then
         !------------------------------------------------------------------------
         ! 2.0 Calculate gradient with respect the pressure at the observed height: 
         !------------------------------------------------------------------------

         to = -888888.0
         qo = -888888.0

         ! 2.3 Zero out the adjoint variables:

         !----------------------------------------------------------------
         ! 3.0 Surface pressure gradient at the observed height
         !----------------------------------------------------------------

         ! 4.0 Get the surface pressure gradient at the model surface height (hsm)
         ! 4.1 Both observed to and qo available:
         if (synop(n)%t%qc >= 0 .and. synop(n)%q%qc >= 0) then
            to = tsm(1,n) + synop(n)%t%inv
            qo = qsm(1,n) + synop(n)%q%inv
            call da_sfc_pre_adj (j_synop_y(n)%p, psm_prime(1,n), j_synop_y(n)%t, j_synop_y(n)%q, &
               psm(1,n), tsm(1,n), qsm(1,n), hsm(1,n), synop(n)%h, to, qo)

            ! 4.2 only observed to available:
         else if (synop(n)%t%qc >= 0 .and. synop(n)%q%qc < 0) then
            to = tsm(1,n) + synop(n)%t%inv
            call da_sfc_pre_adj (j_synop_y(n)%p, psm_prime(1,n), j_synop_y(n)%t, j_synop_y(n)%q, &
               psm(1,n), tsm(1,n), qsm(1,n), hsm(1,n), synop(n)%h, to)

            ! 4.3 Both observed to and qo missing:
         else
            call da_sfc_pre_adj (j_synop_y(n)%p, psm_prime(1,n), j_synop_y(n)%t, j_synop_y(n)%q, &
               psm(1,n), tsm(1,n), qsm(1,n), hsm(1,n), synop(n)%h)
         end if
      end if
      t(1,n)=j_synop_y(n)%t
      q(1,n)=j_synop_y(n)%q
      u(1,n)=j_synop_y(n)%u
      v(1,n)=j_synop_y(n)%v
   end do

   ! 2.2 convert the jo_grad_y to the model grids (t2, q2, u10, v10, and psfc)
   call da_interp_lin_2d_adj(jo_grad_x%t2,   iv%info(obs_index), 1, t)
   call da_interp_lin_2d_adj(jo_grad_x%q2,   iv%info(obs_index), 1, q)
   call da_interp_lin_2d_adj(jo_grad_x%u10,  iv%info(obs_index), 1, u)
   call da_interp_lin_2d_adj(jo_grad_x%v10,  iv%info(obs_index), 1, v)
   call da_interp_lin_2d_adj(jo_grad_x%psfc, iv%info(obs_index), 1, psm_prime)

   deallocate (hsm)
   deallocate (tsm)
   deallocate (qsm)
   deallocate (psm)
   deallocate (psm_prime)
   deallocate (u)
   deallocate (v)
   deallocate (t)
   deallocate (q)

   if (trace_use) call da_trace_exit("da_transform_xtopsfc_adj")

end subroutine da_transform_xtopsfc_adj


subroutine da_transform_xtowtq (grid)

   !--------------------------------------------------------------------------
   ! Purpose: TBD
   !--------------------------------------------------------------------------
 
   implicit none

   type (domain), intent(inout) :: grid

   integer :: i, j
   real    :: height

   if (trace_use) call da_trace_entry("da_transform_xtowtq")

   !------------------------------------------------------------------------
   ! [1.0] Calculate surface variable(1-m wind, 2-m moisture and temperature)
   !------------------------------------------------------------------------

   !------------------------------------------------------------------------
   ! [2.0] Calculate sfc variable perturbations at the cross point
   !------------------------------------------------------------------------

   do j=jts, jte
      do i=its, ite
         grid%xa%tgrn(i,j) = 0.0
         height = grid%xb%h(i,j,kts) - grid%xb%terr(i,j)                 
         if (height <= 0.0) then
            message(1) = "Negative height found"
            write(unit=message(2),FMT='(2I6,A,F10.2,A,F10.2)') &
               i,j,' ht = ',grid%xb%h(i,j,kts) ,' terr =  ',grid%xb%terr(i,j)
            call da_error("da_transform_xtowtq.inc",32,message(1:2))
         end if
         call da_sfc_wtq_lin(grid%xb%psfc(i,j), grid%xb%tgrn(i,j), &
            grid%xb%p(i,j,kts), grid%xb%t(i,j,kts), grid%xb%q(i,j,kts), &
            grid%xb%u(i,j,kts), grid%xb%v(i,j,kts), &
            grid%xb%regime(i,j), &
            grid%xa%psfc(i,j), grid%xa%tgrn(i,j), &
            grid%xa%p(i,j,kts), grid%xa%t(i,j,kts), grid%xa%q(i,j,kts), &
            grid%xa%u(i,j,kts), grid%xa%v(i,j,kts), &
            height      , grid%xb%rough(i,j), grid%xb%xland(i,j), grid%xb%ds, & ! Modified by Eric Chiang (AUG 2010)
            grid%xa%u10(i,j), grid%xa%v10(i,j), &
            grid%xa%t2(i,j),  grid%xa%q2(i,j))
      end do
   end do

   if (trace_use) call da_trace_exit("da_transform_xtowtq")

end subroutine da_transform_xtowtq
 

subroutine da_transform_xtowtq_adj (grid)

   !--------------------------------------------------------------------------
   ! Purpose: TBD
   !--------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout)    :: grid

   integer :: i, j, is, js, ie, je
   real    :: height

   if (trace_use) call da_trace_entry("da_transform_xtowtq_adj")

   is = its
   js = jts

   ie = ite
   je = jte

   if (test_transforms) then
      is = its-1
      js = jts-1

      ie = ite+1
      je = jte+1

      if (is < ids) is = ids
      if (js < jds) js = jds

      if (ie > ide) ie = ide
      if (je > jde) je = jde
   end if

   ! Adjoint from Gridded 10-m wind and 2-m moisture and temperature
   ! to the model adjoint variables

   do j=js, je
      do i=is, ie
         grid%xa%tgrn(i,j)=0.0

         height = grid%xb%h(i,j,kts) - grid%xb%terr(i,j)                 
         if (height <= 0.0) then
            message(1) = "Negative height found"
            write(unit=message(2),FMT='(2I6,A,F10.2,A,F10.2)') &
               i,j,' ht = ',grid%xb%h(i,j,kts) ,' terr =  ',grid%xb%terr(i,j)
            call da_error("da_transform_xtowtq_adj.inc",48,message(1:2))
         end if
         call da_sfc_wtq_adj(grid%xb%psfc(i,j), grid%xb%tgrn(i,j), &
            grid%xb%p(i,j,kts), grid%xb%t(i,j,kts), grid%xb%q(i,j,kts), &
            grid%xb%u(i,j,kts), grid%xb%v(i,j,kts), &
            grid%xb%regime(i,j),  &
            grid%xa%psfc(i,j), grid%xa%tgrn(i,j), &
            grid%xa%p(i,j,kts), grid%xa%t(i,j,kts), grid%xa%q(i,j,kts), &
            grid%xa%u(i,j,kts), grid%xa%v(i,j,kts), &
            height      , grid%xb%rough(i,j), grid%xb%xland(i,j), grid%xb%ds, & ! Modified by Eric Chiang (AUG 2010) 
            grid%xa%u10(i,j),grid%xa%v10(i,j), &
            grid%xa%t2 (i,j),grid%xa%q2 (i,j))

         grid%xa%tgrn(i,j)=0.0
      end do
   end do

   if (trace_use) call da_trace_exit("da_transform_xtowtq_adj")

end subroutine da_transform_xtowtq_adj


subroutine da_sfc_pre (psfcm, psm, tsm, qsm, hsm, ho, to, qvo)

   !-----------------------------------------------------------------------
   ! Purpose: Correct pressure between two levels. 
   !
   ! Reference: make use of the hydrosatic equation:
   !
   !  P2 = P1 * exp [-G/R * (z2-z1) / (tv1 + tv2)/2)
   !
   ! Where:
   !  z1  = height at level 1
   !  z1  = height at level 2
   !  tv1 = temperature at level 1
   !  tv2 = temperature at level 2
   !  P1  = Pressure at level 1
   !  P2  = Pressure at level 2
   !---------------------------------------------------------------------------

   implicit none

   real, intent (out)   :: psfcm   ! model pressure at ho
   real, intent (in)    :: psm, tsm, qsm

   real, intent (in)           :: hsm, ho
   real, intent (in), optional :: to, qvo

   real                 :: tvo, tvsm, tv, dz, arg0, arg

   real, parameter      :: GASR =  gas_constant
   real, parameter      :: G = gravity

   if (trace_use) call da_trace_entry("da_sfc_pre")

   ! 1.  MODEL AND OBSERVATION VIRTUAL TEMPERATURE
   ! ---------------------------------------------

   tvsm = tsm  * (1.0 + 0.608 * qsm)
   if (present(to) .and. present(qvo)) then
      tvo = to  * (1.0 + 0.608 * qvo)
   else if (present(to) .and. .not.present(qvo)) then
      tvo = to
   else
      tvo = tvsm
   end if

   tv  = 0.5 * (tvsm + tvo)

   ! 2. HEIGHT DifFERENCE BEWTEEN MODEL SURFACE AND OBSERVATIONS
   ! ------------------------------------------------------------

   dz = hsm - ho
   arg0 = dz * g / gasr

   ! 3.  EXTRAPOLATE PRESSURE OBS TO MODEL SURFACE
   ! ---------------------------------------------

   arg = arg0    / tv 

   psfcm = psm * exp (arg)

   if (trace_use) call da_trace_exit("da_sfc_pre")

end subroutine da_sfc_pre


subroutine da_sfc_pre_lin (psfcm_prime, psm_prime, tsm_prime, qsm_prime, psm, tsm, qsm, hsm, ho, to, qvo)

   !---------------------------------------------------------------------------
   ! Purpose: Correct pressure between two levels. 
   !
   ! Reference: make use of the hydrosatic equation:
   !
   !  P2 = P1 * exp [-G/R * (z2-z1) / (tv1 + tv2)/2)
   !
   ! Where:
   !  z1  = height at level 1
   !  z1  = height at level 2
   !  tv1 = temperature at level 1
   !  tv2 = temperature at level 2
   !  P1  = Pressure at level 1
   !  P2  = Pressure at level 2
   !---------------------------------------------------------------------------

   implicit none

   ! Perturbation:
   real, intent (out)   :: psfcm_prime          ! model pressure at ho
   real, intent (in)    :: psm_prime, tsm_prime, qsm_prime            ! model surface p, t, q 
   ! Basic state:
   real, intent (in)    :: psm, tsm, qsm ! model pressure at ho and
                                                ! model surface p, t, q 
   real, intent (in)           :: hsm, ho
   real, intent (in), optional :: to, qvo

   ! working array:
   real                 :: tvo, tvsm, tv, dz, arg0
   real                 :: tvsm_prime, tvo_prime, tv_prime, arg, arg_prime

   real, parameter      :: GASR =  gas_constant
   real, parameter      :: G = gravity

   if (trace_use) call da_trace_entry("da_sfc_pre_lin")

   ! 1.  MODEL AND OBSERVATION VIRTUAL TEMPERATURE
   ! ---------------------------------------------

   tvsm_prime = tsm_prime * (1.0 + 0.608 * qsm) &
              + qsm_prime * tsm * 0.608
   tvsm = tsm  * (1.0 + 0.608 * qsm)

   if (present(to) .and. present(qvo)) then
      tvo_prime = 0.0
      tvo = to  * (1.0 + 0.608 * qvo)
   else if (present(to) .and. .not.present(qvo)) then
      tvo_prime = 0.0
      tvo = to
   else
      tvo_prime = tvsm_prime
      tvo = tvsm
   end if

   tv_prime = 0.5 * (tvsm_prime + tvo_prime)
   tv  = 0.5 * (tvsm + tvo)

   ! 2. HEIGHT DifFERENCE BEWTEEN MODEL SURFACE AND OBSERVATIONS
   ! ------------------------------------------------------------

   dz = hsm - ho
   arg0 = dz * g / gasr

   ! 3.  EXTRAPOLATE PRESSURE OBS TO MODEL SURFACE
   ! ---------------------------------------------

   arg_prime = - arg0 * tv_prime / (tv * tv)
   arg = arg0    / tv 

   psfcm_prime = exp(arg) *(psm_prime + psm * arg_prime)

   if (trace_use) call da_trace_exit("da_sfc_pre_lin")

end subroutine da_sfc_pre_lin


subroutine da_sfc_pre_adj (psfcm_prime, psm_prime, tsm_prime, qsm_prime, &
   psm, tsm, qsm, hsm, ho, to, qvo)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   !---------------------------------------------------------------------------
   !
   ! Correct pressure between two levels. 
   !
   ! Reference: make use of the hydrosatic equation:
   !
   !  P2 = P1 * exp [-G/R * (z2-z1) / (tv1 + tv2)/2)
   !
   ! Where:
   !  z1  = height at level 1
   !  z1  = height at level 2
   !  tv1 = temperature at level 1
   !  tv2 = temperature at level 2
   !  P1  = Pressure at level 1
   !  P2  = Pressure at level 2
   !---------------------------------------------------------------------------

   implicit none

   ! Perturbation:
   real, intent (in)     :: psfcm_prime          ! model pressure at ho
   real, intent (inout)  :: psm_prime, tsm_prime, qsm_prime            ! model surface p, t, q 
   ! Basic state:
   real, intent (in)     :: psm, tsm, qsm        ! model pressure at ho and
                                                 ! model surface p, t, q 
   ! Constant variables:
   real, intent (in)           :: hsm, ho
   real, intent (in), optional :: to, qvo
   ! working array:
   real                 :: tvo, tvsm, tv, dz, arg0
   real                 :: tvsm_prime, tvo_prime, tv_prime, arg, arg_prime

   real, parameter      :: GASR =  gas_constant
   real, parameter      :: G = gravity

   if (trace_use) call da_trace_entry("da_sfc_pre_adj")

   !---------------------------------------------------------------------------
   ! 1.0 Basic state
   ! --------------------------------------------------------------------------

   ! 1.1  MODEL AND OBSERVATION VIRTUAL TEMPERATURE
   ! ---------------------------------------------

   tvsm = tsm  * (1.0 + 0.608 * qsm)

   if (present(to) .and. present(qvo)) then
      tvo = to  * (1.0 + 0.608 * qvo)
   else if (present(to) .and. .not.present(qvo)) then
      tvo = to
   else
      tvo = tvsm
   end if

   ! 1.2  Mean virtual temperature
   ! ----------------------------

   tv  = 0.5 * (tvsm + tvo)

   ! 1.3  Compute (g/RTv) * dZ
   ! --------------------------

   dz = hsm - ho
   arg0 = dz * g / gasr     
   arg =  arg0    / tv

   ! ---------------------------------------------------------------------------|
   ! 2.0 Adjoint
   ! ---------------------------------------------------------------------------|

   ! 2.1  psfcm_prime ==> psm_prime, arg_prime
   ! -----------------------------------------

   arg_prime = exp(arg) * psm * psfcm_prime
   psm_prime = exp(arg) * psfcm_prime + psm_prime

   ! 2.2 arg_prim ==> tv_prime
   ! -------------------------

   tv_prime = - arg0 * arg_prime / (tv * tv)

   ! 2.3 tv_prime ==> tvsm_prime, tvo_prime
   ! --------------------------------------

   tvsm_prime = 0.5 * tv_prime
   tvo_prime  = 0.5 * tv_prime

   ! 2.4 tvo_prime ==> tsm_prime
   ! ---------------------------

   if (present(to) .and. present(qvo)) then
      tvo_prime = 0.0
   else if (present(to) .and. .not.present(qvo)) then
      tvo_prime = 0.0
   else
      tvsm_prime = tvo_prime + tvsm_prime
   end if

   ! 2.5 tvsm_prime ==>  tsm_prime, qsm_prime
   ! ----------------------------------------

   tsm_prime = tvsm_prime * (1.0 + 0.608 * qsm) + tsm_prime
   qsm_prime = tvsm_prime * tsm * 0.608 + qsm_prime

   if (trace_use) call da_trace_exit("da_sfc_pre_adj")

end subroutine da_sfc_pre_adj


subroutine da_moist_phys_adj(grid)

   !---------------------------------------------------------------------------
   !  Purpose: Partition of the hydrometeors via the moist egrid%xplicit scheme.
   !           A warm rain process is used in this subroutine. 
   !           This is the adjoint code of the scheme.
   !
   !  Method: The warm rain process is according to Hsie and Anthes (1984)
   !          and Dudhia (1989)
   !
   !  Assumptions: 1) Model level stored top down.
   !---------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)               :: grid

   real, dimension(ims:ime,jms:jme,kms:kme) :: T_OLD,T_NEW
   real, dimension(ims:ime,jms:jme,kms:kme) :: Q_OLD,Q_NEW
   real, dimension(ims:ime,jms:jme,kms:kme) :: QCW_OLD,QCW_NEW
   real, dimension(ims:ime,jms:jme,kms:kme) :: QRN_OLD,QRN_NEW

   real, dimension(kms:kme)                 :: EES, QVSS
   real, dimension(kms:kme)                 :: EES9, QVSS9

   real, dimension(its:ite,jts:jte,kms:kme) :: DT
   real, dimension(kms:kme)                   :: QVT,QCT,QRT,TTT
   real, dimension(kms:kme)                   :: QVT9,QCT9,QRT9,TTT9
   real, dimension(kms:kme) :: SCR2,SCR3,SCR4,SCR5,SCR6
   real, dimension(kms:kme) :: DUM31
   real, dimension(kms:kme) :: PRA,PRC,PRD,PRE
   real, dimension(kms:kme) :: SCR31,SCR42,SCR71
   real, dimension(kms:kme) :: DUM112,DUM113,DUM211,DUM411
   real, dimension(kms:kme) :: PRA2,PRC2

   real, dimension(kms:kme) :: SCR29,SCR39,SCR49,SCR59,SCR69
   real, dimension(kms:kme) :: DUM319
   real, dimension(kms:kme) :: PRA9,PRC9,PRD9,PRE9
   real, dimension(kms:kme) :: SCR319,SCR429,SCR719
   real, dimension(kms:kme) :: DUM1129,DUM1139,DUM2119,DUM4119
   real, dimension(kms:kme) :: TMP

   real, parameter :: QCTH  = 0.5E-3 
   real, parameter :: QRTH  = 1.0e-6  
   real, parameter :: alpha = 0.001
   real, parameter :: beta  = 0.0486
   real, parameter :: gamma = 0.002

   integer :: i, j, k
   real    :: tmp1, tmp2
   real    :: qrth1

   if (trace_use) call da_trace_entry("da_moist_phys_adj")

   qrth1 = (QRTH*1.0e3)**0.875

   QRN_NEW(its:ite,jts:jte,kts:kte) = grid%xa % qrn (its:ite,jts:jte,kts:kte) 
   QRN_OLD(its:ite,jts:jte,kts:kte) = grid%xa % qrn (its:ite,jts:jte,kts:kte)
   QCW_NEW(its:ite,jts:jte,kts:kte) = grid%xa % qcw (its:ite,jts:jte,kts:kte)
   QCW_OLD(its:ite,jts:jte,kts:kte) = grid%xa % qcw (its:ite,jts:jte,kts:kte)
   Q_NEW(its:ite,jts:jte,kts:kte) = grid%xa % q (its:ite,jts:jte,kts:kte)
   Q_OLD(its:ite,jts:jte,kts:kte) = grid%xa % q (its:ite,jts:jte,kts:kte)
   T_NEW(its:ite,jts:jte,kts:kte) = grid%xa % t (its:ite,jts:jte,kts:kte)
   T_OLD(its:ite,jts:jte,kts:kte) = grid%xa % t (its:ite,jts:jte,kts:kte)

   call da_filter_adj(grid,t_new)
   call da_filter_adj(grid,q_new)
   call da_filter_adj(grid,qcw_new)
   call da_filter_adj(grid,qrn_new)

   grid%xa % qrn (its:ite,jts:jte,kts:kte) = QRN_NEW(its:ite,jts:jte,kts:kte)
   QRN_OLD(its:ite,jts:jte,kts:kte) = QRN_OLD(its:ite,jts:jte,kts:kte) - QRN_NEW(its:ite,jts:jte,kts:kte)
   grid%xa % qcw (its:ite,jts:jte,kts:kte) = QCW_NEW(its:ite,jts:jte,kts:kte)
   QCW_OLD(its:ite,jts:jte,kts:kte) = QCW_OLD(its:ite,jts:jte,kts:kte) - QCW_NEW(its:ite,jts:jte,kts:kte)
   grid%xa % q (its:ite,jts:jte,kts:kte) = Q_NEW(its:ite,jts:jte,kts:kte)
   Q_OLD(its:ite,jts:jte,kts:kte) = Q_OLD(its:ite,jts:jte,kts:kte) - Q_NEW(its:ite,jts:jte,kts:kte)
   grid%xa % t (its:ite,jts:jte,kts:kte) = T_NEW(its:ite,jts:jte,kts:kte)
   T_OLD(its:ite,jts:jte,kts:kte) = T_OLD(its:ite,jts:jte,kts:kte) - T_NEW(its:ite,jts:jte,kts:kte)

   DT(its:ite,jts:jte,kts:kte) = grid%xb%delt(its:ite,jts:jte,kts:kte)

   do j=jts,jte
      do i=its,ite
         do K=kts,kte

            if (DT(i,j,k) <= 0.0) cycle

            if( grid%xb%t(I,J,K) > TO )then
               EES(K)=SVP1*EXP(SVP2*(grid%xb%t(I,J,K)-SVPT0)/(grid%xb%t(I,J,K)-SVP3))
            else
               EES(K)=.611*EXP(22.514-6.15E3/grid%xb%t(I,J,K))
            end if

            QVSS(K)=622.0*EES(K)/(grid%xb%p(I,J,K)-EES(K))

            SCR4(K)=grid%xb%q(I,J,K)/QVSS(K)

            if(grid%xb%qcw(I,J,K) > 0.0) then
               SCR2(K)=grid%xb%qcw(I,J,K)
            else
               SCR2(K)=0.0
            end if
            if(grid%xb%qrn(I,J,K) > 1.0e-25) then
               SCR3(K)=grid%xb%qrn(I,J,K)
            else
               SCR3(K)=1.0e-25
            end if
            SCR5(K)=grid%xb%q(I,J,K)/SCR4(K)

            SCR6(K)=grid%xb%p(I,J,K)/(gas_constant*grid%xb%t(I,J,K))

            DUM31(K)=3.1484E6-XLV1*grid%xb%t(I,J,K)

            ! autoc

            if ( SCR2(k) >= QCTH ) then
               PRC(k) = alpha * ( SCR2(k) - QCTH )
            else
               PRC(k) = 0.0
            end if

            PRC2(K)=PRC(K)

            ! accre

            if ( SCR2(k) > 0.0 .and. SCR3(k) > QRTH ) then
               PRA(k) = gamma * SCR2(k) * (SCR3(k)*1.0e3)**0.875
            else if ( SCR2(k) > 0.0 .and. SCR3(k) <= QRTH ) then
               PRA(k) = gamma * SCR2(k) * (QRTH*1.0e3)**0.875 
            else
               PRA(k) = 0.0
            end if

            PRA2(K)=PRA(K)

            ! evapo

            if (scr3(k) > qrth .and. grid%xb%q(I,J,k) < scr5(k)) then
               pre(k) = beta * (grid%xb%q(I,J,k)-scr5(k) ) * (scr6(k)*scr3(k))**0.65
            else if (scr3(k) <= qrth .and. scr3(k) > 0.0 .and. grid%xb%q(I,J,k) < scr5(k)) then
               pre(k) = beta * (grid%xb%q(I,J,k)-scr5(k)) * (scr6(k)*qrth)**0.65
            else
               pre(k) = 0.0
            end if

            if ( pre(k) < -scr3(k)/dt(i,j,k) ) then
               pre(k) = -scr3(k) / dt(i,j,k)
            end if

            !  Readjust

            DUM112(K)=(PRC(k)+PRA(k))*DT(i,j,k)
            if (DUM112(K) > SCR2(k)) then
               PRC(k)=SCR2(K)*PRC(K)/DUM112(K)
               PRA(k)=SCR2(K)*PRA(K)/DUM112(K)
            end if
            QVT(K)=-PRE(K)
            QCT(K)=-PRC(K)-PRA(K)
            QRT(K)=PRC(K)+PRA(K)+PRE(K)
            if(grid%xb%t(I,J,K).GT.TO)then
               DUM411(K)=DUM31(K)
            else
               DUM411(K)=XLS
            end if
            PRD(K)=cp*(1.0+0.887*grid%xb%q(I,J,K))
            TTT(K)=-DUM411(K)*QVT(K)/PRD(K)

            DUM113(K)=grid%xb%q(I,J,K)+DT(i,j,k)*QVT(K)
            if(DUM113(K) > 1.0E-12 ) then
               SCR42(K)=DUM113(K)
            else
               SCR42(K)=1.0E-12
            end if
            DUM211(K)=grid%xb%qcw(I,J,K)+QCT(K)*DT(i,j,k)
            if(DUM211(K) > 0.0) then
               SCR31(K)=DUM211(K)
            else
               SCR31(K)=0.0
            end if
            SCR71(K)=grid%xb%t(I,J,K)+TTT(K)*DT(i,j,k)
         end do

         call da_condens_adj(DT(i,j,:),SCR31,SCR42,SCR71,DUM31,PRD,       &
                             QVT,QCT,QRT,TTT,                      &
                             grid%xb%p(I,J,:),grid%xb%t(I,J,:),grid%xb%q(I,J,:),  &
                             grid%xb%qcw(I,J,:),grid%xb%qrn(I,J,:),          &
                             SCR319,SCR429,SCR719,DUM319,PRD9,     &
                             QVT9,QCT9,QRT9,TTT9,                  &
                             grid%xa%p(I,J,:),grid%xa%t(I,J,:),grid%xa%q(I,J,:),  &
                             grid%xa%qcw(I,J,:),grid%xa%qrn(I,J,:),kts,kte)

         do K=kts, kte
            if (DT(i,j,k) <= 0.0) cycle

            !  Readjust

            grid%xa%t(I,J,K)=grid%xa%t(I,J,K)+SCR719(K)
            TTT9(K)=TTT9(K)+DT(i,j,k)*SCR719(K)
            DUM2119(K)=0.0
            if(DUM211(K) > 0.0) then
               DUM2119(K)=SCR319(K)
            end if
            grid%xa%qcw(I,J,K)=grid%xa%qcw(I,J,K)+DUM2119(K)
            QCT9(K)=QCT9(K)+DT(i,j,k)*DUM2119(K)
            DUM1139(K)=0.0
            if(DUM113(K) > 1.0e-12 ) then
               DUM1139(K)=SCR429(K)
            end if
            grid%xa%q(I,J,K)=grid%xa%q(I,J,K)+DUM1139(K)
            QVT9(K)=QVT9(K)+DT(i,j,k)*DUM1139(K)
            PRD9(K)=PRD9(K)+DUM411(K)*QVT(K)/(PRD(K)*PRD(K))*TTT9(K)
            QVT9(K)=QVT9(K)-DUM411(K)/PRD(K)*TTT9(K)
            DUM4119(K)=-QVT(K)/PRD(K)*TTT9(K)
            grid%xa%q(I,J,K)=grid%xa%q(I,J,K)+cp*0.887*PRD9(K)

            if(grid%xb%t(I,J,K).GT.TO)then
               DUM319(K)=DUM319(K)+DUM4119(K)
            end if
            PRC9(K)=QRT9(K)
            PRA9(K)=QRT9(K)
            PRE9(K)=QRT9(K)
            PRC9(K)=PRC9(K)-QCT9(K)
            PRA9(K)=PRA9(K)-QCT9(K)
            PRE9(K)=PRE9(K)-QVT9(K)

            if (DUM112(K) > SCR2(k)) then
               DUM1129(K)=-SCR2(K)*PRA2(K)/(DUM112(K)*DUM112(K))*PRA9(K)
               SCR29(K)=PRA2(K)/DUM112(K)*PRA9(K)
               PRA9(K)=PRA9(K)*SCR2(K)/DUM112(K)
               DUM1129(K)=DUM1129(K)-SCR2(K)*PRC2(K)/(DUM112(K)*DUM112(K))*PRC9(K)
               SCR29(K)=SCR29(K)+PRC2(K)/DUM112(K)*PRC9(K)
               PRC9(K)=PRC9(K)*SCR2(K)/DUM112(K)
            else
               SCR29(K)=0.0
               DUM1129(K)=0.0        
            end if
            PRC9(K)=PRC9(K)+DT(i,j,k)*DUM1129(K)
            PRA9(K)=PRA9(K)+DT(i,j,k)*DUM1129(K)

            ! evapo

            if ( SCR3(K) > QRTH .and. grid%xb%q(I,J,k) < SCR5(k) ) then
               PRE(k) = beta * ( grid%xb%q(I,J,k)-SCR5(K) ) * ( SCR6(k)*SCR3(K) )**0.65
            else if ( SCR3(K) <= QRTH .and. SCR3(k) > 0.0 .and. grid%xb%q(I,J,k) < SCR5(k) ) then
               PRE(k) = beta * ( grid%xb%q(I,J,k)-SCR5(K) ) * ( SCR6(k)*QRTH )**0.65
            else
               PRE(k) = 0.0
            end if

            SCR39(k) = 0.0
            if ( PRE(k) < -SCR3(k)/DT(i,j,k) ) then
               SCR39(k) = -PRE9(k) / DT(i,j,k)
               PRE9(k)  = 0.0
            end if

            SCR59(k) = 0.0
            SCR69(k) = 0.0
            if ( SCR3(k) > QRTH .and. grid%xb%q(I,J,k) < SCR5(k) ) then
               TMP1  = beta * ( grid%xb%q(I,J,k)-SCR5(k) ) * 0.65 * ( SCR6(k)*SCR3(k) )**(-0.35)
               TMP2 = beta * ( SCR6(k)*SCR3(k) )**0.65

               grid%xa%q(I,J,k) = grid%xa%q(I,J,k) + TMP2 * PRE9(k)
               SCR59(k) = -TMP2 * PRE9(k)
               SCR39(k) = SCR39(k) + TMP1 * SCR6(k) * PRE9(k)
               SCR69(k) = TMP1 * SCR3(k) * PRE9(k)
            else if (SCR3(k) <= QRTH .and. SCR3(k) > 0.0 .and. grid%xb%q(I,J,k) < SCR5(k) ) then
               TMP1  = beta * ( grid%xb%q(I,J,k)-SCR5(k) ) * 0.65 * ( SCR6(k)*QRTH )**(-0.35)
               TMP2 = beta * ( SCR6(k)*QRTH )**0.65

               grid%xa%q(I,J,k) = grid%xa%q(I,J,k) + TMP2 * PRE9(k)
               SCR59(k) = -TMP2 * PRE9(k)
               SCR69(k) = TMP1 * QRTH * PRE9(k)
            end if
            ! accre

            if ( SCR2(k) > 0.0 .and. SCR3(k) > QRTH ) then
               SCR39(K) = SCR39(K) + gamma * 0.875 * SCR2(k) * (SCR3(K)*1.0e3)**(-0.125 ) * 1.0e3 * PRA9(K)
               SCR29(K) = SCR29(K) + gamma * (SCR3(K)*1.0e3)**0.875 * PRA9(K)
            else if (SCR2(k) > 0.0 .and. SCR3(k) <= QRTH ) then
               SCR29(k) = SCR29(k) + gamma * (QRTH1 * PRA9(k))
            end if
      
            ! autoc

            if ( scr2(k) >= qcth ) then
               scr29(k) = scr29(k) + alpha * prc9(k)
            end if

            grid%xa%t(I,J,K)=grid%xa%t(I,J,K)-XLV1*DUM319(K)

            grid%xa%p(I,J,K)=grid%xa%p(I,J,K)+SCR69(K)/(gas_constant*grid%xb%t(I,J,K))
            grid%xa%t(I,J,K)=grid%xa%t(I,J,K)-grid%xb%p(I,J,K)/  &
                        (gas_constant*grid%xb%t(I,J,K)**2)*SCR69(K)
            grid%xa%q(I,J,K)=grid%xa%q(I,J,K)+SCR59(K)/SCR4(K)
            SCR49(K)=-grid%xb%q(I,J,K)/SCR4(K)**2*SCR59(K)

            if(grid%xb%qrn(I,J,K) > 1.0e-25) then
               grid%xa%qrn(I,J,K)=grid%xa%qrn(I,J,K)+SCR39(K)
            end if
            if(grid%xb%qcw(I,J,K) > 0.0) then
               grid%xa%qcw(I,J,K)=grid%xa%qcw(I,J,K)+SCR29(K)
            end if

            grid%xa%q(I,J,K)=grid%xa%q(I,J,K)+SCR49(K)/QVSS(K)
            QVSS9(K)=-grid%xb%q(I,J,K)/QVSS(K)**2*SCR49(K)
            TMP(K)=622.0/((grid%xb%p(I,J,K)-EES(K))**2)
            EES9(K)=TMP(K)*grid%xb%p(I,J,K)*QVSS9(K)
            grid%xa%p(I,J,K)=grid%xa%p(I,J,K)-TMP(K)*EES(K)*QVSS9(K)
            if( grid%xb%t(I,J,K) > TO )then
               grid%xa%t(I,J,K)=grid%xa%t(I,J,K)+EES(K)*SVP2*(SVPT0-SVP3)/ ((grid%xb%t(I,J,K)-SVP3)*(grid%xb%t(I,J,K)-SVP3))*EES9(K)
            else
               grid%xa%t(I,J,K)=grid%xa%t(I,J,K)+EES(K)*6.15E3/(grid%xb%t(I,J,K)* grid%xb%t(I,J,K))*EES9(K)
            end if

         end do
      end do
   end do

   grid%xa % qt (its:ite,jts:jte,kts:kte) = grid%xa % qt (its:ite,jts:jte,kts:kte) + grid%xa % q(its:ite,jts:jte,kts:kte)
   grid%xa % qcw(its:ite,jts:jte,kts:kte) = grid%xa % qcw(its:ite,jts:jte,kts:kte) - grid%xa % q(its:ite,jts:jte,kts:kte)
   grid%xa % qrn(its:ite,jts:jte,kts:kte) = grid%xa % qrn(its:ite,jts:jte,kts:kte) - grid%xa % q(its:ite,jts:jte,kts:kte)

   grid%xa % qrn (its:ite,jts:jte,kts:kte) = grid%xa % qrn (its:ite,jts:jte,kts:kte) + QRN_OLD(its:ite,jts:jte,kts:kte)
   grid%xa % qcw (its:ite,jts:jte,kts:kte) = grid%xa % qcw (its:ite,jts:jte,kts:kte) + QCW_OLD(its:ite,jts:jte,kts:kte)
   grid%xa % q   (its:ite,jts:jte,kts:kte) = grid%xa % q (its:ite,jts:jte,kts:kte)   + Q_OLD(its:ite,jts:jte,kts:kte)
   grid%xa % t   (its:ite,jts:jte,kts:kte) = grid%xa % t (its:ite,jts:jte,kts:kte)   + T_OLD(its:ite,jts:jte,kts:kte)

   if (trace_use) call da_trace_exit("da_moist_phys_adj")

end subroutine da_moist_phys_adj


subroutine da_moist_phys_lin(grid)

   !---------------------------------------------------------------------------
   !  Purpose: Partition of the hydrometeors via the moist explicit scheme.
   !           A warm rain process is used in this subroutine.
   !           This is the tangent linear code of the scheme.
   !
   !  Method: The warm rain process is according to Hsie and Anthes (1984)
   !          and Dudhia (1989)
   !
   !  Assumptions: 1) Model level stored top down.
   !---------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)               :: grid

   real :: T_OLD(ims:ime,jms:jme,kms:kme),T_NEW(ims:ime,jms:jme,kms:kme)
   real :: Q_OLD(ims:ime,jms:jme,kms:kme),Q_NEW(ims:ime,jms:jme,kms:kme)
   real :: QCW_OLD(ims:ime,jms:jme,kms:kme),QCW_NEW(ims:ime,jms:jme,kms:kme)
   real :: QRN_OLD(ims:ime,jms:jme,kms:kme),QRN_NEW(ims:ime,jms:jme,kms:kme)

   real    :: EES(kms:kme)
   real    :: QVSS(kms:kme)
   real    :: EES9(kms:kme)
   real    :: QVSS9(kms:kme)
   real    :: DT(its:ite,jts:jte,kms:kme)
   real    :: QVT(kms:kme)
   real    :: QCT(kms:kme)
   real    :: QRT(kms:kme)
   real    :: TTT(kms:kme)
   real    :: QVT9(kms:kme)
   real    :: QCT9(kms:kme)
   real    :: QRT9(kms:kme)
   real    :: TTT9(kms:kme)
   real    :: SCR2(kms:kme)
   real    :: SCR3(kms:kme)
   real    :: SCR4(kms:kme)
   real    :: SCR5(kms:kme)
   real    :: SCR6(kms:kme)
   real    :: DUM31(kms:kme)
   real    :: PRA(kms:kme)
   real    :: PRC(kms:kme)
   real    :: PRD(kms:kme)
   real    :: PRE(kms:kme)
   real    :: SCR31(kms:kme)
   real    :: SCR42(kms:kme)
   real    :: SCR71(kms:kme)
   real    :: DUM112(kms:kme)
   real    :: DUM113(kms:kme)
   real    :: DUM211(kms:kme)
   real    :: DUM411(kms:kme)
   real    :: SCR29(kms:kme)
   real    :: SCR39(kms:kme)
   real    :: SCR49(kms:kme)
   real    :: SCR59(kms:kme)
   real    :: SCR69(kms:kme)
   real    :: DUM319(kms:kme)
   real    :: PRA9(kms:kme)
   real    :: PRC9(kms:kme)
   real    :: PRD9(kms:kme)
   real    :: PRE9(kms:kme)
   real    :: SCR319(kms:kme)
   real    :: SCR429(kms:kme)
   real    :: SCR719(kms:kme)
   real    :: DUM1129(kms:kme)
   real    :: DUM4119(kms:kme)
   real    :: TMP(kms:kme)

   integer, parameter :: qcth = 0.5e-3
   integer, parameter :: qrth = 1.0e-6
   integer, parameter :: alpha = 0.001
   integer, parameter :: beta  = 0.0486
   integer, parameter :: gamma = 0.002

   integer :: i, j, k

   real    :: qrth1

   if (trace_use) call da_trace_entry("da_moist_phys_lin")

   qrth1 = (QRTH*1.0e3)**0.875

   T_OLD  (its:ite,jts:jte,kts:kte) = grid%xa%t  (its:ite,jts:jte,kts:kte)
   Q_OLD  (its:ite,jts:jte,kts:kte) = grid%xa%q  (its:ite,jts:jte,kts:kte)
   QCW_OLD(its:ite,jts:jte,kts:kte) = grid%xa%qcw(its:ite,jts:jte,kts:kte)
   QRN_OLD(its:ite,jts:jte,kts:kte) = grid%xa%qrn(its:ite,jts:jte,kts:kte)

   !  Preparation

   grid%xa%q(its:ite,jts:jte,kts:kte) =grid%xa%qt(its:ite,jts:jte,kts:kte) - grid%xa%qcw(its:ite,jts:jte,kts:kte) - grid%xa %qrn(its:ite,jts:jte,kts:kte)
   DT(its:ite,jts:jte,kts:kte) = grid%xb%delt(its:ite,jts:jte,kts:kte)

   do j=jts,jte
      do i=its,ite
         do K=kts,kte

            if (dt(i,j,k) <= 0.0) cycle

            if ( grid%xb%t(I,J,K) > TO )then
               EES(K)=SVP1*EXP(SVP2*(grid%xb%t(I,J,K)-SVPT0)/(grid%xb%t(I,J,K)-SVP3))
               EES9(K)=EES(K)*SVP2*(SVPT0-SVP3)/((grid%xb%t(I,J,K)-SVP3) * (grid%xb%t(I,J,K)-SVP3))*grid%xa%t(I,J,K)
            else
               EES(K)=.611*EXP(22.514-6.15E3/grid%xb%t(I,J,K))
               EES9(K)=EES(K)*6.15E3/(grid%xb%t(I,J,K)*grid%xb%t(I,J,K))*grid%xa%t(I,J,K)
            end if

            TMP(K)=622.0/((grid%xb%p(I,J,K)-EES(K))**2)
            QVSS9(K)=TMP(K)*grid%xb%p(I,J,K)*EES9(K) - TMP(K)*EES(K)*grid%xa%p(I,J,K)
            QVSS(K)=622.0*EES(K)/(grid%xb%p(I,J,K)-EES(K))

            SCR49(K)=grid%xa%q(I,J,K)/QVSS(K)-grid%xb%q(I,J,K)/QVSS(K)**2*QVSS9(K)
            SCR4(K)=grid%xb%q(I,J,K)/QVSS(K)

            if (grid%xb%qcw(I,J,K) > 0.0) then
               SCR29(K)=grid%xa%qcw(I,J,K)
               SCR2(K)=grid%xb%qcw(I,J,K)
            else
               SCR29(K)=0.0
               SCR2(K)=0.0
            end if
            if (grid%xb%qrn(I,J,K) > 1.0e-25) then
               SCR39(K)=grid%xa%qrn(I,J,K)
               SCR3(K)=grid%xb%qrn(I,J,K)
            else
               SCR39(K)=0.0
               SCR3(K)=1.0E-25
            end if
            SCR59(K)=grid%xa%q(I,J,K)/SCR4(K)-grid%xb%q(I,J,K)/SCR4(K)**2*SCR49(K)
            SCR5(K)=grid%xb%q(I,J,K)/SCR4(K)

            SCR69(K)=grid%xa%p(I,J,K)/(gas_constant*grid%xb%t(I,J,K))-grid%xb%p(I,J,K)/  &
                     (gas_constant*grid%xb%t(I,J,K)**2)*grid%xa%t(I,J,K)
            SCR6(K)=grid%xb%p(I,J,K)/(gas_constant*grid%xb%t(I,J,K))

            DUM319(K)=-XLV1*grid%xa%t(I,J,K) 
            DUM31(K)=3.1484E6-XLV1*grid%xb%t(I,J,K)
 
            ! Auto conversion

            if (scr2(k) >= qcth) then
               prc9(k) = alpha * scr29(k)
               prc(k) = alpha * (scr2(k) - qcth)
            else
               prc9(k) = 0.0
               prc(k) = 0.0
            end if

            ! Accretion

            if (SCR2(k) > 0.0 .and. SCR3(k) > QRTH ) then
               PRA9(K) = gamma * 0.875 * SCR2(k) * (SCR3(K)*1.0e3)**(-0.125) * 1.0e3 * SCR39(K)  &
                            + gamma * SCR29(k) * (SCR3(K)*1.0e3)**0.875
               PRA(k) = gamma * SCR2(k) * (SCR3(k)*1.0e3)**0.875
            else if (SCR2(k) > 0.0 .and. SCR3(k) <= QRTH ) then
               PRA9(K) = gamma * SCR29(k) * qrth1
               PRA(k) = gamma * SCR2(k) * qrth1
            else
               PRA9(K) = 0.0
               PRA(k) = 0.0
            end if

         end do

         call da_evapo_lin(DT(i,j,:),SCR3,SCR5,grid%xb%q(I,J,:),PRE,SCR6,  &
                           SCR39,SCR59,grid%xa%q(I,J,:),PRE9,SCR69, &
                           kts,kte,kms,kme)

         do K=kts, kte

            if (dt(i,j,k) <= 0.0) cycle

            !  Readjust

            DUM112(K)=(PRC(k)+PRA(k))*dt(i,j,k)
            if (DUM112(K) > SCR2(k)) then
               DUM1129(K)=(PRC9(k)+PRA9(k))*dt(i,j,k)
               PRC9(K)=SCR29(K)*PRC(K)/DUM112(K)  &
                      +PRC9(K)*SCR2(K)/DUM112(K)  &
                      -SCR2(K)*PRC(K)/(DUM112(K)*DUM112(K))*DUM1129(K)
               PRC(k)=SCR2(K)*PRC(K)/DUM112(K)
               PRA9(K)=SCR29(K)*PRA(K)/DUM112(K)  &
                      +PRA9(K)*SCR2(K)/DUM112(K)  &
                      -SCR2(K)*PRA(K)/(DUM112(K)*DUM112(K))*DUM1129(K)
               PRA(k)=SCR2(K)*PRA(K)/DUM112(K)
            end if
            QVT9(K)=-PRE9(K)
            QVT(K)=-PRE(K)
            QCT9(K)=-PRC9(K)-PRA9(K)
            QCT(K)=-PRC(K)-PRA(K)
            QRT9(K)=PRC9(K)+PRA9(K)+PRE9(K)
            QRT(K)=PRC(K)+PRA(K)+PRE(K)
            if (grid%xb%t(I,J,K).GT.TO)then
               DUM4119(K)=DUM319(K)
               DUM411(K)=DUM31(K)
            else
               DUM4119(K)=0.0
               DUM411(K)=XLS
            end if
            PRD9(K)=cp*0.887*grid%xa%q(I,J,K)
            PRD(K)=cp*(1.0+0.887*grid%xb%q(I,J,K))
            TTT9(K)=-DUM4119(K)*QVT(K)/PRD(K)  &
                   -QVT9(K)*DUM411(K)/PRD(K)  &
                   +DUM411(K)*QVT(K)/(PRD(K)*PRD(K))*PRD9(K)
            TTT(K)=-DUM411(K)*QVT(K)/PRD(K)

            DUM113(K)=grid%xb%q(I,J,K)+dt(i,j,k)*QVT(K)
            if (DUM113(K) > 1.0e-12 ) then
               SCR429(K)=grid%xa%q(I,J,K)+dt(i,j,k)*QVT9(K)
               SCR42(K)=DUM113(K)
            else
               SCR429(K)=0.0
               SCR42(K)=1.0e-12
            end if
            DUM211(K)=grid%xb%qcw(I,J,K)+QCT(K)*dt(i,j,k)
            if (DUM211(K) > 0.0) then
               SCR319(K)=grid%xa%qcw(I,J,K)+QCT9(K)*dt(i,j,k)
               SCR31(K)=DUM211(K)
            else
               SCR319(K)=0.0
               SCR31(K)=0.0
            end if
            SCR719(K)=grid%xa%t(I,J,K)+TTT9(K)*dt(i,j,k)
            SCR71(K)=grid%xb%t(I,J,K)+TTT(K)*dt(i,j,k)
         end do

         call da_condens_lin(DT(i,j,:),SCR31,SCR42,SCR71,DUM31,PRD,         &
                             QVT,QCT,QRT,TTT,                        &
                             grid%xb%p(I,J,:),grid%xb%t(I,J,:),grid%xb%q(I,J,:),    &
                             grid%xb%qcw(I,J,:),grid%xb%qrn(I,J,:),            &
                             SCR319,SCR429,SCR719,DUM319,PRD9,       &
                             QVT9,QCT9,QRT9,TTT9,                    &
                             grid%xa%p(I,J,:),grid%xa%t(I,J,:),grid%xa%q(I,J,:),    &
                             grid%xa%qcw(I,J,:),grid%xa%qrn(I,J,:),kts,kte)
      end do
   end do

   T_NEW  (its:ite,jts:jte,kds:kde) = grid%xa%t   (its:ite,jts:jte,kds:kde) - T_OLD  (its:ite,jts:jte,kds:kde)
   Q_NEW  (its:ite,jts:jte,kds:kde) = grid%xa%q   (its:ite,jts:jte,kds:kde) - Q_OLD  (its:ite,jts:jte,kds:kde)
   QCW_NEW(its:ite,jts:jte,kds:kde) = grid%xa%qcw (its:ite,jts:jte,kds:kde) - QCW_OLD(its:ite,jts:jte,kds:kde)
   QRN_NEW(its:ite,jts:jte,kds:kde) = grid%xa%qrn (its:ite,jts:jte,kds:kde) - QRN_OLD(its:ite,jts:jte,kds:kde)

   call da_filter(grid, t_new)
   call da_filter(grid, q_new)
   call da_filter(grid, qcw_new)
   call da_filter(grid, qrn_new)

   grid%xa%t   (its:ite,jts:jte,kds:kde) = T_NEW  (its:ite,jts:jte,kds:kde) + T_OLD  (its:ite,jts:jte,kds:kde)
   grid%xa%q   (its:ite,jts:jte,kds:kde) = Q_NEW  (its:ite,jts:jte,kds:kde) + Q_OLD  (its:ite,jts:jte,kds:kde)
   grid%xa%qcw (its:ite,jts:jte,kds:kde) = QCW_NEW(its:ite,jts:jte,kds:kde) + QCW_OLD(its:ite,jts:jte,kds:kde)
   grid%xa%qrn (its:ite,jts:jte,kds:kde) = QRN_NEW(its:ite,jts:jte,kds:kde) + QRN_OLD(its:ite,jts:jte,kds:kde)

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_XA_CLOUD.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_XA_CLOUD_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE

   if (trace_use) call da_trace_exit("da_moist_phys_lin")

end subroutine da_moist_phys_lin


subroutine da_condens_adj(DT,SCR31,SCR42,SCR71,DUM31,PRD,  &
                   QVT,QCT,QRT,TTT,P_B,T_B,QV_B,QCW_B,QRN_B,  &
                   SCR319,SCR429,SCR719,DUM319,PRD9,  &
                   QVT9,QCT9,QRT9,TTT9,P_A,T_A,QV_A,QCW_A,QRN_A,kts,kte)

   !-----------------------------------------------------------------------
   ! Purpose: Condensation
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)                     :: kts, kte
   real, dimension(kts:kte), intent(in)    :: DT,SCR31,SCR42,SCR71,PRD,DUM31
   real, dimension(kts:kte), intent(in)    :: P_B,T_B,QV_B,QCW_B,QRN_B
   real, dimension(kts:kte), intent(inout) :: SCR319,SCR429,SCR719,PRD9
   real, dimension(kts:kte), intent(inout) :: P_A,T_A,QV_A,QCW_A,QRN_A,DUM319

   real, dimension(kts:kte), intent(in)    :: QVT,QCT,QRT,TTT
   real, dimension(kts:kte), intent(inout) :: QRT9,QCT9,QVT9,TTT9


   real, dimension(kts:kte)  :: DUM2139
   real, dimension(kts:kte)  :: TMP,DUM114,DUM2129,SCR89,DUM212,DUM115
   real, dimension(kts:kte)  :: PRC5,PRC59,DUM1149,SCR61,SCR8,DUM213
   real, dimension(kts:kte)  :: SCR619
   integer                   :: k

   !  initilization

   do K=kts,kte
      DUM2129(K) = 0.0
      SCR89 (K) = 0.0
      PRC59 (K) = 0.0
   end do

   do K=kts, kte

   if (DT(k) <= 0.0) cycle

      DUM114(K)=1.0e3*SVP1*EXP(SVP2*(SCR71(K)-SVPT0)/(SCR71(K)-SVP3))

      if(SCR71(K) > TO) then
         DUM212(K)=DUM31(K)*DUM31(K)/(gas_constant_v*PRD(K))
      else
         DUM212(K)=XLS*DUM31(K)/(gas_constant_v*PRD(K))
      end if
      PRC5(K)=.622*DUM114(K)/(P_B(K)-DUM114(K))

      if(SCR42(K) < PRC5(K) .AND. SCR71(K) < TO) then
         SCR61(K)=0.0
      else
         SCR8(K)=(SCR42(K)-PRC5(K))/(1.0+DUM212(K)*PRC5(K)/  &
                 (SCR71(K)*SCR71(K)))

         DUM115(K)=SCR31(K)+SCR8(K)
         if (DUM115(K) >= 0.0)then
            SCR61(K)=SCR8(K)/DT(k)
         else
            SCR61(K)=-SCR31(K)/DT(k)
         end if
      end if
      if(SCR71(K) > TO)then
         DUM213(K)=DUM31(K)/PRD(K)
      else
         DUM213(K)=XLS/PRD(K)
      end if

      TTT9(K)=DT(K)*T_A(K)
      SCR619(K)=DT(K)*DUM213(K)*T_A(K)
      DUM2139(K)=DT(K)*SCR61(K)*T_A(K)
      if(QRN_B(K) < 1.0e-25) QRN_A(K)=0.0
      QRT9(K)=DT(K)*QRN_A(K)

      DUM319(K)=0.0
      if(SCR71(K) > TO)then
         DUM319(K)=DUM2139(K)/PRD(K)
         PRD9(K)=-DUM31(K)/(PRD(K)*PRD(K))*DUM2139(K)
      else
         PRD9(K)=-XLS/(PRD(K)*PRD(K))*DUM2139(K)
      end if
      if(QCW_B(K) < 1.0e-25) QCW_A(K)=0.0
      QCT9(K)=DT(K)*QCW_A(K)
      SCR619(K)=SCR619(K)+DT(K)*QCW_A(K)
      if(QV_B(K) < 1.0e-25) QV_A(K)=0.0
      QVT9(K)=DT(K)*QV_A(K)
      SCR619(K)=SCR619(K)-DT(K)*QV_A(K)

      SCR319(K)=0.0
      SCR429(K)=0.0
      SCR719(K)=0.0
      if(SCR42(K) >= PRC5(K) .OR. SCR71(K) >= TO) then
         if(DUM115(K) >= 0.0)then
            SCR89(K)=SCR89(K)+SCR619(K)/DT(k)
         else
            SCR319(K)=-SCR619(K)/DT(k)
         end if

         TMP(K)=1.0/(1.0+DUM212(K)*PRC5(K)/(SCR71(K)*SCR71(K)))
         SCR719(K)=TMP(K)*TMP(K)*2.0*DUM212(K)*PRC5(K)  &
                   *(SCR42(K)-PRC5(K))/(SCR71(K)*SCR71(K)*SCR71(K))*SCR89(K)
         DUM2129(K)=DUM2129(K)-TMP(K)*TMP(K)*(SCR42(K)-PRC5(K))*PRC5(K)/  &
                     (SCR71(K)*SCR71(K))*SCR89(K)
         PRC59(K)=PRC59(K)-TMP(K)*(1.0+(SCR42(K)-PRC5(K))*DUM212(K)/  &
                     (SCR71(K)*SCR71(K))*TMP(K))*SCR89(K)
         SCR429(K)=TMP(K)*SCR89(K)
      end if

      TMP(K)=.622/(P_B(K)-DUM114(K))**2
      DUM1149(K)=TMP(K)*P_B(K)*PRC59(K)
      P_A(K)=P_A(K)-TMP(K)*DUM114(K)*PRC59(K)
      if(SCR71(K) > TO) then
         PRD9(K)=PRD9(K)-DUM31(K)*DUM31(K)/   &
                 (gas_constant_v*PRD(K)*PRD(K))*DUM2129(K)
         DUM319(K)=DUM319(K)+2.0*DUM31(K)/(gas_constant_v*PRD(K))*DUM2129(K)
      else
         PRD9(K)=PRD9(K)-XLS*DUM31(K)/(gas_constant_v*PRD(K)*PRD(K))*DUM2129(K)
         DUM319(K)=DUM319(K)+XLS/(gas_constant_v*PRD(K))*DUM2129(K)
      end if
      DUM114(K)=1.0e3*SVP1*EXP(SVP2*(SCR71(K)-SVPT0)/(SCR71(K)-SVP3))
      SCR719(K)=SCR719(K)+DUM114(K)*SVP2*(SVPT0-SVP3)/  &
                (SCR71(K)-SVP3)**2*DUM1149(K)

   end do

end subroutine da_condens_adj
subroutine da_condens_lin(DT,SCR31,SCR42,SCR71,DUM31,PRD,  &
                   QVT,QCT,QRT,TTT,P_B,T_B,QV_B,QCW_B,QRN_B,  &
                   SCR319,SCR429,SCR719,DUM319,PRD9,  &
                   QVT9,QCT9,QRT9,TTT9,P_A,T_A,QV_A,QCW_A,QRN_A,kts,kte)

   !-----------------------------------------------------------------------
   ! Purpose: Condensation
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)       :: kts, kte
   ! real                      :: SVP1,SVP2,SVP3,SVPT0,TO,gas_constant_v,XLS
   real, dimension(kts:kte), intent(in)    :: DT,SCR31,SCR42,SCR71,PRD,DUM31
   real, dimension(kts:kte), intent(in)    :: P_B,T_B,QV_B,QCW_B,QRN_B
   real, dimension(kts:kte), intent(in)    :: SCR319,SCR429,SCR719,PRD9,DUM319
   real, dimension(kts:kte), intent(in)    :: P_A
   real, dimension(kts:kte), intent(inout) :: T_A,QV_A,QCW_A,QRN_A

   real, dimension(kts:kte), intent(in)  :: QVT,QCT,QRT,TTT
   real, dimension(kts:kte), intent(in)  :: QVT9,QCT9,QRT9,TTT9

   real, dimension(kts:kte)  :: TMP,PRC59,DUM2139,PRC5,DUM213,SCR619,SCR89,SCR8
   real, dimension(kts:kte)  :: SCR61
   real, dimension(kts:kte)  :: DUM114,DUM1149,DUM2129,DUM115,DUM212,DUM1159

   integer                   :: k

   do K=kts, kte

      if (DT(k) <= 0.0) cycle

      DUM114(K)=1.0E3*SVP1*EXP(SVP2*(SCR71(K)-SVPT0)/(SCR71(K)-SVP3))
      DUM1149(k)=DUM114(K)*SVP2*(SVPT0-SVP3)/(SCR71(K)-SVP3)**2*SCR719(K)

      if(SCR71(K) > TO) then
         DUM2129(K)=2.0*DUM31(K)/(gas_constant_v*PRD(K))*DUM319(K)  &
                  -DUM31(K)*DUM31(K)/(gas_constant_v*PRD(K)*PRD(K))*PRD9(K)
         DUM212(K)=DUM31(K)*DUM31(K)/(gas_constant_v*PRD(K))
      else
         DUM2129(K)=XLS*DUM319(K)/(gas_constant_v*PRD(K))  &
                  -XLS*DUM31(K)/(gas_constant_v*PRD(K)*PRD(K))*PRD9(K)
         DUM212(K)=XLS*DUM31(K)/(gas_constant_v*PRD(K))
      end if
      TMP(K)=0.622/(P_B(K)-DUM114(K))**2
      PRC59(K)=TMP(K)*P_B(K)*DUM1149(K)-TMP(K)*DUM114(K)*P_A(K)
      PRC5(K)=0.622*DUM114(K)/(P_B(K)-DUM114(K))

      if(SCR42(K) < PRC5(K) .AND. SCR71(K) < TO) then
         SCR619(K)=0.0
         SCR61(K)=0.0
      else
         TMP(K)=1./(1.0+DUM212(K)*PRC5(K)/(SCR71(K)*SCR71(K)))
         SCR89(K)=TMP(K)*SCR429(K)  &
                 -TMP(K)*(1.0+(SCR42(K)-PRC5(K))*DUM212(K)/  &
                     (SCR71(K)*SCR71(K))*TMP(K))*PRC59(K)  &
                 -TMP(K)*TMP(K)*(SCR42(K)-PRC5(K))*PRC5(K)/  &
                     (SCR71(K)*SCR71(K))*DUM2129(K)  &
         ! WHY?
         !error      -TMP(K)*TMP(K)*2.0*DUM212(K)*PRC5(K)  &
         !error      *(SCR42(K)-PRC5(K))/(SCR71(K)*SCR71(K))*SCR719(K)
                 +TMP(K)*TMP(K)*2.0*DUM212(K)*PRC5(K)  &
                      *(SCR42(K)-PRC5(K))/(SCR71(K)*SCR71(K)*SCR71(K))*SCR719(K)
         SCR8(K)=(SCR42(K)-PRC5(K))/(1.0+DUM212(K)*PRC5(K)/  &
                 (SCR71(K)*SCR71(K)))

         DUM1159(K)=SCR319(K)+SCR89(K)
         DUM115(K)=SCR31(K)+SCR8(K)
         if(DUM115(K) >= 0.0)then
            SCR619(K)=SCR89(K)/DT(k)
            SCR61(K)=SCR8(K)/DT(k)
         else
            SCR619(K)=-SCR319(K)/DT(k)
            SCR61(K)=-SCR31(K)/DT(k)
         end if
      end if
      QV_A(K)=QV_A(K)+(QVT9(K)-SCR619(K))*DT(K)
      if(QV_B(K) < 1.0E-25) QV_A(K)=0.0
      QCW_A(K)=QCW_A(K)+(QCT9(K)+SCR619(K))*DT(K)
      if(QCW_B(K) < 1.0E-25) QCW_A(K)=0.0
      if(SCR71(K) > TO)then
         DUM2139(K)=DUM319(K)/PRD(K)-DUM31(K)/(PRD(K)*PRD(K))*PRD9(K)
         DUM213(K)=DUM31(K)/PRD(K)
      else
         DUM2139(K)=-XLS/(PRD(K)*PRD(K))*PRD9(K)
         DUM213(K)=XLS/PRD(K)
      end if

      QRN_A(K)=QRN_A(K)+DT(K)*QRT9(K)
      if(QRN_B(K) < 1.0E-25) QRN_A(K)=0.0
      T_A(K)=T_A(K)+DT(K)*(TTT9(K)+SCR619(K)*DUM213(K)+SCR61(K)*DUM2139(K))

   end do

end subroutine da_condens_lin
subroutine da_evapo_lin(dt,scr3,scr5,qv_b,pre,scr6, scr39,scr59,qv_a,pre9,scr69, kts,kte,kms,kme)

   !-----------------------------------------------------------------------
   ! Purpose: Rainwater evaporation
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)  :: kts, kte, kms, kme
   real,    intent(in)  :: dt(kms:kme)
   real,    intent(in)  :: scr3(kms:kme)
   real,    intent(in)  :: scr5(kms:kme)
   real,    intent(in)  :: scr6(kms:kme)
   real,    intent(in)  :: qv_b(kms:kme)
   real,    intent(out) :: pre(kms:kme)
   real,    intent(out) :: pre9(kms:kme)
   real,    intent(in)  :: scr39(kms:kme)
   real,    intent(in)  :: scr59(kms:kme)
   real,    intent(in)  :: scr69(kms:kme)
   real,    intent(in)  :: qv_a(kms:kme)

   integer :: k
   real    :: beta, qrth
   real    :: tmp, tmp2

   if (trace_use) call da_trace_entry("da_evapo_lin")

   qrth = 1.0e-6
   beta = 0.0486   ! original

   do k = kts, kte
      if (dt(k) <= 0.0) cycle

      if ( scr3(k) > qrth .and. qv_b(k) < scr5(k) ) then
         tmp  = beta * ( qv_b(k)-scr5(k) )* 0.65 * ( scr6(k)*scr3(k) )**(-0.35)
         tmp2 = beta * ( scr6(k)*scr3(k) )**0.65
         pre9(k) = tmp * ( scr69(k)*scr3(k)+scr6(k)*scr39(k) ) + &
                   tmp2 * ( qv_a(k)-scr59(k) )
         pre(k)  = beta * ( qv_b(k)-scr5(k) ) * ( scr6(k)*scr3(k) )**0.65
      else if ( scr3(k) <= qrth .and. scr3(k) > 0.0 .and. qv_b(k) < scr5(k) ) then
         tmp  = beta * ( qv_b(k)-scr5(k) ) * 0.65 * ( scr6(k)*qrth )**(-0.35)
         tmp2 = beta * ( scr6(k)*qrth )**0.65
         pre9(k) = tmp * ( scr69(k)*qrth ) + tmp2 * ( qv_a(k)-scr59(k) )
         pre(k)  = beta * ( qv_b(k)-scr5(k) ) * ( scr6(k)*qrth )**0.65
      else
         pre9(k) = 0.0
         pre(k) = 0.0
      end if

      if ( pre(k) < -scr3(k)/dt(k) ) then
         pre9(k) = -scr39(k) / dt(k)
         pre(k)  = -scr3(k) / dt(k)
      end if
   end do

   if (trace_use) call da_trace_exit("da_evapo_lin")

end subroutine da_evapo_lin

subroutine da_filter(grid, var)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(domain), intent(inout) :: grid
   real,         intent(inout) :: var(ims:ime,jms:jme,kms:kme)

   integer :: i, j, k
   real    :: w(3)

   data  w/0.25,0.5,0.25/

   if (trace_use) call da_trace_entry("da_filter")

   ! Copy var for transpose.

   grid%xp%v1z(its:ite,jts:jte,kts:kte) = &
      var(its:ite,jts:jte,kts:kte)

   ! Apply (i',j',k -> i,j',k') transpose (v1z -> v1x).

   call da_transpose_z2x (grid)

   ! Perform x-direction filter:

   do k = grid%xp%ktsx, grid%xp%ktex
      do j=grid%xp%jtsx, grid%xp%jtex
         ! Forward 
         do i=ids+1, ide-1
            grid%xp%v1x(i,j,k) = w(1)*grid%xp%v1x(i-1,j,k) + w(2)*grid%xp%v1x(i,j,k) + &
               w(3)*grid%xp%v1x(i+1,j,k)
         end do

         ! Backward
         do i=ide-1,ids+1,-1
            grid%xp%v1x(i,j,k) = w(1)*grid%xp%v1x(i-1,j,k) + w(2)*grid%xp%v1x(i,j,k) + &
               w(3)*grid%xp%v1x(i+1,j,k)
         end do
      end do
   end do


   ! Apply (i,j',k' -> i',j,k') transpose (v1x -> v1y).

   call da_transpose_x2y (grid)

   ! Perform y-direction filter:

   do k=grid%xp%ktsy, grid%xp%ktey
      do i=grid%xp%itsy, grid%xp%itey
         ! Forward
         do j=jds+1, jde-1
            grid%xp%v1y(i,j,k) = w(1)*grid%xp%v1y(i,j-1,k) + w(2)*grid%xp%v1y(i,j,k) + &
               w(3)*grid%xp%v1y(i,j+1,k)
         end do

         ! Backward
         do j=jde-1,jds+1,-1
            grid%xp%v1y(i,j,k) = w(1)*grid%xp%v1y(i,j-1,k) + w(2)*grid%xp%v1y(i,j,k) + &
               w(3)*grid%xp%v1y(i,j+1,k)
         end do
      end do
   end do

   ! Apply (i',j,k' -> i',j',k) transpose (v1y -> v1z).

   call da_transpose_y2z (grid)

   var(its:ite,jts:jte,kts:kte) = grid%xp%v1z(its:ite,jts:jte,kts:kte)

   if (trace_use) call da_trace_exit("da_filter")

end subroutine da_filter


subroutine da_filter_adj(grid, var)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(domain), intent(inout) :: grid
   real,         intent(inout) :: var(ims:ime,jms:jme,kms:kme)

   integer :: i, j, k

   real    :: w(3)

   data  w/0.25,0.5,0.25/

   if (trace_use) call da_trace_entry("da_filter_adj")

   ! Copy var for transpose.

   grid%xp%v1z(its:ite,jts:jte,kts:kte) = var(its:ite,jts:jte,kts:kte)

   ! Apply (i',j',k -> i',j,k') transpose (v1z -> v1y).
   call da_transpose_z2y (grid)

   ! Perform y-direction filter:
   do k=grid%xp%ktsy, grid%xp%ktey
      do i=grid%xp%itey, grid%xp%itsy, -1
         ! Backward
         do j=jds+1, jde-1
            grid%xp%v1y(i,j-1,k) = grid%xp%v1y(i,j-1,k) + w(1)*grid%xp%v1y(i,j,k)
            grid%xp%v1y(i,j+1,k) = grid%xp%v1y(i,j+1,k) + w(3)*grid%xp%v1y(i,j,k)
            grid%xp%v1y(i,j  ,k) = w(2)*grid%xp%v1y(i,j,k)
         end do

         ! Forward
         do j=jde-1,jds+1,-1
            grid%xp%v1y(i,j-1,k) = grid%xp%v1y(i,j-1,k) + w(1)*grid%xp%v1y(i,j,k)
            grid%xp%v1y(i,j+1,k) = grid%xp%v1y(i,j+1,k) + w(3)*grid%xp%v1y(i,j,k)
            grid%xp%v1y(i,j  ,k) = w(2)*grid%xp%v1y(i,j,k)
         end do
      end do
   end do

   ! Apply (i',j,k' -> i,j',k') transpose (v1y -> v1x).
   call da_transpose_y2x (grid)

   ! Perform x-direction filter:
   do k=grid%xp%ktsx, grid%xp%ktex
      do j=grid%xp%jtex, grid%xp%jtsx, -1
         ! Backward
         do i=ids+1, ide-1
            grid%xp%v1x(i-1,j,k) = grid%xp%v1x(i-1,j,k) + w(1)*grid%xp%v1x(i,j,k)
            grid%xp%v1x(i+1,j,k) = grid%xp%v1x(i+1,j,k) + w(3)*grid%xp%v1x(i,j,k)
            grid%xp%v1x(i,j  ,k) = w(2)*grid%xp%v1x(i,j,k)
         end do

         ! Forward 
         do i=ide-1,ids+1,-1
            grid%xp%v1x(i-1,j,k) = grid%xp%v1x(i-1,j,k) + w(1)*grid%xp%v1x(i,j,k)
            grid%xp%v1x(i+1,j,k) = grid%xp%v1x(i+1,j,k) + w(3)*grid%xp%v1x(i,j,k)
            grid%xp%v1x(i,j  ,k) = w(2)*grid%xp%v1x(i,j,k)
         end do
      end do
   end do

   ! Apply (i,j',k' -> i',j',k) transpose (v1x -> v1z).
   call da_transpose_x2z (grid)

   var(its:ite,jts:jte,kts:kte) = grid%xp%v1z(its:ite,jts:jte,kts:kte)

   if (trace_use) call da_trace_exit("da_filter_adj")

end subroutine da_filter_adj


subroutine da_wdt(h,w,terr,dt)

   !----------------------------------------------------------------------
   ! Purpose: Calculate DT
   !----------------------------------------------------------------------

   implicit none

   real, intent(in)  :: h(kts:kte)
   real, intent(out) :: dt(kts:kte)
   real, intent(in)  :: w(kts:kte+1)
   real, intent(in)  :: terr

   integer :: k

   if (trace_use) call da_trace_entry("da_wdt")

   do k=kte,kts+1,-1
      if (w(k) >= 0.1) then
         dt(k)=(h(k)-h(k-1))/w(k)
      else
         dt(k)=0.0
      end if
   end do

   if (w(kts) >= 0.1) then
      dt(kts)=(h(kts)-terr)/w(kts)
   else
      dt(kts)=0.0
   end if

   if (trace_use) call da_trace_exit("da_wdt")

end subroutine da_wdt


subroutine da_integrat_dz(grid)

   !---------------------------------------------------------------------------
   ! Non-linear PW forward operator.
   ! ===============================
   !
   ! Purpose: To calculate the IWV from the model QV and PP, TT.
   !
   ! Method:  IWV = sum {QV * RHO * dZ}
   !
   !           Unit: Qv (Kg/Kg), RHO(Kg/M^3), dZ(M)
   !                 PW (cm)
   !
   ! input     : QV, PP, TT
   !
   ! output    : PW
   !
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   integer :: i, j, K, ij 

   real    :: pw

   if (trace_use) call da_trace_entry("da_integrat_dz")

   ! weighted sum of vertical column 

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (ij, i, j, pw)
   do ij = 1, grid%num_tiles

   do j=jts, jte
      do i=its, ite
         pw = 0.0
         do k=kts, kte
            pw = pw + (grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k)) * grid%xb%q(i,j,k)*grid%xb%rho(i,j,k)
         end do

         grid%xb%tpw(i,j) = 0.1*pw
      end do
   end do

   end do
   !$OMP END PARALLEL DO

   if (trace_use) call da_trace_exit("da_integrat_dz")

end subroutine da_integrat_dz


subroutine da_uv_to_sd_lin(spd,dir,u,v,ub,vb)
   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)     :: u, v
   real, intent(inout)  :: ub, vb
   real, intent(out)    :: spd, dir
   
   if (trace_use_dull) call da_trace_entry("da_uv_to_sd_lin")

   if ( (ub*ub+vb*vb) == 0.0 ) then      ! Avoid division by zero
      spd = 0.0
      dir = 0.0
      return
      if (trace_use_dull) call da_trace_exit("da_uv_to_sd_lin")
   end if

   if (abs(ub - 0.0) <= 0.1) ub = (ub/abs(ub))*0.1
   if (abs(vb - 0.0) <= 0.1) vb = (vb/abs(vb))*0.1

   spd = (ub*u+vb*v)/sqrt(ub*ub+vb*vb)
   dir = (vb*u-ub*v)/(ub*ub+vb*vb) * 180.0/pi

   if (trace_use_dull) call da_trace_exit("da_uv_to_sd_lin")

end subroutine da_uv_to_sd_lin 
subroutine da_uv_to_sd_adj(spd,dir,u,v,ub,vb)
   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)     ::  spd, dir
   real, intent(inout)  ::  ub, vb
   real, intent(out)    ::  u, v 
 
   if (trace_use_dull) call da_trace_entry("da_uv_to_sd_adj")

   u = 0.0
   v = 0.0

   if( (ub*ub+vb*vb) == 0.0 ) then      ! Avoid division by zero
      return
      if (trace_use_dull) call da_trace_exit("da_uv_to_sd_adj")
   end if

   if (abs(ub - 0.0) <= 0.1) ub = (ub/abs(ub))*0.1
   if (abs(vb - 0.0) <= 0.1) vb = (vb/abs(vb))*0.1


   u = ub/sqrt(ub*ub+vb*vb)*spd + vb/(ub*ub+vb*vb)*dir*180.0/pi
   v = vb/sqrt(ub*ub+vb*vb)*spd - ub/(ub*ub+vb*vb)*dir*180.0/pi

   if (trace_use_dull) call da_trace_exit("da_uv_to_sd_adj")

end subroutine da_uv_to_sd_adj 

end module da_physics

