












module da_interpolation

   use da_control, only : stdout, trace_use, trace_use_frequent, missing_r, &
      anal_type_verify, v_interp_h, v_interp_p,ims,ime,jms,jme,kms,kme, &
      kts,kte, trace_use_dull, interp_option
   use da_define_structures, only : infa_type
   use da_tools, only : da_togrid
   use da_tracing, only : da_trace_entry, da_trace_exit
   use da_reporting, only: da_error, da_warning, message

   implicit none

contains

subroutine da_to_zk(obs_v, mdl_v, v_interp_optn, zk)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                  intent(in)  :: v_interp_optn
   real,                     intent(in)  :: obs_v
   real, dimension(kms:kme), intent(in)  :: mdl_v
   real,                     intent(out) :: zk

   integer                :: k

   if (trace_use_dull) call da_trace_entry("da_to_zk")

   zk = missing_r

   if (v_interp_optn == v_interp_p) then
      if (obs_v > mdl_v(kts) .or. obs_v < mdl_v(kte)) then
         if (anal_type_verify) then
            ! Guo (02/06/2006), for VERifY, allow the extrapolation to obtain the zk:
            if (obs_v > mdl_v(kts)) then
               ! below the lowest level:
               zk = real(kts+1) - &
                  (obs_v-mdl_v(kts+1))/(mdl_v(kts)-mdl_v(kts+1))
            else if (obs_v < mdl_v(kte)) then
               ! above the highest level:
               zk = real(kte-1) + &
                  (obs_v-mdl_v(kte-1))/(mdl_v(kte)-mdl_v(kte-1))
            end if
         else
            if (trace_use_dull) call da_trace_exit("da_to_zk")
            return
         end if
      else
         do k = kts,kte-1
            if(obs_v <= mdl_v(k) .and. obs_v >= mdl_v(k+1)) then
               zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
               exit
            end if
         end do
      end if
   else if(v_interp_optn == v_interp_h) then
      if (obs_v < mdl_v(kts) .or. obs_v > mdl_v(kte)) then
         if (anal_type_verify) then
            ! Guo (02/06/2006), for VERifY, allow the extrapolation to obtain the zk:
            if (obs_v < mdl_v(kts)) then
               ! below the lowest level:
               zk = real(kts+1) - &
                  (obs_v-mdl_v(kts+1))/(mdl_v(kts)-mdl_v(kts+1))
            else if (obs_v > mdl_v(kte)) then
               ! above the highest level:
               zk = real(kte-1) + &
                  (obs_v-mdl_v(kte-1))/(mdl_v(kte)-mdl_v(kte-1))
            end if
         else
            if (trace_use_dull) call da_trace_exit("da_to_zk")
            return
         end if
      else
         do k = kts,kte-1
            if(obs_v >= mdl_v(k) .and. obs_v <= mdl_v(k+1)) then
               zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
               exit
            end if
         end do
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_to_zk")

end subroutine da_to_zk


subroutine da_to_zk_new(obs_v, mdl_v, v_interp_optn, num,nlevels,zk)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real,             intent(in)  :: obs_v(nlevels)
   real,             intent(in)  :: mdl_v(kms:kme,num)
   integer,          intent(in)  :: v_interp_optn
   integer,          intent(in)  :: num
   integer,          intent(in)  :: nlevels
   real,             intent(out) :: zk(nlevels,num)

   integer                :: k,n,kk
   real    :: r_ktsplus, r_kteminus

   if (trace_use) call da_trace_entry("da_to_zk_new")

   zk(:,:) = missing_r

   r_ktsplus  = real(kts+1)
   r_kteminus = real(kte-1)

   if (v_interp_optn == v_interp_p) then
      if (anal_type_verify) then
         do n=1,num
            do k=1,nlevels
               if (obs_v(k) > mdl_v(kts,n)) then
                  ! below the lowest level:
                  zk(k,n) = r_ktsplus - &
                     (obs_v(k)-mdl_v(kts+1,n))/(mdl_v(kts,n)-mdl_v(kts+1,n))
               else if (obs_v(k) < mdl_v(kts,n)) then
                  ! above the highest level:
                  zk(k,n) = r_kteminus + &
                     (obs_v(k)-mdl_v(kte-1,n))/(mdl_v(kte,n)-mdl_v(kte-1,n))
               else
                  do kk = kts,kte-1
                     if (obs_v(k) <= mdl_v(kk,n) .and. obs_v(k) >= mdl_v(kk+1,n)) then
                        zk(k,n) = real(kk) + (mdl_v(kk,n) - obs_v(k))/(mdl_v(kk,n) - mdl_v(kk+1,n))
                        exit
                     end if
                  end do
               end if
            end do
         end do
      else
         do n=1,num
            do k=1,nlevels
               do kk = kts,kte-1
                  if (obs_v(k) <= mdl_v(kk,n) .and. obs_v(k) >= mdl_v(kk+1,n)) then
                     zk(k,n) = real(kk) + (mdl_v(kk,n) - obs_v(k))/(mdl_v(kk,n) - mdl_v(kk+1,n))
                     exit
                  end if
               end do
            end do
         end do
      end if
   else if(v_interp_optn == v_interp_h) then
      if (anal_type_verify) then
         do n=1,num
            do k=1,nlevels
               if (obs_v(k) < mdl_v(kts,n)) then
                   ! below the lowest level:
                   zk(k,n) = r_ktsplus - &
                     (obs_v(k)-mdl_v(kts+1,n))/(mdl_v(kts,n)-mdl_v(kts+1,n))
               else if (obs_v(k) > mdl_v(kts,n)) then
                  ! above the highest level:
                  zk(k,n) = r_kteminus + &
                     (obs_v(k)-mdl_v(kte-1,n))/(mdl_v(kte,n)-mdl_v(kte-1,n))
               else
                  do kk = kts,kte-1
                     if (obs_v(k) >= mdl_v(kk,n) .and. obs_v(k) <= mdl_v(kk+1,n)) then
                        zk(k,n) = real(kk) + (mdl_v(kk,n) - obs_v(k))/(mdl_v(kk,n) - mdl_v(kk+1,n))
                        exit
                     end if
                  end do
               end if
            end do
         end do
      else
         do n=1,num
            do k=1,nlevels
               do kk = kts,kte-1
                  if (obs_v(k) >= mdl_v(kk,n) .and. obs_v(k) <= mdl_v(kk+1,n)) then
                     zk(k,n) = real(kk) + (mdl_v(kk,n) - obs_v(k))/(mdl_v(kk,n) - mdl_v(kk+1,n))
                     exit
                  end if
               end do
            end do
         end do
      end if
   end if

   if (trace_use) call da_trace_exit("da_to_zk_new")

end subroutine da_to_zk_new



subroutine da_interp_2d_partial(fm2d, info, k, n1, n2, fo2d, mask2d)

   !-----------------------------------------------------------------------
   ! Interpolation option (1=linear)
   !                      (2=quadratic)
   !                      (3=masked average)
   !-----------------------------------------------------------------------

   implicit none

   real,            intent(in)  :: fm2d(ims:ime,jms:jme)   ! Input variable
   type(infa_type), intent(in)  :: info
   integer,         intent(in)  :: k                       ! level
   integer,         intent(in)  :: n1,n2                   ! Range of obs
   real,            intent(out) :: fo2d(n1:n2)             ! Output variable 
   integer, optional, intent(in):: mask2d(ims:ime,jms:jme) ! Field mask (1=valid, 0=invalid)
   
   if (trace_use_frequent) call da_trace_entry("da_interp_2d_partial")
   
   select case(interp_option)
   case (1) ;
      call da_interp_lin_2d_partial(fm2d, info, k, n1, n2, fo2d)
   case (2) ;
      call da_interp_quad_2d_partial(fm2d, info, k, n1, n2, fo2d)
   case (3) ;
      if (.not. present(mask2d)) call da_error("da_interp_2d_partial.inc",26,(/"da_interp_2d_partial: no mask given for masked interpolation"/))
      call da_interp_msk_avg_2d_partial(fm2d, info, k, n1, n2, fo2d, mask2d)
   case default;
      write(unit=message(1),fmt='(a,i4)') 'Interpolation option not supported: ', interp_option
      call da_warning("da_interp_2d_partial.inc",30,message(1:1))
   end select   
      
   if (trace_use_frequent) call da_trace_exit("da_interp_2d_partial")

end subroutine da_interp_2d_partial
subroutine da_interp_lin_2d_partial(fm2d, info, k, n1, n2, fo2d)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real,            intent(in)  :: fm2d(ims:ime,jms:jme) ! Input variable
   type(infa_type), intent(in)  :: info
   integer,         intent(in)  :: k                     ! level
   integer,         intent(in)  :: n1,n2                 ! Range of obs
   real,            intent(out) :: fo2d(n1:n2)           ! Output variable 
   
   integer :: n

   if (trace_use_frequent) call da_trace_entry("da_interp_lin_2d_partial")

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( n )
   do n=n1,n2
      fo2d(n) = info%dym(k,n)*(info%dxm(k,n)*fm2d(info%i(k,n),info%j(k,n))   + info%dx(k,n)*fm2d(info%i(k,n)+1,info%j(k,n))) &
              + info%dy(k,n) *(info%dxm(k,n)*fm2d(info%i(k,n),info%j(k,n)+1) + info%dx(k,n)*fm2d(info%i(k,n)+1,info%j(k,n)+1))
   end do
   !$OMP END PARALLEL DO

   if (trace_use_frequent) call da_trace_exit("da_interp_lin_2d_partial")

end subroutine da_interp_lin_2d_partial


subroutine da_interp_lin_2d(fm2d, info, k, fo2d)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real,            intent(in)  :: fm2d(ims:ime,jms:jme) ! Input variable
   type(infa_type), intent(in)  :: info
   integer,         intent(in)  :: k                     ! level
   real,            intent(out) :: fo2d(info%n1:info%n2)           ! Output variable 
   
   integer :: n

   if (trace_use_frequent) call da_trace_entry("da_interp_lin_2d")

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (n) 
   do n=info%n1,info%n2
      fo2d(n) = info%dym(k,n)*(info%dxm(k,n)*fm2d(info%i(k,n),info%j(k,n))   + info%dx(k,n)*fm2d(info%i(k,n)+1,info%j(k,n))) &
              + info%dy(k,n) *(info%dxm(k,n)*fm2d(info%i(k,n),info%j(k,n)+1) + info%dx(k,n)*fm2d(info%i(k,n)+1,info%j(k,n)+1))
   end do
   !$OMP END PARALLEL DO

   if (trace_use_frequent) call da_trace_exit("da_interp_lin_2d")

end subroutine da_interp_lin_2d


subroutine da_interp_lin_2d_adj_partial(fm2d, info, k1,k2, fo2d)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real,            intent(inout) :: fm2d(ims:ime,jms:jme,k1:k2) ! 'Input' variable modified by adjoint
   type(infa_type), intent(in)    :: info
   integer,         intent(in)    :: k1,k2                 ! Range of levels
   real,            intent(in)    :: fo2d(k1:k2,info%n1:info%n2)           ! 'Output' variable unchanged by adjoint

   integer :: n,k

   if (trace_use) call da_trace_entry("da_interp_lin_2d_adj_partial")

   do n=info%n1,info%n2
      do k=k1,k2
         fm2d(info%i(k,n)  ,info%j(k,n),k)   = info%dym(k,n) * info%dxm(k,n) * fo2d(k,n) + fm2d(info%i(k,n)  ,info%j(k,n),k)
         fm2d(info%i(k,n)+1,info%j(k,n),k)   = info%dym(k,n) * info%dx(k,n)  * fo2d(k,n) + fm2d(info%i(k,n)+1,info%j(k,n),k)
         fm2d(info%i(k,n)  ,info%j(k,n)+1,k) = info%dy(k,n)  * info%dxm(k,n) * fo2d(k,n) + fm2d(info%i(k,n)  ,info%j(k,n)+1,k)
         fm2d(info%i(k,n)+1,info%j(k,n)+1,k) = info%dy(k,n)  * info%dx(k,n)  * fo2d(k,n) + fm2d(info%i(k,n)+1,info%j(k,n)+1,k)
      end do
   end do

   if (trace_use) call da_trace_exit("da_interp_lin_2d_adj_partial")

end subroutine da_interp_lin_2d_adj_partial


subroutine da_interp_lin_2d_adj(fm2d, info, k, fo2d)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real,            intent(inout) :: fm2d(ims:ime,jms:jme) ! 'Input' variable modified by adjoint
   type(infa_type), intent(in)    :: info
   integer,         intent(in)    :: k                     ! level
   real,            intent(in)    :: fo2d(info%n1:info%n2)           ! 'Output' variable unchanged by adjoint

   integer :: n

   if (trace_use) call da_trace_entry("da_interp_lin_2d_adj")

   do n=info%n1,info%n2
      fm2d(info%i(k,n)  ,info%j(k,n))   = info%dym(k,n) * info%dxm(k,n) * fo2d(n) + fm2d(info%i(k,n)  ,info%j(k,n))
      fm2d(info%i(k,n)+1,info%j(k,n))   = info%dym(k,n) * info%dx(k,n)  * fo2d(n) + fm2d(info%i(k,n)+1,info%j(k,n))
      fm2d(info%i(k,n)  ,info%j(k,n)+1) = info%dy(k,n)  * info%dxm(k,n) * fo2d(n) + fm2d(info%i(k,n)  ,info%j(k,n)+1)
      fm2d(info%i(k,n)+1,info%j(k,n)+1) = info%dy(k,n)  * info%dx(k,n)  * fo2d(n) + fm2d(info%i(k,n)+1,info%j(k,n)+1)
   end do

   if (trace_use) call da_trace_exit("da_interp_lin_2d_adj")

end subroutine da_interp_lin_2d_adj


subroutine da_interp_lin_3d(fm3d, info, fo3d)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   real,            intent(in)    :: fm3d(ims:ime,jms:jme,kms:kme) ! Input variable 
   type(infa_type), intent(in)    :: info       
   real,            intent(inout) :: fo3d(1:info%max_lev,info%n1:info%n2)           ! Output variable 

   integer :: n, k
   real    :: fmz(kms:kme)

   if (trace_use) call da_trace_entry("da_interp_lin_3d")

   fo3d(:,:) = 0.0

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( n, fmz, k )
   do n=info%n1,info%n2
      fmz(:)=0.0

      do k = kts,kte
         fmz(k) = &
              info%dym(k,n) * (info%dxm(k,n)*fm3d(info%i(k,n), info%j(k,n), k) &
            + info%dx (k,n) * fm3d(info%i(k,n)+1,info%j(k,n), k)) &
            + info%dy (k,n) * (info%dxm(k,n)*fm3d(info%i(k,n), info%j(k,n)+1, k) &
            + info%dx (k,n) * fm3d(info%i(k,n)+1, info%j(k,n)+1, k))
      end do
      do k = 1, info%levels(n)
         if (info%k(k,n) > 0) then
            fo3d(k,n) = info%dzm(k,n)*fmz(info%k(k,n)) + info%dz(k,n)*fmz(info%k(k,n)+1)
         end if
      end do
   end do
   !$OMP END PARALLEL DO

   if (trace_use) call da_trace_exit("da_interp_lin_3d")

end subroutine da_interp_lin_3d
subroutine da_interp_lin_3d_adj(fm3d, info, fo3d)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   real,            intent(inout) :: fm3d(ims:ime,jms:jme,kms:kme)
   type(infa_type), intent(in)    :: info
   real,            intent(in)    :: fo3d(1:info%max_lev,info%n1:info%n2) 

   integer :: n,k
   real    :: fmz(kms:kme)

   if (trace_use) call da_trace_entry("da_interp_lin_3d_adj")

   do n=info%n1,info%n2
      fmz = 0.0
      do k = 1, info%levels(n)
         if (info%k(k,n) > 0) then
            fmz(info%k(k,n))   = info%dzm(k,n) * fo3d(k,n) + fmz(info%k(k,n))
            fmz(info%k(k,n)+1) = info%dz (k,n) * fo3d(k,n) + fmz(info%k(k,n)+1)
         end if
      end do
  
      do k=kts,kte
        fm3d(info%i(k,n)  ,info%j(k,n),k)   = info%dym(k,n)*info%dxm(k,n)*fmz(k) + fm3d(info%i(k,n)  ,info%j(k,n)  ,k)
         fm3d(info%i(k,n)+1,info%j(k,n),k)   = info%dym(k,n)*info%dx (k,n)*fmz(k) + fm3d(info%i(k,n)+1,info%j(k,n)  ,k)
         fm3d(info%i(k,n)  ,info%j(k,n)+1,k) = info%dy (k,n)*info%dxm(k,n)*fmz(k) + fm3d(info%i(k,n)  ,info%j(k,n)+1,k)
         fm3d(info%i(k,n)+1,info%j(k,n)+1,k) = info%dy (k,n)*info%dx (k,n)*fmz(k) + fm3d(info%i(k,n)+1,info%j(k,n)+1,k)
      end do
   end do

   if (trace_use) call da_trace_exit("da_interp_lin_3d_adj")

end subroutine da_interp_lin_3d_adj
subroutine da_interp_quad_2d_partial(fm2d, info, k, n1, n2, fo2d)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real,            intent(in)  :: fm2d(ims:ime,jms:jme) ! Input variable
   type(infa_type), intent(in)  :: info
   integer,         intent(in)  :: k                     ! level
   integer,         intent(in)  :: n1,n2                 ! Range of obs
   real,            intent(out) :: fo2d(n1:n2)           ! Output variable 
   
   integer :: n

   if (trace_use_frequent) call da_trace_entry("da_interp_quad_2d_partial")

   do n=n1,n2
      fo2d(n) = info%dym(k,n)*(info%dxm(k,n)*fm2d(info%i(k,n),info%j(k,n))   + info%dx(k,n)*fm2d(info%i(k,n)+1,info%j(k,n))) &
              + info%dy(k,n) *(info%dxm(k,n)*fm2d(info%i(k,n),info%j(k,n)+1) + info%dx(k,n)*fm2d(info%i(k,n)+1,info%j(k,n)+1))
   end do

   if (trace_use_frequent) call da_trace_exit("da_interp_quad_2d_partial")

end subroutine da_interp_quad_2d_partial


subroutine da_interp_msk_avg_2d_partial(fm2d, info, k, n1, n2, fo2d, mask2d)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real,            intent(in)  :: fm2d(ims:ime,jms:jme)  ! Input variable
   type(infa_type), intent(in)  :: info
   integer,         intent(in)  :: k                      ! level
   integer,         intent(in)  :: n1,n2                  ! Range of obs
   real,            intent(out) :: fo2d(n1:n2)            ! Output variable 
   integer, optional, intent(in):: mask2d(ims:ime,jms:jme)! Input variable mask (1=valid, 0=invalid)
   
   integer :: n

   if (trace_use_frequent) call da_trace_entry("da_interp_msk_avg_2d_partial")

   do n=n1,n2
      fo2d(n) = info%dym(k,n)*(info%dxm(k,n)*fm2d(info%i(k,n),info%j(k,n))   + info%dx(k,n)*fm2d(info%i(k,n)+1,info%j(k,n))) &
              + info%dy(k,n) *(info%dxm(k,n)*fm2d(info%i(k,n),info%j(k,n)+1) + info%dx(k,n)*fm2d(info%i(k,n)+1,info%j(k,n)+1))
   end do

   if (trace_use_frequent) call da_trace_exit("da_interp_msk_avg_2d_partial")

end subroutine da_interp_msk_avg_2d_partial



end module da_interpolation

