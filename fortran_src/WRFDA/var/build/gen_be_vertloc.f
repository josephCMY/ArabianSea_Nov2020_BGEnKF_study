












program gen_be_vertloc

   use da_gen_be, only : da_eof_decomposition

   implicit none

   integer, parameter    :: ni = 1000 
   integer, parameter    :: stdout = 6 
   character (len=3)     :: cnk                       
   real*8, parameter       :: rscale = 0.1
   real*8, parameter       :: var_threshold = 0.99

   integer               :: i, k, k1, m,nk,nk2
   integer               :: nm
   integer               :: ktarget
   real*8                  :: ni1_inv
   real*8                  :: nk_inv, nk3_inv
   real*8                  :: kscale, kscale_invsq
   real*8                  :: kdist
   real*8                  :: r
   real*8                  :: totvar, totvar_inv, cumvar
   real*8,allocatable      :: rho(:,:)
   real*8,allocatable      :: e(:,:)
   real,allocatable      :: cov(:,:)
   real*8,allocatable      :: eval(:)
   real*8,allocatable      :: evec(:,:)
   real*8,allocatable      :: v(:)
   real*8,allocatable      :: x(:)
   
  cnk=""
  call getarg( 1, cnk )
  read(cnk,'(i3)')nk2
  nk=nk2-1
  write(stdout,'(a,i6)') ' vertical level = ', nk2

  allocate(rho(1:nk,1:nk))
  allocate(e(1:ni,1:nk))
  allocate(cov(1:nk,1:nk))
  allocate(eval(1:nk))
  allocate(evec(1:nk,1:nk))
  allocate(v(1:nk))
  allocate(x(1:nk))
   ni1_inv = 1.0 / real (ni - 1)
   nk_inv = 1.0 / real (nk)
   nk3_inv = 1.0 / real (nk-3)


   do k = 1, nk
      kscale = 10.0 * real(k) * nk_inv  


      kscale_invsq = 1.0 / ( kscale * kscale )
      do k1 = k, nk
         kdist = k1 - k
         rho(k,k1) = exp ( -real(kdist * kdist) * kscale_invsq )
         rho(k1,k) = rho(k,k1)
      end do
   end do
   cov = rho




























   call da_eof_decomposition( nk, cov, evec, eval )


   totvar = sum(eval(1:nk))
   do k = 1, nk
     if ( eval(k) < 0.0 ) eval(k) = 0.0
   end do
   totvar_inv = 1.0 / sum(eval(1:nk))
   eval(:) = eval(:) * totvar * totvar_inv 


   nm = 0
   totvar_inv = 1.0 / sum(eval(1:nk))
   do k = 1, nk
      cumvar = sum(eval(1:k)) * totvar_inv
      if ( nm == 0 .and. cumvar >= var_threshold ) nm = k
   end do

   write(stdout,'(a,f15.5,i6)')' Threshold, truncation = ', var_threshold, nm

   open(10, file = 'be.vertloc.dat', status = 'unknown', form='unformatted')
   write(10) nk
   eval(:) = sqrt(eval(:)) 
   write(10)eval(1:nk)
   write(10)evec(1:nk,1:nk)

   close(10)





   ktarget = nk
   x(:) = 0.0
   x(ktarget) = 1.0


   do m = 1, nm
      v(m) = eval(m) * sum(evec(1:nk,m) * x(1:nk))
   end do


   do k = 1, nk
      x(k) = sum(evec(k,1:nm) * eval(1:nm) * v(1:nm))

   end do
   
  deallocate(rho)
  deallocate(e)
  deallocate(cov)
  deallocate(eval)
  deallocate(evec)
  deallocate(v)
  deallocate(x)

end program gen_be_vertloc

