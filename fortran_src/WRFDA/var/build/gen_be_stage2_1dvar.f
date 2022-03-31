












program gen_be_stage2_1dvar










   use da_control, only : filename_len
   use da_gen_be, only : da_eof_decomposition, da_eof_decomposition_test
   use da_tools_serial, only : da_get_unit,da_advance_cymdh

   implicit none

   real, parameter     :: t_factor = 1.0             
   real, parameter     :: q_factor = 1000.0          
   real, parameter     :: ps_factor = 0.01           
   real, parameter     :: u_factor = 1.0             

   character*10        :: start_date, end_date       
   character*10        :: date, new_date             
   character*10        :: variable                   
   character(len=filename_len)        :: dat_dir                    
   character*80        :: expt                       
   character(len=filename_len)        :: filename                   
   character*3         :: ce                         
   integer             :: ni, nj, nk, nkdum          
   integer             :: nk1                        
   integer             :: nkbe                       
   integer             :: i, j, k, member, k2        
   integer             :: kbe, k2be                  
   integer             :: b                          
   integer             :: sdate, cdate, edate        
   integer             :: interval                   
   integer             :: ne                         
   integer             :: bin_type                   
   integer             :: num_bins                   
   integer             :: num_bins2d                 
   real                :: lat_min, lat_max           
   real                :: binwidth_lat               
   real                :: hgt_min, hgt_max           
   real                :: binwidth_hgt               
   real                :: coeffa, coeffb             
   logical             :: first_time                 
   logical             :: ldum1, ldum2               
   logical             :: testing_eofs               

   real, allocatable   :: temp(:,:,:)                
   real, allocatable   :: q(:,:,:)                   
   real, allocatable   :: ps(:,:)                    
   real, allocatable   :: u10(:,:)                   
   real, allocatable   :: v10(:,:)                   
   integer, allocatable:: bin(:,:,:)                 
   integer, allocatable:: bin2d(:,:)                 
   integer, allocatable:: bin_pts2d(:)               
   real, allocatable   :: covar_tt(:,:,:)            
   real, allocatable   :: covar_tq(:,:,:)            
   real, allocatable   :: covar_tps(:,:)             
   real, allocatable   :: covar_tu10(:,:)            
   real, allocatable   :: covar_tv10(:,:)            
   real, allocatable   :: covar_qq(:,:,:)            
   real, allocatable   :: covar_qps(:,:)             
   real, allocatable   :: covar_qu10(:,:)            
   real, allocatable   :: covar_qv10(:,:)            
   real, allocatable   :: covar_psps(:)              
   real, allocatable   :: covar_psu10(:)             
   real, allocatable   :: covar_psv10(:)             
   real, allocatable   :: covar_u10u10(:)            
   real, allocatable   :: covar_u10v10(:)            
   real, allocatable   :: covar_v10v10(:)            
   real, allocatable   :: be(:,:,:)                  
   real, allocatable   :: be_inv(:,:,:)              
   real, allocatable   :: be_loc(:,:)                
   real, allocatable   :: be_loc_inv(:,:)            
   real, allocatable   :: be_test(:,:)               

   namelist / gen_be_stage2_1dvar_nl / start_date, end_date, interval, &
                                       ne, testing_eofs, expt, dat_dir

   integer :: ounit,iunit,namelist_unit


   write(6,'(a)')' [1] Initialize namelist variables and other scalars.'


   call da_get_unit(ounit)
   call da_get_unit(iunit)
   call da_get_unit(namelist_unit)

   start_date = '2004030312'
   end_date = '2004033112'
   interval = 24
   ne = 1
   testing_eofs = .false.
   expt = 'gen_be_stage2_1dvar'
   dat_dir = '/mmmtmp1/dmbarker'

   open(unit=namelist_unit, file='gen_be_stage2_1dvar_nl.nl', &
        form='formatted', status='old', action='read')
   read(namelist_unit, gen_be_stage2_1dvar_nl)
   close(namelist_unit)

   read(start_date(1:10), fmt='(i10)')sdate
   read(end_date(1:10), fmt='(i10)')edate
   write(6,'(4a)')' Computing multivariate covariances for 1D-Var'
   write(6,'(4a)') ' Time period is ', start_date, ' to ', end_date
   write(6,'(a,i8,a)')' Interval between dates = ', interval, 'hours.'
   write(6,'(a,i8)')' Number of ensemble members at each time = ', ne

   date = start_date
   cdate = sdate


   write(6,'(2a)')' [2] Read fields, and calculate correlations'


   write(6,'(/a)')' Unbiased perturbations read in have been scaled by the following factors:'
   write(6,'(a,f15.5)')' T scaling factor = ', t_factor
   write(6,'(a,f15.5)')' q scaling factor = ', q_factor
   write(6,'(a,f15.5)')' ps scaling factor = ', ps_factor
   write(6,'(a,f15.5)')' u scaling factor = ', u_factor

   first_time = .true.

   do while ( cdate <= edate )
      write(6,'(a,a)')'    Processing data for date ', date

      do member = 1, ne

         write(ce,'(i3.3)')member


         variable = 't'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         nk1 = nk + 1

         if ( first_time ) then
            write(6,'(a,3i8)')'    i, j, k dimensions are ', ni, nj, nk1
            write(6,'(a,3i8)')'    Note k dimension expanded by one to include T2m, q2m in array.'

            allocate( bin(1:ni,1:nj,1:nk) ) 
            allocate( temp(1:ni,1:nj,1:nk1) )
            allocate( q(1:ni,1:nj,1:nk1) )
            allocate( bin2d(1:ni,1:nj) )
            allocate( ps(1:ni,1:nj) )
            allocate( u10(1:ni,1:nj) )
            allocate( v10(1:ni,1:nj) )
         end if

         read(iunit)temp(1:ni,1:nj,2:nk+1)
         close(iunit)

         if ( first_time ) then

            filename = 'bin.data'
            open (iunit, file = filename, form='unformatted')
            read(iunit)bin_type
            read(iunit)lat_min, lat_max, binwidth_lat
            read(iunit)hgt_min, hgt_max, binwidth_hgt
            read(iunit)num_bins, num_bins2d
            read(iunit)bin(1:ni,1:nj,1:nk) 
            read(iunit)bin2d(1:ni,1:nj)
            close(iunit)

            allocate( bin_pts2d(1:num_bins2d) )
            allocate( covar_tt(1:nk1,1:nk1,1:num_bins2d) )     
            allocate( covar_tq(1:nk1,1:nk1,1:num_bins2d) )     
            allocate( covar_tps(1:nk1,1:num_bins2d) )
            allocate( covar_tu10(1:nk1,1:num_bins2d) )
            allocate( covar_tv10(1:nk1,1:num_bins2d) )
            allocate( covar_qq(1:nk1,1:nk1,1:num_bins2d) )     
            allocate( covar_qps(1:nk1,1:num_bins2d) )
            allocate( covar_qu10(1:nk1,1:num_bins2d) )
            allocate( covar_qv10(1:nk1,1:num_bins2d) )
            allocate( covar_psps(1:num_bins2d) )
            allocate( covar_psu10(1:num_bins2d) )
            allocate( covar_psv10(1:num_bins2d) )
            allocate( covar_u10u10(1:num_bins2d) )
            allocate( covar_u10v10(1:num_bins2d) )
            allocate( covar_v10v10(1:num_bins2d) )
            bin_pts2d(:) = 0
            covar_tt(:,:,:) = 0.0
            covar_tq(:,:,:) = 0.0
            covar_tps(:,:) = 0.0
            covar_tu10(:,:) = 0.0
            covar_tv10(:,:) = 0.0
            covar_qq(:,:,:) = 0.0
            covar_qps(:,:) = 0.0
            covar_qu10(:,:) = 0.0
            covar_qv10(:,:) = 0.0
            covar_psps(:) = 0.0
            covar_psu10(:) = 0.0
            covar_psv10(:) = 0.0
            covar_u10u10(:) = 0.0
            covar_u10v10(:) = 0.0
            covar_v10v10(:) = 0.0
            first_time = .false.
         end if


         variable = 't2'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)ldum1, ldum2 
         read(iunit)temp(1:ni,1:nj,1) 
         close(iunit)


         variable = 'q'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)q(1:ni,1:nj,2:nk+1)
         close(iunit)


         variable = 'q2'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)ldum1, ldum2 
         read(iunit)q(1:ni,1:nj,1) 
         close(iunit)


         variable = 'ps'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)ldum1, ldum2 
         read(iunit)ps
         close(iunit)


         variable = 'u10'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)ldum1, ldum2 
         read(iunit)u10
         close(iunit)


         variable = 'v10'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)ldum1, ldum2 
         read(iunit)v10
         close(iunit)



         do j = 1, nj
            do i = 1, ni
               b = bin2d(i,j)
               bin_pts2d(b) = bin_pts2d(b) + 1
               coeffa = 1.0 / real(bin_pts2d(b))
               coeffb = real(bin_pts2d(b)-1) * coeffa

               do k = 1, nk1


                  do k2 = 1, k
                     covar_tt(k,k2,b) = coeffb * covar_tt(k,k2,b) + &
                                        coeffa * temp(i,j,k) * temp(i,j,k2)
                  end do


                  do k2 = 1, nk1
                     covar_tq(k,k2,b) = coeffb * covar_tq(k,k2,b) + &
                                        coeffa * temp(i,j,k) * q(i,j,k2)
                  end do


                  covar_tps(k,b) = coeffb * covar_tps(k,b) + coeffa * temp(i,j,k) * ps(i,j)


                  covar_tu10(k,b) = coeffb * covar_tu10(k,b) + coeffa * temp(i,j,k) * u10(i,j)


                  covar_tv10(k,b) = coeffb * covar_tv10(k,b) + coeffa * temp(i,j,k) * v10(i,j)


                  do k2 = 1, k
                     covar_qq(k,k2,b) = coeffb * covar_qq(k,k2,b) + &
                                        coeffa * q(i,j,k) * q(i,j,k2)
                  end do


                  covar_qps(k,b) = coeffb * covar_qps(k,b) + coeffa * q(i,j,k) * ps(i,j)


                  covar_qu10(k,b) = coeffb * covar_qu10(k,b) + coeffa * q(i,j,k) * u10(i,j)


                  covar_qv10(k,b) = coeffb * covar_qv10(k,b) + coeffa * q(i,j,k) * v10(i,j)
               end do


               covar_psps(b) = coeffb * covar_psps(b) + coeffa * ps(i,j) * ps(i,j)
               covar_psu10(b) = coeffb * covar_psu10(b) + coeffa * ps(i,j) * u10(i,j)
               covar_psv10(b) = coeffb * covar_psv10(b) + coeffa * ps(i,j) * v10(i,j)
               covar_u10u10(b) = coeffb * covar_u10u10(b) + coeffa * u10(i,j) * u10(i,j)
               covar_u10v10(b) = coeffb * covar_u10v10(b) + coeffa * u10(i,j) * v10(i,j)
               covar_v10v10(b) = coeffb * covar_v10v10(b) + coeffa * v10(i,j) * v10(i,j)

            end do
         end do
      end do  


      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate
   end do     


   write(6,'(2a)')' [3] Pack into total background error covariance and write '


   nkbe = 2 * nk1 + 3
   allocate( be(1:nkbe,1:nkbe,1:num_bins2d) )
   be = 0.0

   do b = 1, num_bins2d


      do k = 1, nk1
         kbe = k
         do k2 = 1, k
            k2be = k2
            be(kbe,k2be,b) = covar_tt(k,k2,b)
         end do
      end do


      do k = 1, nk1
         kbe = nk1 + k
         do k2 = 1, nk1
            k2be = k2
            be(kbe,k2be,b) = covar_tq(k,k2,b)
         end do
      end do


      kbe = 2 * nk1 + 1
      do k2 = 1, nk1
         k2be = k2
         be(kbe,k2be,b) = covar_tps(k2,b)
      end do


      kbe = 2 * nk1 + 2
      do k2 = 1, nk1
         k2be = k2
         be(kbe,k2be,b) = covar_tu10(k2,b)
      end do


      kbe = 2 * nk1 + 3
      do k2 = 1, nk1
         k2be = k2
         be(kbe,k2be,b) = covar_tv10(k2,b)
      end do


      do k = 1, nk1
         kbe = nk1 + k
         do k2 = 1, k
            k2be = nk1 + k2
            be(kbe,k2be,b) = covar_qq(k,k2,b)
         end do
      end do


      kbe = 2 * nk1 + 1
      do k2 = 1, nk1
         k2be = nk1 + k2
         be(kbe,k2be,b) = covar_qps(k2,b)
      end do


      kbe = 2 * nk1 + 2
      do k2 = 1, nk1
         k2be = nk1 + k2
         be(kbe,k2be,b) = covar_qu10(k2,b)
      end do


      kbe = 2 * nk1 + 3
      do k2 = 1, nk1
         k2be = nk1 + k2
         be(kbe,k2be,b) = covar_qv10(k2,b)
      end do


      kbe = 2 * nk1 + 1
      k2be = 2 * nk1 + 1
      be(kbe,k2be,b) = covar_psps(b)


      kbe = 2 * nk1 + 2
      k2be = 2 * nk1 + 1
      be(kbe,k2be,b) = covar_psu10(b)


      kbe = 2 * nk1 + 3
      k2be = 2 * nk1 + 1
      be(kbe,k2be,b) = covar_psv10(b)


      kbe = 2 * nk1 + 2
      k2be = 2 * nk1 + 2
      be(kbe,k2be,b) = covar_u10u10(b)


      kbe = 2 * nk1 + 3
      k2be = 2 * nk1 + 2
      be(kbe,k2be,b) = covar_u10v10(b)


      kbe = 2 * nk1 + 3
      k2be = 2 * nk1 + 3
      be(kbe,k2be,b) = covar_v10v10(b)


      do k = 1, nkbe
         do k2 = k+1, nkbe
            be(k,k2,b) = be(k2,k,b)
         end do
      end do

   end do


   filename = 'gen_be_stage2_1dvar.dat'
   open (ounit, file = filename, form='unformatted')
   write(ounit)ni, nj, nk1
   write(ounit)num_bins2d
   write(ounit)be
   close(ounit)


   write(6,'(2a)')' [4] Calculate and write inverse covariances '


   allocate( be_inv(1:nkbe,1:nkbe,1:num_bins2d) )
   allocate( be_loc(1:nkbe,1:nkbe) )
   allocate( be_loc_inv(1:nkbe,1:nkbe) )
   allocate( be_test(1:nkbe,1:nkbe) )
   be_test = 0.0

   do b = 1, num_bins2d

      be_loc(:,:) = be(:,:,b)
      call da_get_be_inverse( testing_eofs, nkbe, be_loc, be_loc_inv )
      be_inv(:,:,b) = be_loc_inv(:,:)


      be_test = matmul(be_loc,be_loc_inv)


      do k = 1, nkbe
         do k2 = 1, nkbe
            write(52,'(f15.5)')be_test(k,k2)
         end do
      end do
   end do

   filename = 'gen_be_stage2_1dvar.inv.dat'
   open (ounit, file = filename, form='unformatted')
   write(ounit)ni, nj, nk1
   write(ounit)num_bins2d
   write(ounit)be_inv
   close(ounit)


   write(6,'(2a)')' [4] Read in again, and output for graphics '


   write(6,'(a,i10)')' num_bins2d = ', num_bins2d

   be = 0.0
   be_inv = 0.0


   filename = 'gen_be_stage2_1dvar.dat'
   open (iunit, file = filename, form='unformatted', status = 'old')
   read(iunit)ni, nj, nk1
   read(iunit)num_bins2d
   read(iunit)be
   close(iunit)

   do b = 1, num_bins2d
      do k = 1, nkbe
         do k2 = 1, nkbe
            write(50,'(f15.5)')be(k,k2,b)
         end do
      end do
   end do


   filename = 'gen_be_stage2_1dvar.inv.dat'
   open (iunit, file = filename, form='unformatted', status = 'old')
   read(iunit)ni, nj, nk1
   read(iunit)num_bins2d
   read(iunit)be_inv
   close(iunit)

   do b = 1, num_bins2d
      do k = 1, nkbe
         do k2 = 1, nkbe
            write(51,'(f15.5)')be_inv(k,k2,b)
         end do
      end do
   end do

contains

subroutine da_get_be_inverse( testing_eofs, nk, be, be_inv )
   


   
   implicit none

   real, parameter     :: variance_threshold = 1e-6  

   logical, intent(inout) :: testing_eofs            
   integer, intent(in)    :: nk                      
   real, intent(in)       :: be(1:nk,1:nk)           
   real, intent(out)      :: be_inv(1:nk,1:nk)       

   integer                :: m, k, k2                
   integer                :: mmax                    

   real                   :: summ                    
   real                   :: total_var_inv           
   real                   :: cumul_variance          

   real*8                 :: evec(1:nk,1:nk)         
   real*8                 :: eval(1:nk)              
   real                   :: eval_inv(1:nk)          
   real                   :: LamInvET(1:nk,1:nk)     


   call da_eof_decomposition( nk, be, evec, eval )


   if ( testing_eofs ) then
      call da_eof_decomposition_test( nk, be, evec, eval )
      testing_eofs = .false.
   end if


   summ = 0.0
   do m = 1, nk
      summ = summ + eval(m)
   end do
   total_var_inv = 1.0 / summ

   cumul_variance = 0.0
   mmax = nk
   do m = 1, nk
      cumul_variance = cumul_variance + eval(m) * total_var_inv
      if ( eval(m) < 0.0 .or. cumul_variance > 1.0 - variance_threshold ) then
         mmax = m - 1
         exit
      end if
   end do


   LamInvET(:,:) = 0.0
   eval_inv(1:mmax) = 1.0 / eval(1:mmax)
   do k = 1, nk
      do m = 1, mmax
         LamInvET(m,k) = evec(k,m) * eval_inv(m)
      end do
   end do


   do k = 1, nk
      do k2 = 1, k
         summ = 0.0
         do m = 1, nk
            summ = summ + evec(k,m) * LamInvET(m,k2)
         end do
         be_inv(k,k2) = summ
      end do
   end do


   do k = 1, nk
      do k2 = k+1, nk 
         be_inv(k,k2) = be_inv(k2,k)
      end do
   end do

end subroutine da_get_be_inverse

end program gen_be_stage2_1dvar
