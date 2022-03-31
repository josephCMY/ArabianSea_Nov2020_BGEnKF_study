












program gen_be_cov2d3d_contrib









   use da_control, only : stderr, stdout, filename_len
   use da_tools_serial, only : da_get_unit,da_advance_cymdh

   implicit none

   character*10        :: start_date, end_date       
   character*10        :: date, new_date             
   character*80        :: variable
   character*80        :: regcoeff_var
   character*10        :: variable1                  
   character*10        :: variable2                  

   real, allocatable   :: regcoeff_psi_chi(:)        
   real, allocatable   :: regcoeff_psi_t(:,:,:)      
   real, allocatable   :: regcoeff_psi_ps(:,:)       
   real, allocatable   :: regcoeff_psi_rh(:,:,:)     
   real, allocatable   :: regcoeff_chi_u_t(:,:,:)    
   real, allocatable   :: regcoeff_chi_u_ps(:,:)     
   real, allocatable   :: regcoeff_chi_u_rh(:,:,:)   
   real, allocatable   :: regcoeff_t_u_rh(:,:,:)     
   real, allocatable   :: regcoeff_ps_u_rh(:,:)      
   
   character(len=filename_len)        :: filename                   
   character*3         :: ce                         
   integer             :: ni, nj, nk, nkdum                 
   integer             :: i, j, k, member            
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

   real, allocatable   :: field1(:,:)              
   real, allocatable   :: field2(:,:,:)              
   real, allocatable   :: field2_balanced(:,:,:)   
   real*4, allocatable   :: field(:,:)   
   integer, allocatable:: bin(:,:,:)                 
   integer, allocatable:: bin2d(:,:)                 
   integer, allocatable:: bin_pts(:)                 
   real, allocatable   :: covar(:)                   
   real, allocatable   :: var(:)                     

   real, allocatable   :: regcoeff_field1_field2(:,:) 

   namelist / gen_be_cov3d_nl / start_date, end_date, interval, &
                                ne, variable1, variable2

   integer :: ounit,iunit,namelist_unit

   stderr = 0
   stdout = 6


   write(6,'(a)')' [1] Initialize namelist variables and other scalars.'


   call da_get_unit(ounit)
   call da_get_unit(iunit)
   call da_get_unit(namelist_unit)


   start_date = '2004030312'
   end_date = '2004033112'
   interval = 24
   ne = 1
   variable1 = 'ps_u' 
   variable2 = 'rh'

   open(unit=namelist_unit, file='gen_be_cov_contrib_nl.nl', &
        form='formatted', status='old', action='read')
   read(namelist_unit, gen_be_cov3d_nl)
   close(namelist_unit)

   read(start_date(1:10), fmt='(i10)')sdate
   read(end_date(1:10), fmt='(i10)')edate
   write(6,'(4a)')' Computing covariance for fields ', variable1 , ' and ', variable2
   write(6,'(4a)') ' Time period is ', start_date, ' to ', end_date
   write(6,'(a,i8,a)')' Interval between dates = ', interval, 'hours.'
   write(6,'(a,i8)')' Number of ensemble members at each time = ', ne

   date = start_date
   cdate = sdate

   member = 1
   write(ce,'(i3.3)')member


   filename = trim(variable2)//'/'//date(1:10)
   filename = trim(filename)//'.'//trim(variable2)//'.e'//ce
   open (iunit, file = filename, form='unformatted')
   read(iunit)ni, nj, nk
   close(iunit)


   allocate( bin(1:ni,1:nj,1:nk) )
   allocate( bin2d(1:ni,1:nj) )
   allocate( field1(1:ni,1:nj) )
   allocate( field2(1:ni,1:nj,1:nk) )
   allocate( field2_balanced(1:ni,1:nj,1:nk) )
   allocate( field(1:ni,1:nj) )


   filename = 'bin.data'
   open (iunit, file = filename, form='unformatted')
   read(iunit)bin_type
   read(iunit)lat_min, lat_max, binwidth_lat
   read(iunit)hgt_min, hgt_max, binwidth_hgt
   read(iunit)num_bins, num_bins2d
   read(iunit)bin(1:ni,1:nj,1:nk)
   read(iunit)bin2d(1:ni,1:nj)
   close(iunit)

   allocate  (regcoeff_psi_chi(1:num_bins))
   allocate  (regcoeff_psi_t(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_psi_ps(1:nk,1:num_bins2d))
   allocate  (regcoeff_psi_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_t(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_ps(1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_t_u_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_ps_u_rh(1:nk,1:num_bins2d))

   allocate( regcoeff_field1_field2(1:nk,1:num_bins2d) )

   allocate( bin_pts(1:num_bins) )
   allocate( covar(1:num_bins) )
   allocate( var(1:num_bins) )
   bin_pts(:) = 0
   covar(:) = 0.0
   var(:) = 0.0


   write(6,'(2a)')' [2] Read fields, and calculate correlations'



   filename = 'gen_be_stage2.dat'
   open(iunit, file = filename, form='unformatted',status='old')
   do k = 1, 2
   read (iunit)
   end do

   regcoeff_var='regcoeff_'//trim(variable1)//'_'//trim(variable2)
   do k = 1 , 9
   read (iunit) variable
   select case( trim(adjustl(variable)) )

   case ('regcoeff_psi_chi')
   read (iunit) regcoeff_psi_chi
   case ('regcoeff_psi_t')
   read (iunit) regcoeff_psi_t
   case ('regcoeff_psi_ps')
   read (iunit) regcoeff_psi_ps
   case ('regcoeff_psi_rh')
   read (iunit) regcoeff_psi_rh
   case ('regcoeff_chi_u_t')
   read (iunit) regcoeff_chi_u_t
   case ('regcoeff_chi_u_ps')
   read (iunit) regcoeff_chi_u_ps
   case ('regcoeff_chi_u_rh')
   read (iunit) regcoeff_chi_u_rh
   case ('regcoeff_t_u_rh')
   read (iunit) regcoeff_t_u_rh
   case ('regcoeff_ps_u_rh')
   read (iunit) regcoeff_ps_u_rh
   if(trim(adjustl(variable)) == trim(adjustl(regcoeff_var)) )&
   regcoeff_field1_field2 = regcoeff_ps_u_rh
   exit
   case default;
      write(6,fmt='(a)')'Read problem gen_be_stage2.dat'
      write(6,'(a,a)')' Trying to read regression coefficients:',trim(adjustl(variable))
      stop
   end select
   end do

   close(iunit)
   deallocate (regcoeff_psi_chi)
   deallocate (regcoeff_psi_t)
   deallocate (regcoeff_psi_ps)
   deallocate (regcoeff_psi_rh)
   deallocate (regcoeff_chi_u_t)
   deallocate (regcoeff_chi_u_ps)
   deallocate (regcoeff_chi_u_rh)
   deallocate (regcoeff_t_u_rh)
   deallocate (regcoeff_ps_u_rh)


   do while ( cdate <= edate )

      do member = 1, ne

         write(ce,'(i3.3)')member


         filename = trim(variable1)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable1)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)field1
         close(iunit)


         filename = trim(variable2)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable2)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)field2
         close(iunit)


         do j = 1, nj
            do i = 1, ni
               b = bin2d(i,j)
               do k = 1, nk
                  field2_balanced(i,j,k) = regcoeff_field1_field2(k,b) * field1(i,j)
               end do
            end do
         end do
         

         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  b = bin(i,j,k)
                  bin_pts(b) = bin_pts(b) + 1
                  coeffa = 1.0 / real(bin_pts(b))
                  coeffb = real(bin_pts(b)-1) * coeffa
                  covar(b) = coeffb * covar(b) + coeffa * field2_balanced(i,j,k) * field2(i,j,k)
                  var(b) = coeffb * var(b) + coeffa * field2(i,j,k) * field2(i,j,k)
               end do
            end do
         end do

      end do  


      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate

   end do     

   filename = trim(variable1)//'.'//trim(variable2)//'.dat'
   open (ounit, file = filename, status='unknown')
   do k = 1, nk
      do j = 1, nj
         do i = 1, ni
            b = bin(i,j,k) 
            if ( var(b) /= 0.0 ) then
               write(ounit,'(f16.8)')covar(b) / var(b)
            else
               write(ounit,'(f16.8)')0.0
            end if
         end do
      end do
   end do
   close(ounit)

   filename = trim(variable1)//'.'//trim(variable2)//'.bin'
   open (ounit, file = filename, form = 'unformatted')
   do k = 1, nk
      do j = 1, nj
         do i = 1, ni
            b = bin(i,j,k) 
            if ( var(b) /= 0.0 ) then

               field(i,j) = covar(b) / var(b)
            else

               field(i,j) = 0.
            end if
         end do
      end do
      write(ounit) field
   end do
   close(ounit)


end program gen_be_cov2d3d_contrib

