












program gen_be_hist








   use da_control, only : stderr, stdout, filename_len
   use da_tools_serial, only : da_get_unit,da_advance_cymdh

   implicit none

   character*10        :: start_date, end_date       
   character*10        :: date, new_date             
   character*10        :: variable                   
   character(len=filename_len)        :: dat_dir     
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
   logical             :: first_time                 
   real, allocatable   :: field(:,:,:)               
   integer, allocatable:: bin(:,:,:)                 
   integer, allocatable:: bin2d(:,:)                 
   integer, allocatable:: bin_pts(:)                 
   integer, allocatable:: hist(:,:)                  
   real, allocatable   :: var_field(:)               
   real, allocatable   :: stdev_field(:)               
   real, allocatable   :: class_hist(:,:)            
   integer, allocatable   :: rain_class(:,:)         
   character*10           :: rainclvar               
   integer                :: Nstdev, N_dim_hist      
   integer             :: intcl                      

   namelist / gen_be_hist_nl / start_date, end_date, interval, &
                                ne, variable, Nstdev, N_dim_hist

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
   variable = 'psi'
   dat_dir = '/mmmtmp1/dmbarker'
   Nstdev = 5
   N_dim_hist = 20

   open(unit=namelist_unit, file='gen_be_hist_nl.nl', &
        form='formatted', status='old', action='read')
   read(namelist_unit, gen_be_hist_nl)
   close(namelist_unit)

   read(start_date(1:10), fmt='(i10)')sdate
   read(end_date(1:10), fmt='(i10)')edate
   write(6,'(2a)')' Computing error histogram for field ', variable
   write(6,'(4a)') ' Time period is ', start_date, ' to ', end_date
   write(6,'(a,i8,a)')' Interval between dates = ', interval, 'hours.'
   write(6,'(a,i8)')' Number of ensemble members at each time = ', ne
   write(6,'(2(a,i8))')' Parameters of the histogram Nstdev = ', Nstdev,&
        ' and N_dim_hist = ',N_dim_hist

   date = start_date
   cdate = sdate


   write(6,'(a)')' [2] Accumulate values of errors per bin and vert. level'

   first_time = .true.

   do while ( cdate <= edate )
      write(6,'(a,a)')'    Processing data for date ', date

      do member = 1, ne

         write(ce,'(i3.3)')member


         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
   
        if ( first_time ) then
            write(6,'(a,3i8)')'    i, j, k dimensions are ', ni, nj, nk
            allocate( bin(1:ni,1:nj,1:nk) )
            allocate( bin2d(1:ni,1:nj) )
            allocate( field(1:ni,1:nj,1:nk) )

            filename = 'bin.data'
            open (iunit+1, file = filename, form='unformatted')
            read (iunit+1) bin_type
            read (iunit+1) lat_min, lat_max, binwidth_lat
            read (iunit+1) hgt_min, hgt_max, binwidth_hgt
            read (iunit+1) num_bins, num_bins2d
            read (iunit+1) bin(1:ni,1:nj,1:nk)
            read (iunit+1) bin2d(1:ni,1:nj)
            close(iunit+1)
         end if

         read(iunit)field
         close(iunit)

         if ( first_time ) then
            write(6,'(a)')' Setup of histogram parameters'
            allocate( class_hist(1:nk,1:N_dim_hist) )
            allocate( var_field(1:nk) )
            allocate( stdev_field(1:nk) )

            class_hist  = 0.0
            var_field   = 0.0
            stdev_field = 0.0

            
            do k=1, nk
               var_field(k)   = 1.0/real(ni*nj-1.0)*sum(field(1:ni,1:nj,k)**2) &
                    -real(ni*nj)/real(ni*nj-1.0)*(sum(field(1:ni,1:nj,k))/real(ni*nj))**2
               stdev_field(k) = sqrt(var_field(k))
               write(6,'(a,i4,3(a,e13.5))')' var_field(',k,')=  ', var_field(k),&
                    ' lower=',-Nstdev*stdev_field(k),' upper=',Nstdev*stdev_field(k)
               do i=1,N_dim_hist
                  class_hist(k,i)=-Nstdev*stdev_field(k)+2*Nstdev*stdev_field(k)*real(i-1)/real(N_dim_hist-1)

               end do
            end do

            
            allocate( hist(1:N_dim_hist,1:num_bins))
            hist(:,:) = 0
            
            if (bin_type==7) allocate( rain_class(1:ni,1:nj) )
            
            first_time = .false.
         end if

         
         if (bin_type==7) then
            
            rainclvar = 'raincl'
            filename = trim(rainclvar)//'/'//date(1:10)
            filename = trim(filename)//'.'//trim(rainclvar)//'.e'//ce//'.01'
            open (iunit, file = filename, form='unformatted')
            read(iunit)ni, nj, nkdum
            read(iunit)rain_class
            close(iunit)
  
            
            
            filename = 'bin.data'
            open (iunit+1, file = filename, form='unformatted')
            read(iunit+1)bin_type
            read(iunit+1)lat_min, lat_max, binwidth_lat
            read(iunit+1)hgt_min, hgt_max, binwidth_hgt
            read(iunit+1)num_bins, num_bins2d
            read(iunit+1)bin(1:ni,1:nj,1:nk)
            read(iunit+1)bin2d(1:ni,1:nj)
            close(iunit+1)
            
            do j = 1, nj
               do i = 1, ni
                  bin2d(i,j)=rain_class(i,j)*num_bins2d/4+bin2d(i,j)
               end do
            end do
            
            do k = 1, nk
               do j = 1, nj
                  do i = 1, ni
                     bin(i,j,k)=rain_class(i,j)*num_bins/4+bin(i,j,k)
                  end do
               end do
            end do
         end if
        
          

         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  b = bin(i,j,k)
                  
                  intcl=int(1+0.5*(N_dim_hist-1.0)*(1+field(i,j,k)/(Nstdev*stdev_field(k))))
                  if (intcl.ge.1 .and. intcl.le.N_dim_hist) then
                     hist(intcl,b)=hist(intcl,b)+1
                  else
                     
                     if (abs(field(i,j,k))>50*Nstdev*stdev_field(k)) then
                        write(6,'(3(a,i3),2(a,e20.5e3))') ' WARNING Gross error -> err_field(',i," ",j," ",k,&
                             ") =",field(i,j,k)," whereas stdev =",stdev_field(k)
                     end if
                  end if
               end do
            end do
         end do
      
      end do  


      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate
   end do     


   write(6,'(a)')' [3] Write out computed histogram'


   
   filename = 'hist.'//trim(variable)//'.dat'
   open (ounit, file = filename, form='unformatted')
   write(ounit)nk,num_bins,N_dim_hist
   
   do k=1,nk
      write(ounit)class_hist(k,:)
   end do
   do b=1,num_bins
      write(ounit)hist(:,b)
   end do
   close(ounit)

end program gen_be_hist

