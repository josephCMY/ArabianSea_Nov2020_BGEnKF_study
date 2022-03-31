












program gen_be_stage1_1dvar














   use da_control, only : filename_len,gaussian_lats
   use da_tools_serial, only : da_get_unit, da_advance_cymdh
   use da_gen_be, only : da_create_bins

   implicit none

   real, parameter     :: t_factor = 1.0             
   real, parameter     :: q_factor = 1000.0          
   real, parameter     :: ps_factor = 0.01           
   real, parameter     :: u_factor = 1.0             

   character*10        :: start_date, end_date       
   character*10        :: date, new_date             
   character*10        :: variable                   
   character*3         :: be_method                  
   character*3         :: ce                         
   character(len=filename_len)        :: dat_dir                    
   character*80        :: expt                       
   character(len=filename_len)        :: filename                   
   integer             :: ni, nj, nk                 
   integer             :: member, b, i, j, k         
   integer             :: sdate, cdate, edate        
   integer             :: interval                   
   integer             :: ne                         
   integer             :: bin_type                   
   integer             :: num_bins                   
   integer             :: num_bins2d                 
   integer             :: bin_pts_total              
   real                :: lat_min, lat_max           
   real                :: binwidth_lat               
   real                :: hgt_min, hgt_max           
   real                :: binwidth_hgt               
   real                :: coeffa, coeffb             
   logical             :: first_time                 
   logical             :: remove_mean                

   real, allocatable   :: t_prime(:,:,:)             
   real, allocatable   :: q_prime(:,:,:)             
   real, allocatable   :: dummy(:,:,:)               
   real, allocatable   :: ps_prime(:,:)              
   real, allocatable   :: t2_prime(:,:)              
   real, allocatable   :: q2_prime(:,:)              
   real, allocatable   :: u10_prime(:,:)             
   real, allocatable   :: v10_prime(:,:)             
   real, allocatable   :: height(:,:,:)              
   real, allocatable   :: latitude(:,:)              
   integer, allocatable:: bin(:,:,:)                 
   integer, allocatable:: bin2d(:,:)                 
   integer, allocatable:: bin_pts(:)                 
   integer, allocatable:: bin_pts2d(:)               
   real, allocatable   :: t_mean(:)                  
   real, allocatable   :: q_mean(:)                  
   real, allocatable   :: ps_mean(:)                 
   real, allocatable   :: t2_mean(:)                 
   real, allocatable   :: q2_mean(:)                 
   real, allocatable   :: u10_mean(:)                
   real, allocatable   :: v10_mean(:)                
   real, allocatable   :: t_ms(:)                   
   real, allocatable   :: q_ms(:)                   
   real, allocatable   :: ps_ms(:)                  
   real, allocatable   :: t2_ms(:)                  
   real, allocatable   :: q2_ms(:)                  
   real, allocatable   :: u10_ms(:)                 
   real, allocatable   :: v10_ms(:)                 

   namelist / gen_be_stage1_nl / start_date, end_date, interval, &
                                 be_method, ne, bin_type, &
                                 lat_min, lat_max, binwidth_lat, &
                                 hgt_min, hgt_max, binwidth_hgt, &
                                 remove_mean, gaussian_lats, expt, dat_dir

  integer :: ounit,iunit,namelist_unit


   write(6,'(a)')' [1] Initialize namelist variables and other scalars.'


   call da_get_unit(ounit)
   call da_get_unit(iunit)
   call da_get_unit(namelist_unit)

   start_date = '2004030312'
   end_date = '2004033112'
   interval = 24
   be_method = 'NMC'
   ne = 1
   bin_type = 5         
   lat_min = -90.0
   lat_max = 90.0
   binwidth_lat = 10.0
   hgt_min = 0.0
   hgt_max = 20000.0
   binwidth_hgt = 1000.0
   remove_mean = .true.
   gaussian_lats = .false.
   expt = 'gen_be_stage1'
   dat_dir = '/data2/hcshin/youn/DIFF63'

   open(unit=namelist_unit, file='gen_be_stage1_nl.nl', &
        form='formatted', status='old', action='read')
   read(namelist_unit, gen_be_stage1_nl)
   close(namelist_unit)

   read(start_date(1:10), fmt='(i10)')sdate
   read(end_date(1:10), fmt='(i10)')edate
   write(6,'(4a)')' Computing statistics for dates ', start_date, ' to ', end_date
   write(6,'(a,i8,a)')' Interval between dates = ', interval, 'hours.'
   write(6,'(a,i8)')' Number of ensemble members at each time = ', ne

   date = start_date
   cdate = sdate
   first_time = .true.


   write(6,'(a)')' [2] Read fields from 1D-Var files, and calculate mean fields'


   do while ( cdate <= edate )
      do member = 1, ne
         if ( member < 10 ) write(ce,'(I1)') member
         if ( member >= 10 .and. member < 100 ) write(ce,'(I2)') member


         write(6,'(a,a)')'    Processing 1D-Var data for date ', date

         filename = 'diff.'//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'_'//date(9:10)//':00:00'
         if ( be_method == 'ENS' ) filename = trim(filename)//'.'//trim(ce)

         open (iunit, file = trim(dat_dir)//'/'//filename, form='unformatted')
         read(iunit)date, ni, nj, nk

         if ( first_time ) then
            write(6,'(a,3i8)')'    i, j, k dimensions are ', ni, nj, nk

            allocate( t_prime(1:ni,1:nj,1:nk) )
            allocate( q_prime(1:ni,1:nj,1:nk) )
            allocate( dummy(1:ni,1:nj,1:nk) )
            allocate( ps_prime(1:ni,1:nj) )
            allocate( t2_prime(1:ni,1:nj) )
            allocate( q2_prime(1:ni,1:nj) )
            allocate( u10_prime(1:ni,1:nj) )
            allocate( v10_prime(1:ni,1:nj) )
            allocate( height(1:ni,1:nj,1:nk) )
            allocate( latitude(1:ni,1:nj) )
            allocate( bin(1:ni,1:nj,1:nk) )
            allocate( bin2d(1:ni,1:nj) )

            allocate( t_mean(1:num_bins) )
            allocate( q_mean(1:num_bins) )
            allocate( ps_mean(1:num_bins2d) )
            allocate( t2_mean(1:num_bins2d) )
            allocate( q2_mean(1:num_bins2d) )
            allocate( u10_mean(1:num_bins2d) )
            allocate( v10_mean(1:num_bins2d) )
            allocate( t_ms(1:num_bins) )
            allocate( q_ms(1:num_bins) )
            allocate( ps_ms(1:num_bins2d) )
            allocate( t2_ms(1:num_bins2d) )
            allocate( q2_ms(1:num_bins2d) )
            allocate( u10_ms(1:num_bins2d) )
            allocate( v10_ms(1:num_bins2d) )
            allocate( bin_pts(1:num_bins) )
            allocate( bin_pts2d(1:num_bins2d) )
            t_mean(:) = 0.0
            q_mean(:) = 0.0
            ps_mean(:) = 0.0
            t2_mean(:) = 0.0
            q2_mean(:) = 0.0
            u10_mean(:) = 0.0
            v10_mean(:) = 0.0
            t_ms(:) = 0.0
            q_ms(:) = 0.0
            ps_ms(:) = 0.0
            t2_ms(:) = 0.0
            q2_ms(:) = 0.0
            u10_ms(:) = 0.0
            v10_ms(:) = 0.0
            bin_pts(:) = 0
            bin_pts2d(:) = 0

         end if

         read(iunit)dummy 
         read(iunit)dummy 
         read(iunit)t_prime 
         read(iunit)dummy 
         read(iunit)ps_prime
         read(iunit)height
         read(iunit)latitude
         close(iunit)

         if ( first_time ) then
            call da_create_bins( ni, nj, nk, bin_type, num_bins, num_bins2d, bin, bin2d, &
                                 lat_min, lat_max, binwidth_lat, &
                                 hgt_min, hgt_max, binwidth_hgt, latitude, height )


            filename = 'bin.data'
            open (ounit, file = filename, form='unformatted')
            write(ounit)bin_type
            write(ounit)lat_min, lat_max, binwidth_lat
            write(ounit)hgt_min, hgt_max, binwidth_hgt
            write(ounit)num_bins, num_bins2d
            write(ounit)bin(1:ni,1:nj,1:nk)
            write(ounit)bin2d(1:ni,1:nj)
            close(ounit)

            first_time = .false.
         end if

         filename = 'q_diff.'//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'_'//date(9:10)//':00:00'
         if ( be_method == 'ENS' ) filename = trim(filename)//'.'//trim(ce)

         open (iunit, file = trim(dat_dir)//'/'//filename, form='unformatted')
         read(iunit)date, ni, nj, nk
         read(iunit)q_prime
         close(iunit)


         filename = 'sfc_var_diff.'//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'_'//date(9:10)//':00:00'
         if ( be_method == 'ENS' ) then
            if ( member < 10 ) write(ce,'(I1)') member
            if ( member >= 10 .and. member < 100 ) write(ce,'(I2)') member
            filename = trim(filename)//'.'//trim(ce)
         endif
         open (iunit, file = trim(dat_dir)//'/'//filename, form='unformatted')
         read(iunit)date, ni, nj
         read(iunit)t2_prime
         read(iunit)q2_prime
         read(iunit)u10_prime
         read(iunit)v10_prime
         close(iunit)


         t_prime = t_prime * t_factor
         q_prime = q_prime * q_factor
         ps_prime = ps_prime * ps_factor
         t2_prime = t2_prime * t_factor
         q2_prime = q2_prime * q_factor
         u10_prime = u10_prime * u_factor
         v10_prime = v10_prime * u_factor






         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  b = bin(i,j,k)
                  bin_pts(b) = bin_pts(b) + 1
                  coeffa = 1.0 / real(bin_pts(b))
                  coeffb = real(bin_pts(b)-1) * coeffa
                  t_mean(b) = coeffb * t_mean(b)  + coeffa * t_prime(i,j,k)
                  q_mean(b) = coeffb * q_mean(b)  + coeffa * q_prime(i,j,k)
                  t_ms(b) = coeffb * t_ms(b)  + coeffa * t_prime(i,j,k) * t_prime(i,j,k)
                  q_ms(b) = coeffb * q_ms(b)  + coeffa * q_prime(i,j,k) * q_prime(i,j,k)
               end do
            end do
         end do


         do j = 1, nj
            do i = 1, ni
               b = bin2d(i,j)
               bin_pts2d(b) = bin_pts2d(b) + 1
               coeffa = 1.0 / real(bin_pts2d(b))
               coeffb = real(bin_pts2d(b)-1) * coeffa
               ps_mean(b) = coeffb * ps_mean(b) + coeffa * ps_prime(i,j)
               t2_mean(b) = coeffb * t2_mean(b) + coeffa * t2_prime(i,j)
               q2_mean(b) = coeffb * q2_mean(b) + coeffa * q2_prime(i,j)
               u10_mean(b) = coeffb * u10_mean(b) + coeffa * u10_prime(i,j)
               v10_mean(b) = coeffb * v10_mean(b) + coeffa * v10_prime(i,j)
               ps_ms(b) = coeffb * ps_ms(b) + coeffa * ps_prime(i,j) * ps_prime(i,j)
               t2_ms(b) = coeffb * t2_ms(b) + coeffa * t2_prime(i,j) * t2_prime(i,j)
               q2_ms(b) = coeffb * q2_ms(b) + coeffa * q2_prime(i,j) * q2_prime(i,j)
               u10_ms(b) = coeffb * u10_ms(b) + coeffa * u10_prime(i,j) * u10_prime(i,j)
               v10_ms(b) = coeffb * v10_ms(b) + coeffa * v10_prime(i,j) * v10_prime(i,j)
            end do
         end do

      end do  


      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate
   end do     

   bin_pts_total = sum(bin_pts(1:num_bins))
   write(6,'(a40,i10,2f18.8)') ' #pts/Mean/ms for model level T: ', bin_pts_total, &
                              sum(t_mean(:)) / ( num_bins * t_factor ), &
                              sqrt(sum(t_ms(:)) / num_bins ) / t_factor
   write(6,'(a40,i10,2f18.8)') ' #pts/Mean/ms for model level q: ', bin_pts_total, &
                              sum(q_mean(:)) / ( num_bins * q_factor ), &
                              sqrt(sum(q_ms(:)) / num_bins ) / q_factor

   bin_pts_total = sum(bin_pts2d(1:num_bins2d))
   write(6,'(a40,i10,2f18.8)') ' #pts/Mean/ms for surface pressure: ', bin_pts_total, &
                              sum(ps_mean(:)) / ( num_bins2d * ps_factor ), &
                              sqrt(sum(ps_ms(:)) / num_bins2d ) / ps_factor
   write(6,'(a40,i10,2f18.8)') ' #pts/Mean/ms/factor for 2m Temperature: ', bin_pts_total, &
                              sum(t2_mean(:)) / ( num_bins2d * t_factor ), &
                              sqrt(sum(t2_ms(:)) / num_bins2d ) / t_factor
   write(6,'(a40,i10,2f18.8)') ' #pts/Mean/ms/factor for 2m Humidity: ', bin_pts_total, &
                              sum(q2_mean(:)) / ( num_bins2d * q_factor ), &
                              sqrt(sum(q2_ms(:)) / num_bins2d ) / q_factor
   write(6,'(a40,i10,2f18.8)') ' #pts/Mean/ms/factor for 10m u-wind: ', bin_pts_total, &
                              sum(u10_mean(:)) / ( num_bins2d * u_factor ), &
                              sqrt(sum(u10_ms(:)) / num_bins2d ) / u_factor
   write(6,'(a40,i10,2f18.8)') ' #pts/Mean/ms/factor for 10m v-wind: ', bin_pts_total, &
                              sum(v10_mean(:)) / ( num_bins2d * u_factor ), &
                              sqrt(sum(v10_ms(:)) / num_bins2d ) / u_factor


   write(6,'(a)')' [2] Read fields again, and remove time/ensemble/area mean'


   date = start_date
   cdate = sdate

   do while ( cdate <= edate )
      do member = 1, ne
         if ( member < 10 ) write(ce,'(I1)') member
         if ( member >= 10 .and. member < 100 ) write(ce,'(I2)') member

         write(6,'(a,a)')'    Removing mean for date ', date

         filename = 'diff.'//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'_'//date(9:10)//':00:00'
         if ( be_method == 'ENS' ) filename = trim(filename)//'.'//trim(ce)
         open (iunit, file = trim(dat_dir)//'/'//filename, form='unformatted')
         read(iunit)date, ni, nj, nk
         read(iunit)dummy
         read(iunit)dummy
         read(iunit)t_prime
         read(iunit)dummy
         read(iunit)ps_prime
         close(iunit)

         filename = 'q_diff.'//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'_'//date(9:10)//':00:00'
         if ( be_method == 'ENS' ) filename = trim(filename)//'.'//trim(ce)
         open (iunit, file = trim(dat_dir)//'/'//filename, form='unformatted')
         read(iunit)date, ni, nj, nk
         read(iunit)q_prime
         close(iunit)

         filename = 'sfc_var_diff.'//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'_'//date(9:10)//':00:00'
         if ( be_method == 'ENS' ) filename = trim(filename)//'.'//trim(ce)
         open (iunit, file = trim(dat_dir)//'/'//filename, form='unformatted')
         read(iunit)date, ni, nj
         read(iunit)t2_prime
         read(iunit)q2_prime
         read(iunit)u10_prime
         read(iunit)v10_prime
         close(iunit)


         t_prime = t_prime * t_factor
         q_prime = q_prime * q_factor
         ps_prime = ps_prime * ps_factor
         t2_prime = t2_prime * t_factor
         q2_prime = q2_prime * q_factor
         u10_prime = u10_prime * u_factor
         v10_prime = v10_prime * u_factor





         if ( remove_mean ) then

            do k = 1, nk
               do j = 1, nj
                  do i = 1, ni
                     b = bin(i,j,k)
                     t_prime(i,j,k) = t_prime(i,j,k) - t_mean(b)
                     q_prime(i,j,k) = q_prime(i,j,k) - q_mean(b)
                  end do
               end do
            end do


            do j = 1, nj
               do i = 1, ni
                  b = bin2d(i,j)
                  ps_prime(i,j) = ps_prime(i,j) - ps_mean(b)
                  t2_prime(i,j) = t2_prime(i,j) - t2_mean(b)
                  q2_prime(i,j) = q2_prime(i,j) - q2_mean(b)
                  u10_prime(i,j) = u10_prime(i,j) - u10_mean(b)
                  v10_prime(i,j) = v10_prime(i,j) - v10_mean(b)
               end do
            end do
         end if





         write(ce,'(i3.3)')member


         variable = 't'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, nk
         write(ounit)t_prime
         close(ounit)


         variable = 'q'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, nk
         write(ounit)q_prime
         close(ounit)


         variable = 'ps' 
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, 1
         write(ounit).true., .false.
         write(ounit)ps_prime
         close(ounit)


         variable = 't2' 
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, 1
         write(ounit).true., .false.
         write(ounit)t2_prime
         close(ounit)


         variable = 'q2' 
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, 1
         write(ounit).true., .false.
         write(ounit)q2_prime
         close(ounit)


         variable = 'u10' 
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, 1
         write(ounit).true., .false.
         write(ounit)u10_prime
         close(ounit)


         variable = 'v10' 
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, 1
         write(ounit).true., .false.
         write(ounit)v10_prime
         close(ounit)

      end do  


      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate
   end do     

   write(6,'(a)')' Unbiased perturbations written out, scaled by the following factors:'
   write(6,'(a,f15.5)')' T scaling factor = ', t_factor
   write(6,'(a,f15.5)')' q scaling factor = ', q_factor
   write(6,'(a,f15.5)')' ps scaling factor = ', ps_factor
   write(6,'(a,f15.5)')' u scaling factor = ', u_factor

end program gen_be_stage1_1dvar

