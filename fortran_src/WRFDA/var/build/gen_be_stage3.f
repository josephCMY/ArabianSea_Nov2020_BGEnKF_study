












program gen_be_stage3

   use da_control, only : stderr, stdout, filename_len, vertical_ip
   use da_gen_be, only : da_eof_decomposition_test, da_eof_decomposition, &
      da_transform_vptovv, da_create_bins
   use da_tools_serial, only : da_get_unit, da_advance_cymdh

   implicit none

   character*10        :: start_date, end_date       
   character*10        :: date, new_date             
   character*10        :: variable                   
   character(len=filename_len)        :: filename                   
   character*2         :: ck                         
   character*3         :: ce                         
   integer             :: ni, nj, nk                 
   integer             :: sdate, cdate, edate        
   integer             :: interval                   
   integer             :: ne                         
   integer             :: i, j, k, k1, k2, b, member 
   integer             :: bin_type                   
   integer             :: num_bins                   
   integer             :: num_bins2d                 
   real                :: lat_min, lat_max           
   real                :: binwidth_lat               
   real                :: hgt_min, hgt_max           
   real                :: binwidth_hgt               
   real                :: inv_nij                    
   real                :: mean_field                 
   real                :: coeffa, coeffb             
   logical             :: first_time                 
   logical             :: testing_eofs               
   logical             :: use_global_eofs            
   logical             :: data_on_levels             
   logical             :: twod_field                 

   integer, allocatable:: bin(:,:,:)                 
   integer, allocatable:: bin2d(:,:)                 
   integer, allocatable:: bin_pts2d(:)               
   real, allocatable   :: latitude(:,:)              
   real, allocatable   :: height(:,:,:)              
   real, allocatable   :: field(:,:,:)               
   real, allocatable   :: field_out(:,:,:)          
   real, allocatable   :: vertical_wgt(:,:,:)        
   real, allocatable   :: bv(:,:,:)                  
   real, allocatable   :: work(:,:)                  
   real, allocatable   :: e_vec_loc(:,:,:)           
   real, allocatable   :: e_val_loc(:,:)             
   real*8, allocatable :: e_vec(:,:)                 
   real*8, allocatable :: e_val(:)                   
   real*8, allocatable :: evec(:,:,:)                
   real*8, allocatable :: eval(:,:)                  

   namelist / gen_be_stage3_nl / start_date, end_date, interval, variable, &
                                 ne, bin_type, &
                                 lat_min, lat_max, binwidth_lat, &
                                 hgt_min, hgt_max, binwidth_hgt, &
                                 testing_eofs, use_global_eofs, data_on_levels

   integer :: ounit,iunit,namelist_unit

   stderr = 0
   stdout = 6

   write(6,'(a)')' [1] Initialize namelist variables and other scalars.'

   call da_get_unit(ounit)
   call da_get_unit(iunit)
   call da_get_unit(namelist_unit)


   vertical_ip = 0
   twod_field = .false.

   start_date = '2004030312'
   end_date = '2004033112'
   interval = 24
   variable = 'psi'
   ne = 1
   bin_type = 5
   lat_min = -90.0
   lat_max = 90.0
   binwidth_lat = 10.0
   hgt_min = 0.0
   hgt_max = 20000.0
   binwidth_hgt = 1000.0
   testing_eofs = .true.
   use_global_eofs = .true.
   data_on_levels = .false.

   open(unit=namelist_unit, file='gen_be_stage3_nl.nl', &
        form='formatted', status='old', action='read')
   read(namelist_unit, gen_be_stage3_nl)
   close(namelist_unit)

   read(start_date(1:10), fmt='(i10)')sdate
   read(end_date(1:10), fmt='(i10)')edate
   write(6,'(4a)')' Computing vertical error statistics for dates ', start_date, ' to ', end_date
   write(6,'(a,i8,a)')' Interval between dates = ', interval, 'hours.'
   write(6,'(a,i8)')' Number of ensemble members at each time = ', ne

   date = start_date
   cdate = sdate
   first_time = .true.

   if ( .not. data_on_levels ) then 

      write(6,'(2a)')' [2] Process fields for variable ', variable

      do while ( cdate <= edate )
         do member = 1, ne
            write(6,'(5a,i4)')'    Processing data for date ', date, ', variable ', trim(variable), &
                              ' and member ', member

            write(ce,'(i3.3)')member

            
            filename = 'fullflds'//'/'//date(1:10)//'.'//'fullflds'//'.e'//ce
            open (iunit, file = filename, form='unformatted')

            read(iunit)ni, nj, nk
            if ( first_time ) then
               write(6,'(a,3i8)')'    i, j, k dimensions are ', ni, nj, nk
               allocate( latitude(1:ni,1:nj) )
               allocate( height(1:ni,1:nj,1:nk) )
               allocate( bin(1:ni,1:nj,1:nk) )
               allocate( bin2d(1:ni,1:nj) )
            end if
            read(iunit)latitude
            read(iunit)height

            close(iunit)

            if ( first_time ) then
               
               call da_create_bins( ni, nj, nk, bin_type, num_bins, num_bins2d, bin, bin2d, &
                                    lat_min, lat_max, binwidth_lat, &
                                    hgt_min, hgt_max, binwidth_hgt, latitude, height )

               allocate( bin_pts2d(1:num_bins2d) )
               bin_pts2d(:) = 0
            end if

            filename = trim(variable)//'/'//date(1:10)
            filename = trim(filename)//'.'//trim(variable)//'.e'//ce
            if ( trim(variable) == 'ps_u' .or. trim(variable) == 'ps' ) then 
               twod_field = .true.
               filename = trim(filename)//'.01'
            end if

            open (iunit, file = filename, form='unformatted')
            read(iunit)ni, nj, nk

            if ( first_time ) then
               inv_nij = 1.0 / real(ni*nj)
               allocate( field(1:ni,1:nj,1:nk) )
               allocate( bv(1:nk,1:nk,1:num_bins2d) )
               bv(:,:,:) = 0.0
               first_time = .false.
            end if

            read(iunit)field
            close(iunit)

            
            do k = 1, nk
               mean_field = sum(field(1:ni,1:nj,k)) * inv_nij
               field(1:ni,1:nj,k) = field(1:ni,1:nj,k) - mean_field
            end do

            do j = 1, nj
               do i = 1, ni
                  b = bin2d(i,j)
                  bin_pts2d(b) = bin_pts2d(b) + 1
                  coeffa = 1.0 / real(bin_pts2d(b))
                  coeffb = real(bin_pts2d(b)-1) * coeffa
                  do k1 = 1, nk
                     do k2 = 1, k1
                        bv(k1,k2,b) = coeffb * bv(k1,k2,b) + coeffa * field(i,j,k1) * field(i,j,k2)
                     end do
                  end do
               end do
            end do
         end do  

         
         call da_advance_cymdh( date, interval, new_date )
         date = new_date
         read(date(1:10), fmt='(i10)')cdate
      end do     

      

      do b = 1, num_bins2d
         do k1 = 1, nk
            do k2 = k1+1, nk 
               bv(k1,k2,b) = bv(k2,k1,b)
            end do
         end do
      end do

      write(6,'(2a)')' [3] Calculate eigenvectors and eigenvalues for variable ', variable

      allocate( work(1:nk,1:nk) )
      allocate( e_vec_loc(1:nk,1:nk,1:num_bins2d) )
      allocate( e_val_loc(1:nk,1:num_bins2d) )
      allocate( e_vec(1:nk,1:nk) )
      allocate( e_val(1:nk) )

      
      do b = 1, num_bins2d
         write(6,'(2(a,i6))')' Calculate eigenvectors and eigenvalues for bin ', b, &
                             ' of ', num_bins2d

         work(1:nk,1:nk) = bv(1:nk,1:nk,b)
         call da_eof_decomposition( nk, work, e_vec, e_val )
         e_vec_loc(1:nk,1:nk,b) = e_vec(1:nk,1:nk)
         e_val_loc(1:nk,b) = e_val(1:nk)
      end do

      
      work(1:nk,1:nk) = 0.0
      do b = 1, num_bins2d
         work(1:nk,1:nk) = work(1:nk,1:nk) + bv(1:nk,1:nk,b)
      end do
      work(1:nk,1:nk) = work(1:nk,1:nk) / real( num_bins2d )

      call da_eof_decomposition( nk, work, e_vec, e_val )

      if ( testing_eofs ) then
         call da_eof_decomposition_test( nk, work, e_vec, e_val )
      end if

      
      filename = 'gen_be_stage3.'//trim(variable)//'.dat'
      open (ounit, file = filename, form='unformatted')
      write(ounit)variable
      write(ounit)nk, num_bins2d
      write(ounit)e_vec
      write(ounit)e_val
      write(ounit)e_vec_loc
      write(ounit)e_val_loc
      close(ounit)

      
      if ( use_global_eofs ) then
         do b = 1, num_bins2d
            e_vec_loc(1:nk,1:nk,b) = e_vec(1:nk,1:nk)
            e_val_loc(1:nk,b) = e_val(1:nk)
         end do
      end if

      
      allocate( evec(1:nj,1:nk,1:nk) )
      allocate( eval(1:nj,1:nk) )

      do j = 1, nj
         do i = 1, ni
            b = bin2d(i,j)
            evec(j,1:nk,1:nk) = e_vec_loc(1:nk,1:nk,b)
            eval(j,1:nk) = sqrt(e_val_loc(1:nk,b))
         end do
      end do

   end if 

   write(6,'(2a)')' [4] Transform perturbations (or not), and output.'

   date = start_date
   cdate = sdate
   first_time = .true.

   if ( .not. twod_field ) then
      do while ( cdate <= edate )
         do member = 1, ne
            write(6,'(5a,i4)')'    Date = ', date, ', variable ', trim(variable), &
                              ' and member ', member

            write(ce,'(i3.3)')member

            filename = trim(variable)//'/'//date(1:10)
            filename = trim(filename)//'.'//trim(variable)//'.e'//ce

            open (iunit, file = filename, form='unformatted')
            read(iunit)ni, nj, nk

            if ( first_time ) then
               if ( data_on_levels) allocate( latitude(1:ni,1:nj) )   
               if ( data_on_levels) allocate( field(1:ni,1:nj,1:nk) ) 
               allocate( field_out(1:ni,1:nj,1:nk) )
               allocate( vertical_wgt(1:ni,1:nj,1:nk) )
               vertical_wgt(1:ni,1:nj,1:nk) = 1.0 
               first_time = .false.
            end if

            read(iunit)field
            close(iunit)

            if ( data_on_levels ) then
               
               field_out(:,:,:) = field(:,:,:)
            else
               
               call da_transform_vptovv( evec, eval, vertical_wgt, &
                  field, field_out, nk, &
                  1, nk, & 
                  1, ni, 1, nj, 1, nk, & 
                  1, ni, 1, nj, 1, nk )  
            end if

            

            do k = 1, nk
               write(ck,'(i2)')k
               if ( k < 10 ) ck = '0'//ck(2:2)
               
               filename = trim(variable)//'/'//date(1:10)//'.'//trim(variable)
               filename = trim(filename)//'.e'//ce//'.'//ck
               open (ounit, file = filename, form='unformatted')
               write(ounit)ni, nj, k
               write(ounit)field_out(1:ni,1:nj,k)
               close(ounit)
            end do
         end do 

         
         call da_advance_cymdh( date, interval, new_date )
         date = new_date
         read(date(1:10), fmt='(i10)')cdate
      end do
   end if

end program gen_be_stage3

