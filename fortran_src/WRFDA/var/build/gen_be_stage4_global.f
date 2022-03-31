












program gen_be_stage4_global

   use da_control, only : stderr, stdout, filename_len, pi, gaussian_lats
   use da_be_spectral, only : da_vv_to_v_spectral, da_initialize_h, &
      da_calc_power
   use da_tools_serial, only : da_get_unit,da_advance_cymdh

   implicit none

   character*10        :: start_date, end_date       
   character*10        :: date, new_date             
   character*10        :: variable                   
   character(len=filename_len)        :: filename                   
   character*2         :: ck                         
   character*3         :: ce                         
   integer             :: ni, nj, nk                 
   integer             :: num_states                 
   integer             :: sdate, cdate, edate        
   integer             :: interval                   
   integer             :: k, member, n               
   integer             :: ne                         
   integer             :: max_wavenumber             
   integer             :: inc                        
   integer             :: lenr                       
   integer             :: lensav                     
   integer             :: lenwrk                     
   integer             :: alp_size                   
   integer             :: r_cvsize                   

   logical             :: first_time                 
   logical             :: testing_spectral           
   real                :: pi_over_180                
   real                :: coeffa, coeffb             
   real                :: variance                   

   real, allocatable   :: field(:,:)                 

   real, allocatable   :: lat(:)                     
   real, allocatable   :: sinlat(:)                  
   real, allocatable   :: coslat(:)                  
   real, allocatable   :: int_wgts(:)                
   real, allocatable   :: lon(:)                     
   real, allocatable   :: sinlon(:)                  
   real, allocatable   :: coslon(:)                  
   real, allocatable   :: wsave(:)                   
   real, allocatable   :: alp(:)                     
   real, allocatable   :: power(:)                   
   real, allocatable   :: total_power(:)             

   real, allocatable   :: rcv(:)                     

   namelist / gen_be_stage4_global_nl / start_date, end_date, interval, variable, gaussian_lats, &
                                        testing_spectral, ne, k

   integer :: ounit,iunit,namelist_unit

   stderr = 0
   stdout = 6

   call da_get_unit(ounit)
   call da_get_unit(iunit)
   call da_get_unit(namelist_unit)

   pi_over_180 = pi / 180.0

   write(unit=6,fmt='(a/)') "[1] Read in 2D perturbation fields for variable , variable"

   start_date = '2004030312'
   end_date = '2004033112'
   interval = 24
   variable = 'psi'
   gaussian_lats = .false.
   ne = 1
   k = 1

   open(unit=namelist_unit, file='gen_be_stage4_global_nl.nl', &
        form='formatted', status='old', action='read')
   read(namelist_unit, gen_be_stage4_global_nl)
   close(namelist_unit)

   write(ck,'(i2)')k
   if ( k < 10 ) ck = '0'//ck(2:2)

   read(start_date(1:10), fmt='(i10)')sdate
   read(end_date(1:10), fmt='(i10)')edate
   write(unit=6,fmt='(4a)') 'Computing horizontal power spectra for period' , start_date, 'to' , end_date
   write(unit=6,fmt='(a,i8,a)') 'Interval between dates =' , interval,' hours'
   write(unit=6,fmt='(a,i8)') 'Number of ensemble members at each time =', ne
   write(unit=6,fmt='(3a,i4)') 'Variable' , variable,' at level' , k

   date = start_date
   cdate = sdate
   first_time = .true.
   num_states = 1

   do while ( cdate <= edate )
      do member = 1, ne
         write(ce,'(i3.3)')member

         
         

         

         filename = trim(variable)//'/'//date(1:10)//'.'//trim(variable)
         filename = trim(filename)//'.e'//ce//'.'//ck
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk 

         if ( first_time ) then
            write(6,'(a,3i8)')'    i, j, k dimensions are ', ni, nj, nk
            allocate( field(1:ni,1:nj) )
         end if
         read(iunit)field
         close(iunit)

         if ( first_time ) then

            write(unit=6,fmt='(a)') "Initialize spectral transforms"

            inc = 1
            lenr = inc * (ni - 1 ) + 1
            lensav = ni + int(log(real(ni))) + 4
            lenwrk = ni

            allocate( lat(1:nj) )
            allocate( sinlat(1:nj) )
            allocate( coslat(1:nj) )
            allocate( int_wgts(1:nj) )
            allocate( lon(1:ni) )
            allocate( sinlon(1:ni) )
            allocate( coslon(1:ni) )
            allocate( wsave(1:lensav) )

            max_wavenumber =  ni / 2 - 1
            allocate ( power(0:max_wavenumber) )
            allocate ( total_power(0:max_wavenumber) )
            total_power(:) = 0.0

            alp_size = ( nj + 1 ) * ( max_wavenumber + 1 ) * ( max_wavenumber + 2 ) / 4
            allocate( alp( 1:alp_size) )

            call da_initialize_h( ni, nj, max_wavenumber, lensav, alp_size, &
                                  wsave, lon, sinlon, coslon, lat, sinlat, coslat, &
                                  int_wgts, alp )

            r_cvsize = ( max_wavenumber + 1 ) * ( max_wavenumber + 2 )
            allocate( rcv( 1:r_cvsize) )

            
            
            
            
            
            first_time = .false.
         end if

         write(unit=6,fmt='(a)') "Perform gridpoint to spectral decomposition"

         call da_vv_to_v_spectral( ni, nj, max_wavenumber, inc, lenr, lensav, lenwrk, &
                                   alp_size, r_cvsize, alp, wsave, int_wgts, rcv, field )

         write(unit=6,fmt='(a)') "Calculate power spectra"

         call da_calc_power( max_wavenumber, r_cvsize, rcv, power )

         coeffa = 1.0 / real(num_states)
         coeffb = real(num_states-1) * coeffa

         do n = 0, max_wavenumber
            total_power(n) = total_power(n) * coeffb + power(n) * coeffa
            
         end do

         num_states = num_states + 1
      end do  

      
      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate

   end do     

   variance = sum( total_power(0:max_wavenumber) )   
   write(6,'(3a,i2,a,1pe15.5)')' Variable = ', trim(variable), ', Vertical Index = ', &
                                k, ', Variance = ', variance

   filename = trim(variable)//'/'//trim(variable)
   filename = trim(filename)//'.'//ck//'.spectrum'
   open (ounit, file = filename, form='unformatted')
   write(ounit)variable
   write(ounit)max_wavenumber, k
   write(ounit)total_power

end program gen_be_stage4_global
