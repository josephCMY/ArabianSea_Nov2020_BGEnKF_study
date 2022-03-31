












program gen_mbe_stage2

















   use da_control, only : stdout, stderr, filename_len, cv_options
   use da_tools_serial, only : da_get_unit, da_advance_cymdh
   use da_gen_be, only : da_filter_regcoeffs

   implicit none
   
   
   real, parameter     :: variance_threshold = 1e-6  

   integer             :: count
   integer             :: ni, nj, nk, nkdum          
   integer             :: i, j, k, k1, k2, k3, m, member 
   integer             :: b                          
   integer             :: sdate, cdate, edate        
   integer             :: interval                   
   integer             :: ne                         
   integer             :: bin_type                   
   integer             :: mmax                       
   integer             :: num_bins                   
   integer             :: num_bins2d                 
   integer             :: ounit,iunit,namelist_unit  
   integer             :: num_passes                 

   
   real                :: count_inv

   real                :: summ, summ1, summ2         
   real                :: lat_min, lat_max           
   real                :: binwidth_lat               
   real                :: hgt_min, hgt_max           
   real                :: binwidth_hgt               
   real                :: rf_scale                   

   character*10        :: start_date, end_date       
   character*10        :: date, new_date             
   character*80        :: variable                   
   character*3         :: ce                         
   character(len=filename_len) :: filename           

   logical             :: testing_eofs               


   integer, allocatable:: bin(:,:,:)                 
   integer, allocatable:: bin2d(:,:)                 
   integer, allocatable:: bin_pts(:)                 
   integer, allocatable:: bin_pts2d(:)               

   real, allocatable   :: psi(:,:,:)                 
   real, allocatable   :: chi(:,:,:)                 
   real, allocatable   :: temp(:,:,:)                
   real, allocatable   :: rh(:,:,:)                  
   real, allocatable   :: ps(:,:)                    

   real, allocatable   :: chi_u_mean(:,:,:)            
   real, allocatable   :: temp_u_mean(:,:,:)           
   real, allocatable   :: rh_u_mean(:,:,:)             
   real, allocatable   :: ps_u_mean(:,:)               

   real, allocatable   :: var_psi_bin3d(:)               
   real, allocatable   :: var_ps_u(:)                
   real, allocatable   :: var_psi(:,:,:)             
   real, allocatable   :: var_temp(:,:,:)            
   real, allocatable   :: var_chi_u(:,:,:)           
   real, allocatable   :: var_temp_u(:,:,:)          

   real, allocatable   :: var_psi_inv(:,:,:)         
   real, allocatable   :: var_temp_inv(:,:,:)        
   real, allocatable   :: var_chi_u_inv(:,:,:)       
   real, allocatable   :: var_temp_u_inv(:,:,:)      

   real, allocatable   :: covar_psi_chi(:)           
   real, allocatable   :: covar_psi_ps(:,:)          
   real, allocatable   :: covar_psi_temp(:,:,:)      
   real, allocatable   :: covar_psi_rh(:,:,:)        
   real, allocatable   :: covar_temp_rh(:,:,:)       

   real, allocatable   :: covar_chi_u_ps(:,:)        
   real, allocatable   :: covar_chi_u_temp(:,:,:)    
   real, allocatable   :: covar_ps_u_rh(:,:)         
   real, allocatable   :: covar_chi_u_rh(:,:,:)      
   real, allocatable   :: covar_temp_u_rh(:,:,:)     
   
   real, allocatable   :: regcoeff_psi_chi(:)        
   real, allocatable   :: regcoeff_psi_ps(:,:)       
   real, allocatable   :: regcoeff_psi_t(:,:,:)      
   real, allocatable   :: regcoeff_psi_rh(:,:,:)     

   real, allocatable   :: regcoeff_chi_u_ps(:,:)     
   real, allocatable   :: regcoeff_chi_u_t(:,:,:)    
   real, allocatable   :: regcoeff_ps_u_rh(:,:)      
   real, allocatable   :: regcoeff_chi_u_rh(:,:,:)   
   real, allocatable   :: regcoeff_t_u_rh(:,:,:)     

   namelist / gen_be_stage2_nl / start_date, end_date, interval, &
                                 ne, testing_eofs, num_passes, rf_scale, cv_options  
    
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
   testing_eofs = .true.
   num_passes = 0
   rf_scale = 1.0
   cv_options = 5

   open(unit=namelist_unit, file='gen_be_stage2_nl.nl', &
        form='formatted', status='old', action='read')
   read(namelist_unit, gen_be_stage2_nl)
   close(namelist_unit)

   read(start_date(1:10), fmt='(i10)')sdate
   read(end_date(1:10), fmt='(i10)')edate
   write(6,'(4a)')' Computing regression coefficients'
   write(6,'(4a)') ' Time period is ', start_date, ' to ', end_date
   write(6,'(a,i8,a)')' Interval between dates = ', interval, 'hours.'
   write(6,'(a,i8)')' Number of ensemble members at each time = ', ne

   date = start_date
   cdate = sdate
   count = 0
   member = 1


   write(6,'(2a)')' [2] Read fields, and calculate correlations with psi and temp'




   write(ce, '(i3.3)') member

   variable = 'fullflds'
   filename = trim(variable)//'/'//date(1:10)
   filename = trim(filename)//'.'//trim(variable)//'.e'//ce
   open (iunit, file = filename, form='unformatted')

   read(iunit)ni, nj, nk
      allocate( bin(1:ni,1:nj,1:nk) )
      allocate( bin2d(1:ni,1:nj) )
      allocate( psi(1:ni,1:nj,1:nk) )
      allocate( chi(1:ni,1:nj,1:nk) )
      allocate( temp(1:ni,1:nj,1:nk) )
      allocate( rh(1:ni,1:nj,1:nk) )
      allocate( ps(1:ni,1:nj) )

      allocate( chi_u_mean(1:ni,1:nj,1:nk) )
      allocate( temp_u_mean(1:ni,1:nj,1:nk) )
      allocate( rh_u_mean(1:ni,1:nj,1:nk) )
      allocate( ps_u_mean(1:ni,1:nj) )




   close(iunit)



   filename = 'bin.data'
   open (iunit, file = filename, form='unformatted')
   read(iunit)bin_type
   read(iunit)lat_min, lat_max, binwidth_lat
   read(iunit)hgt_min, hgt_max, binwidth_hgt
   read(iunit)num_bins, num_bins2d
   read(iunit)bin(1:ni,1:nj,1:nk)
   read(iunit)bin2d(1:ni,1:nj)
   close(iunit)

   allocate( bin_pts(1:num_bins) )
   allocate( bin_pts2d(1:num_bins2d) )

   allocate( var_psi_bin3d(1:num_bins) )
   allocate( var_ps_u(1:num_bins) )
   allocate( var_psi(1:nk,1:nk,1:num_bins2d) )
   allocate( var_temp(1:nk,1:nk,1:num_bins2d) )

   allocate( var_chi_u(1:nk,1:nk,1:num_bins2d) )
   allocate( var_temp_u(1:nk,1:nk,1:num_bins2d) )

   allocate( var_psi_inv(1:nk,1:nk,1:num_bins2d) )
   allocate( var_temp_inv(1:nk,1:nk,1:num_bins2d) )
   allocate( var_chi_u_inv(1:nk,1:nk,1:num_bins2d) )
   allocate( var_temp_u_inv(1:nk,1:nk,1:num_bins2d) )

   allocate( covar_psi_chi(1:num_bins) )
   allocate( covar_psi_ps(1:nk,1:num_bins2d) )
   allocate( covar_psi_temp(1:nk,1:nk,1:num_bins2d) )
   allocate( covar_psi_rh(1:nk,1:nk,1:num_bins2d) )

   allocate( covar_temp_rh(1:nk,1:nk,1:num_bins2d) )

   allocate( covar_chi_u_ps(1:nk,1:num_bins2d) )
   allocate( covar_chi_u_temp(1:nk,1:nk,1:num_bins2d) )
   allocate( covar_chi_u_rh(1:nk,1:nk,1:num_bins2d) )

   allocate( covar_ps_u_rh(1:nk,1:num_bins2d) )
   allocate( covar_temp_u_rh(1:nk,1:nk,1:num_bins2d) )

   allocate( regcoeff_psi_chi(1:num_bins) )
   allocate( regcoeff_psi_ps(1:nk,1:num_bins2d) )
   allocate( regcoeff_psi_t(1:nk,1:nk,1:num_bins2d) )
   allocate( regcoeff_psi_rh(1:nk,1:nk,1:num_bins2d) )

   allocate( regcoeff_chi_u_ps(1:nk,1:num_bins2d) )
   allocate( regcoeff_chi_u_t(1:nk,1:nk,1:num_bins2d) )
   allocate( regcoeff_chi_u_rh(1:nk,1:nk,1:num_bins2d) )

   allocate( regcoeff_ps_u_rh(1:nk,1:num_bins2d) )
   allocate( regcoeff_t_u_rh(1:nk,1:nk,1:num_bins2d) )



   
   bin_pts(:) = 0
   bin_pts2d(:) = 0

   var_psi_bin3d(:) = 0.


   chi_u_mean(:,:,:) = 0.
   temp_u_mean(:,:,:) = 0.
   ps_u_mean(:,:) = 0.
   rh_u_mean(:,:,:) = 0.


   var_ps_u(:) = 0.
   var_psi(:,:,:) = 0.
   var_temp(:,:,:) = 0.
   var_chi_u(:,:,:) = 0.
   var_temp_u(:,:,:) = 0.

   var_psi_inv(:,:,:) = 0.
   var_temp_inv(:,:,:) = 0.
   var_chi_u_inv(:,:,:) = 0.
   var_temp_u_inv(:,:,:) = 0.

   covar_psi_chi(:) = 0.
   covar_psi_ps(:,:) = 0.
   covar_psi_temp(:,:,:) = 0.
   covar_psi_rh(:,:,:) = 0.

   covar_chi_u_ps(:,:) = 0. 
   covar_chi_u_temp(:,:,:) = 0.
   covar_chi_u_rh(:,:,:) = 0.

   covar_ps_u_rh(:,:) = 0.
   covar_temp_u_rh(:,:,:) = 0.

   regcoeff_psi_chi(:) = 0.
   regcoeff_psi_ps(:,:) = 0.
   regcoeff_psi_t(:,:,:) = 0.
   regcoeff_psi_rh(:,:,:) = 0.

   regcoeff_chi_u_ps(:,:) = 0.
   regcoeff_chi_u_t(:,:,:) = 0.
   regcoeff_ps_u_rh(:,:) = 0.

   regcoeff_chi_u_rh(:,:,:) = 0.
   regcoeff_t_u_rh(:,:,:) = 0.


   

   do while ( cdate <= edate )
   

      do member = 1, ne

         write(ce,'(i3.3)')member

         
         
         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  b = bin(i,j,k)
                  bin_pts(b) = bin_pts(b) + 1
               end do
            end do
         end do


         
         do j = 1, nj
            do i = 1, ni
               b = bin2d(i,j)
               bin_pts2d(b) = bin_pts2d(b) + 1
            end do
         end do
   

         variable = 'psi'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)psi
         close(iunit)
   

         variable = 'chi'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)chi
         close(iunit)
            

         variable = 't'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)temp
         close(iunit)


         variable = 'rh'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)rh     
         close(iunit)


         variable = 'ps'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)ps
         close(iunit)





         call compute_var_bin3d (ni, nj, nk, bin, num_bins, bin_pts, psi, var_psi_bin3d)
         call compute_covar_bin3d (ni, nj, nk, bin, num_bins, bin_pts, chi, psi, covar_psi_chi)         


         call compute_covar2d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, ps, psi, covar_psi_ps)


         call compute_covar3d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, temp, psi, covar_psi_temp)


         call compute_covar3d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, rh, psi, covar_psi_rh)





         call compute_var3d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, psi, var_psi)


         call compute_var3d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, temp, var_temp)
         
      end do 
      

      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate

   end do     

   do b = 1, num_bins
      var_psi_bin3d(b) = var_psi_bin3d(b)/real(bin_pts(b))
      covar_psi_chi(b) = covar_psi_chi(b)/real(bin_pts(b))
   end do

   do b = 1, num_bins2d
      do k1 = 1, nk    
	 do k2 = 1, k1
            
            var_psi(k1,k2,b) = var_psi(k1,k2,b)/real(bin_pts2d(b))
            var_temp(k1,k2,b) = var_temp(k1,k2,b)/real(bin_pts2d(b))
            
	    covar_psi_temp(k1,k2,b) = covar_psi_temp(k1,k2,b)/real(bin_pts2d(b))
	    covar_psi_rh(k1,k2,b) = covar_psi_rh(k1,k2,b)/real(bin_pts2d(b))
	    covar_temp_rh(k1,k2,b) = covar_temp_rh(k1,k2,b)/real(bin_pts2d(b))
         end do
	 do k2 = k1+1, nk
            
	    covar_psi_temp(k1,k2,b) = covar_psi_temp(k1,k2,b)/real(bin_pts2d(b))
	    covar_psi_rh(k1,k2,b) = covar_psi_rh(k1,k2,b)/real(bin_pts2d(b))
	    covar_temp_rh(k1,k2,b) = covar_temp_rh(k1,k2,b)/real(bin_pts2d(b))
         end do
         covar_psi_ps(k1,b) = covar_psi_ps(k1,b)/real(bin_pts2d(b))
      end do
   end do

   do b = 1, num_bins2d
      do k1 = 1, nk    
	 do k2 = k1+1, nk
            
            var_psi(k1,k2,b) = var_psi(k2,k1,b)
            var_temp(k1,k2,b) = var_temp(k2,k1,b) 
         end do
      end do
   end do


   write(6,'(2a)')' [3] Calculate eigenvectors, eigenvalues and inverse for psi/psi covariance '


   mmax = nk
   variable = 'psi'

   call compute_var_inv(variance_threshold, nk, num_bins2d, mmax, var_psi, var_psi_inv, testing_eofs, variable)

   mmax = nk
   variable = 'temp'
   call compute_var_inv(variance_threshold, nk, num_bins2d, mmax, var_temp, var_temp_inv, testing_eofs, variable)


   write(6,'(2a)')' [4] Calculate regression coefficients with full psi and T'



   do b = 1, num_bins
      regcoeff_psi_chi(b) = covar_psi_chi(b) / var_psi_bin3d(b)
   end do


   do b = 1, num_bins2d
      do k1 = 1, nk
         summ = 0.0
         do k2 = 1, nk
            summ = summ + covar_psi_ps(k2,b) * var_psi_inv(k2,k1,b)
         end do
         regcoeff_psi_ps(k1,b) = summ
      end do
   end do


   do b = 1, num_bins2d
      do k1 = 1, nk
         do k2 = 1, nk
            summ1 = 0.0
            summ2 = 0.0
            do k3 = 1, nk
               summ1 = summ1 + covar_psi_temp(k1,k3,b) * var_psi_inv(k3,k2,b)
               summ2 = summ2 + covar_psi_rh(k1,k3,b) * var_psi_inv(k3,k2,b)
            end do
            regcoeff_psi_t(k1,k2,b) = summ1
            regcoeff_psi_rh(k1,k2,b) = summ2
         end do
      end do
   end do


   if ( num_passes > 0 ) then


      write(6,'(a,i4,a)')' [5] Apply ', num_passes, ' pass recursive filter to regression coefficients:'

      call da_filter_regcoeffs( ni, nj, nk, num_bins, num_bins2d, num_passes, rf_scale, bin, &
                                regcoeff_psi_chi, regcoeff_psi_ps, regcoeff_psi_t )
   else
      write(6,'(a)')' [5] num_passes = 0. Bypassing recursive filtering.'
   end if



   write(6,'(a)')' [6] Read standard fields, and compute control variable fields:'


   date = start_date
   cdate = sdate
   count = 0


   do while ( cdate <= edate)


      do member = 1, ne

         count = count + 1
         count_inv = 1.0 / real(count)

         write(ce,'(i3.3)')member



         variable = 'psi'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)psi
         close(iunit)

         variable = 'chi'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)chi
         close(iunit)

         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  b = bin(i,j,k)
                  chi(i,j,k) = chi(i,j,k) - regcoeff_psi_chi(b) * psi(i,j,k)
               end do
            end do
         end do

         
         chi_u_mean = ( real( count-1) * chi_u_mean + chi) * count_inv

         variable = 'chi_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, nk
         write(ounit)chi
         close(ounit)

      end do 


      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate
      
   end do  


   write(6,'(a)')' [6.5] Subtract the mean from the control variable'


   
   do while ( cdate <= edate )

      
      do member = 1, ne

         write(ce,'(i3.3)')member

         
         variable = 'chi_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)chi
         close(iunit)

         chi = chi - chi_u_mean
         
         
         variable = 'chi_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, nk
         write(ounit)chi
         close(ounit)

      end do
      
      
      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate

   end do


   write(6,'(2a)')' [7] Calculate correlations with unbalanced chi'


   date = start_date
   cdate = sdate


   do while ( cdate <= edate )
   

      do member = 1, ne

         write(ce,'(i3.3)')member
         

         variable = 'chi_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)chi
         close(iunit)


         variable = 't'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)temp
         close(iunit)


         variable = 'rh'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)rh     
         close(iunit)


         variable = 'ps'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)ps
         close(iunit)





         call compute_covar3d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, temp, chi, covar_chi_u_temp)


         call compute_covar2d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, ps, chi, covar_chi_u_ps)


         call compute_covar3d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, rh, chi, covar_chi_u_rh)


         call compute_var3d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, chi, var_chi_u)
          
      end do  


      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate
      
   end do   

   do b = 1, num_bins2d
      do k = 1, nk    
         covar_chi_u_ps(k,b) = covar_chi_u_ps(k,b)/real(bin_pts2d(b))
      end do
   end do

   do b = 1, num_bins2d
      do k1 = 1, nk   
	 do k2 = 1, k1
            
            var_chi_u(k1,k2,b) = var_chi_u(k1,k2,b)/real(bin_pts2d(b))
            
	    covar_chi_u_temp(k1,k2,b) = covar_chi_u_temp(k1,k2,b)/real(bin_pts2d(b))
	    covar_chi_u_rh(k1,k2,b) = covar_chi_u_rh(k1,k2,b)/real(bin_pts2d(b))
         end do 
	 do k2 = k1+1, nk
	    covar_chi_u_temp(k1,k2,b) = covar_chi_u_temp(k1,k2,b)/real(bin_pts2d(b))
	    covar_chi_u_rh(k1,k2,b) = covar_chi_u_rh(k1,k2,b)/real(bin_pts2d(b))
         end do
      end do
   end do

   do b = 1, num_bins2d
      do k1 = 1, nk   
	 do k2 = k1+1, nk
            
            var_chi_u(k1,k2,b) = var_chi_u(k2,k1,b) 
         end do
      end do
   end do
     

   write(6,'(2a)')' [8] Calculate eigenvectors, eigenvalues and inverse for chi_u/chi_u covariance '

   mmax = nk
   variable='chi_u'
   call compute_var_inv(variance_threshold, nk, num_bins2d, mmax, var_chi_u, var_chi_u_inv, testing_eofs, variable)


   write(6,'(2a)')' [9] Calculate regression coefficients with chi_u'


   date = start_date
   cdate = sdate


   do b = 1, num_bins2d
      do k1 = 1, nk
         summ = 0.0
         do k2 = 1, nk
            summ = summ + covar_chi_u_ps(k2,b) * var_chi_u_inv(k2,k1,b)
         end do
         regcoeff_chi_u_ps(k1,b) = summ
      end do
   end do


   do b = 1, num_bins2d
      do k1 = 1, nk
         do k2 = 1, nk
            summ = 0.0
            summ2 = 0.0
            do k3 = 1, nk
               summ = summ + covar_chi_u_temp(k1,k3,b) * var_chi_u_inv(k3,k2,b)
               summ2 = summ2 + covar_chi_u_rh(k1,k3,b) * var_chi_u_inv(k3,k2,b)
            end do
            regcoeff_chi_u_t(k1,k2,b) = summ
            regcoeff_chi_u_rh(k1,k2,b) = summ2
         end do
      end do
   end do
         

   write(6,'(a)')' [10] Compute control variable fields:'


   date = start_date
   cdate = sdate
   count = 0


   do while ( cdate <= edate)


      do member = 1, ne

         count = count + 1
         count_inv = 1.0 / real(count)

         write(ce,'(i3.3)')member


         variable = 'psi'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)psi
         close(iunit)


         variable = 'chi_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)chi
         close(iunit)


         variable = 'ps'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)ps
         close(iunit)


         variable = 't'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)temp
         close(iunit)


         do j = 1, nj
            do i = 1, ni
               b = bin2d(i,j)
               ps(i,j) = ps(i,j) - SUM(regcoeff_psi_ps(1:nk,b) * psi(i,j,1:nk)) &
                    - SUM(regcoeff_chi_u_ps(1:nk,b) * chi(i,j,1:nk)) 
            end do
         end do
         
         
         ps_u_mean = ( real( count-1) *ps_u_mean + ps) * count_inv
         

         do j = 1, nj
            do i = 1, ni
               b = bin2d(i,j)
               do k = 1, nk
                  temp(i,j,k) = temp(i,j,k) - SUM(regcoeff_psi_t(k,1:nk,b) * psi(i,j,1:nk)) &
                       - SUM(regcoeff_chi_u_t(k,1:nk,b) * chi(i,j,1:nk)) 
               end do
            end do
         end do

         
         temp_u_mean = ( real( count-1) *temp_u_mean + temp) * count_inv

         variable = 'ps_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, 1
         write(ounit)ps
         close(ounit)

         variable = 't_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, nk
         write(ounit)temp
         close(ounit)

      end do 


      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate
      
   end do 


   write(6,'(a)')' [10.5] Subtract the mean from the control variable'


   
   
   do while ( cdate <= edate )

      
      do member = 1, ne

         write(ce,'(i3.3)')member

         
         variable = 't_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)temp
         close(iunit)

         temp = temp - temp_u_mean

         
         variable = 'ps_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)ps
         close(iunit)

         ps = ps - ps_u_mean
         
         
         variable = 't_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, nk
         write(ounit)temp
         close(ounit)

         
         variable = 'ps_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, 1
         write(ounit)ps
         close(ounit)

      end do
      
      
      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate

   end do


   write(6,'(2a)')' [11] Calculate correlations of rh with unbalanced temp and ps'


   date = start_date
   cdate = sdate


   do while ( cdate <= edate )
   

      do member = 1, ne

         write(ce,'(i3.3)')member


         variable = 'ps_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)ps
         close(iunit)


         variable = 't_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)temp
         close(iunit)


         variable = 'rh'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)rh     
         close(iunit)


         call compute_covar2d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, ps, rh, covar_ps_u_rh)
         

         call compute_covar3d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, rh, temp, covar_temp_u_rh)


         call compute_var2d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, ps, var_ps_u)


         call compute_var3d(ni, nj, nk, bin2d, num_bins2d, bin_pts2d, temp, var_temp_u)

      end do  


      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate

   end do   

   do b = 1, num_bins2d
      do k = 1, nk    
         covar_ps_u_rh(k,b) = covar_ps_u_rh(k,b)/real(bin_pts2d(b))
      end do
      var_ps_u(b) = var_ps_u(b)/real(bin_pts2d(b))
   end do

   do b = 1, num_bins2d
      do k1 = 1, nk    
	 do k2 = 1, k1
            var_temp_u(k1,k2,b) = var_temp_u(k1,k2,b)/real(bin_pts2d(b))
            covar_temp_u_rh(k1,k2,b) = covar_temp_u_rh(k1,k2,b)/real(bin_pts2d(b))
         end do
	 do k2 = k1+1, nk
	    covar_temp_u_rh(k1,k2,b) = covar_temp_u_rh(k1,k2,b)/real(bin_pts2d(b))
         end do
      end do
   end do

   do b = 1, num_bins2d
      do k1 = 1, nk    
	 do k2 = k1+1, nk
            
            var_temp_u(k1,k2,b) = var_temp_u(k2,k1,b)
         end do
      end do
   end do

   mmax = nk
   variable='temp_u'
   call compute_var_inv(variance_threshold, nk, num_bins2d, mmax, var_temp_u, var_temp_u_inv, testing_eofs, variable)

   date = start_date
   cdate = sdate

                  
         
         do b = 1, num_bins2d
            do k1 = 1, nk 
               regcoeff_ps_u_rh(k1,b) = covar_ps_u_rh(k1,b) / var_ps_u(b)
            end do
         end do

         
         do b = 1, num_bins2d
            do k1 = 1, nk
               do k2 = 1, nk
                  summ = 0.0
                  do k3 = 1, nk
                     summ = summ + covar_temp_u_rh(k1,k3,b) * var_temp_u_inv(k3,k2,b)
                  end do
                  regcoeff_t_u_rh(k1,k2,b) = summ
               end do
            end do
         end do
         

   write(6,'(a)')' [13] Compute control variable fields: rh_u'


   date = start_date
   cdate = sdate


   do while ( cdate <= edate)


      do member = 1, ne

         write(ce,'(i3.3)')member


         variable = 'psi'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)psi
         close(iunit)


         variable = 'chi_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)chi
         close(iunit)


         variable = 'ps_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce//'.01'
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nkdum
         read(iunit)ps
         close(iunit)


         variable = 't_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)temp
         close(iunit)


         variable = 'rh'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (iunit, file = filename, form='unformatted')
         read(iunit)ni, nj, nk
         read(iunit)rh     
         close(iunit)


         do j = 1, nj
            do i = 1, ni
               b = bin2d(i,j)
               do k = 1, nk
                  rh(i,j,k) = rh(i,j,k) &
                       - SUM(regcoeff_psi_rh(k,1:nk,b) * psi(i,j,1:nk)) &
                       - SUM(regcoeff_chi_u_rh(k,1:nk,b) * chi(i,j,1:nk)) &
                       - SUM(regcoeff_t_u_rh(k,1:nk,b) * temp(i,j,1:nk))
               end do
            end do
         end do

         do j = 1, nj
            do i = 1, ni
               b = bin2d(i,j)
               do k = 1, nk
                  rh(i,j,k) = rh(i,j,k) - regcoeff_ps_u_rh(k,b) * ps(i,j) 
               end do
            end do
         end do


         variable = 'rh_u'
         filename = trim(variable)//'/'//date(1:10)
         filename = trim(filename)//'.'//trim(variable)//'.e'//ce
         open (ounit, file = filename, form='unformatted')
         write(ounit)ni, nj, nk
         write(ounit)rh
         close(ounit)
      
      end do


      call da_advance_cymdh( date, interval, new_date )
      date = new_date
      read(date(1:10), fmt='(i10)')cdate
   end do



   filename = 'gen_be_stage2.dat'
   open (ounit, file = filename, form='unformatted')
   write(ounit)ni, nj, nk
   write(ounit)num_bins, num_bins2d

   variable = 'regcoeff_psi_chi'
   write(ounit)variable
   write(ounit)regcoeff_psi_chi

   variable = 'regcoeff_psi_t'
   write(ounit)variable
   write(ounit)regcoeff_psi_t

   variable = 'regcoeff_psi_ps'
   write(ounit)variable
   write(ounit)regcoeff_psi_ps

   variable = 'regcoeff_psi_rh'
   write(ounit)variable
   write(ounit)regcoeff_psi_rh

   variable = 'regcoeff_chi_u_t'
   write(ounit)variable
   write(ounit)regcoeff_chi_u_t

   variable = 'regcoeff_chi_u_ps'
   write(ounit)variable
   write(ounit)regcoeff_chi_u_ps

   variable = 'regcoeff_chi_u_rh'
   write(ounit)variable
   write(ounit)regcoeff_chi_u_rh

   variable = 'regcoeff_t_u_rh'
   write(ounit)variable
   write(ounit)regcoeff_t_u_rh

   variable = 'regcoeff_ps_u_rh'
   write(ounit)variable
   write(ounit)regcoeff_ps_u_rh

   close(ounit)


   deallocate(regcoeff_psi_chi)
   deallocate(regcoeff_psi_ps)
   deallocate(regcoeff_psi_t)
   deallocate(regcoeff_psi_rh)

   deallocate( regcoeff_chi_u_ps)
   deallocate( regcoeff_chi_u_t)
   deallocate( regcoeff_chi_u_rh)

   deallocate( regcoeff_ps_u_rh)
   deallocate( regcoeff_t_u_rh)


contains

subroutine compute_covar3d (ni, nj, nk, bin2d, num_bins2d, bin_pts2d, field3d_1, field3d_2, covar)

implicit none

integer, intent(in)  :: ni, nj, nk             
integer, intent(in)  :: num_bins2d
integer, intent(in)  :: bin2d(1:ni,1:nj)
integer, intent(in)  :: bin_pts2d(1:num_bins2d)
real, intent(in)     :: field3d_1(1:ni,1:nj,1:nk), field3d_2(1:ni,1:nj,1:nk)
real, intent(inout)  :: covar(1:nk,1:nk,1:num_bins2d)

integer              :: i, j, k1, k2
integer              :: b                  



  do i = 1, ni
    do j = 1, nj
       b = bin2d(i,j)
       do k1 = 1, nk
       	  do k2 = 1, nk
             covar(k1,k2,b) = covar(k1,k2,b) + field3d_1(i,j,k1)*field3d_2(i,j,k2)
          end do
       end do
    end do
  end do

end subroutine compute_covar3d

subroutine compute_covar_bin3d (ni, nj, nk, bin, num_bins, bin_pts, field3d_1, field3d_2, covar)



implicit none

integer, intent(in)    :: ni, nj, nk             
integer, intent(in)    :: num_bins
integer, intent(in)    :: bin(1:ni,1:nj,1:nk)
integer, intent(in)    :: bin_pts(1:num_bins)
real, intent(in)       :: field3d_1(1:ni,1:nj,1:nk), field3d_2(1:ni,1:nj,1:nk) 
real, intent(inout)    :: covar(1:num_bins) 

integer                :: i, j, k
integer                :: b



  do k = 1, nk
    do j = 1, nj
      do i = 1, ni
         b = bin(i,j,k)
	 covar(b) = covar(b) + field3d_1(i,j,k) * field3d_2(i,j,k)
      end do
    end do
  end do

end subroutine compute_covar_bin3d

subroutine compute_var_bin3d(ni, nj, nk, bin, num_bins, bin_pts, field3d, var)

implicit none

integer, intent(in)    :: ni, nj, nk             
integer, intent(in)    :: num_bins
integer, intent(in)    :: bin(1:ni,1:nj,1:nk)
integer, intent(in)    :: bin_pts(1:num_bins)
real, intent(in)       :: field3d(1:ni,1:nj,1:nk)
real, intent(inout)    :: var(1:num_bins) 

integer                :: i, j, k
integer                :: b



  do k = 1, nk
    do j = 1, nj
      do i = 1, ni
         b = bin(i,j,k)
	 var(b) = var(b) + field3d(i,j,k) * field3d(i,j,k)
      end do
    end do
  end do


end subroutine compute_var_bin3d

subroutine compute_covar2d (ni, nj, nk, bin2d, num_bins2d, bin_pts2d, field2d, field3d, covar)

implicit none

integer, intent(in)    :: ni, nj, nk             
integer, intent(in)    :: num_bins2d
integer, intent(in)    :: bin2d(1:ni,1:nj)
integer, intent(in)    :: bin_pts2d(1:num_bins2d)
real, intent(in)       :: field2d(1:ni,1:nj), field3d(1:ni,1:nj,1:nk)
real, intent(inout)      :: covar(1:nk, 1:num_bins2d) 

integer                :: i, j, k
integer                :: b



  do i = 1, ni
    do j = 1, nj
      b = bin2d(i,j)
      do k = 1, nk
         covar(k,b) = covar(k,b) + field2d(i,j) * field3d(i,j,k)
      end do
    end do
  end do

end subroutine compute_covar2d

subroutine compute_var2d (ni, nj, nk, bin2d, num_bins2d, bin_pts2d, field2d, var)

implicit none

integer, intent(in)    :: ni, nj, nk             
integer, intent(in)    :: num_bins2d
integer, intent(in)    :: bin2d(1:ni,1:nj)
integer, intent(in)    :: bin_pts2d(1:num_bins2d)
real, intent(in)       :: field2d(1:ni,1:nj)
real, intent(inout)    :: var(1:num_bins2d) 

integer                :: i, j, k
integer                :: b



  do i = 1, ni
    do j = 1, nj
     	b = bin2d(i,j)
	var(b) = var(b) + field2d(i,j) * field2d(i,j)
     end do
  end do


end subroutine compute_var2d

subroutine compute_var3d (ni, nj, nk, bin2d, num_bins2d, bin_pts2d, field3d, var)

implicit none

integer, intent(in)  :: ni, nj, nk             
integer, intent(in)  :: num_bins2d
integer, intent(in)  :: bin2d(1:ni,1:nj)
integer, intent(in)  :: bin_pts2d(1:num_bins2d)
real, intent(in)     :: field3d(1:ni,1:nj,1:nk)
real, intent(inout)  :: var(1:nk,1:nk,1:num_bins2d)

integer              :: i, j ,k1, k2
integer              :: b                  


  do i = 1, ni
    do j = 1, nj
       b = bin2d(i,j)
       do k1 = 1, nk
          do k2 = 1, k1
             var(k1,k2,b) = var(k1,k2,b) + field3d(i,j,k1)*field3d(i,j,k2)
           end do
       end do
    end do
  end do



end subroutine compute_var3d

subroutine compute_var_inv (variance_threshold, nk, num_bins2d, mmax, var, var_inv, testing_eofs, variable)


   use da_gen_be, only : da_eof_decomposition,da_eof_decomposition_test
   use da_control, only : stdout, stderr, filename_len	

   implicit none

   integer, intent(in) :: nk                         
   integer, intent(in) :: num_bins2d                 

   real, intent(in)    :: variance_threshold         
   integer, intent(inout):: mmax                       
   real, intent(in)    :: var(1:nk,1:nk,1:num_bins2d) 
   real, intent(inout)   :: var_inv(1:nk,1:nk,1:num_bins2d) 
   
   logical, intent(in) :: testing_eofs         
   character*10, intent(in) :: variable                   
   
   character(len=filename_len) :: filename           


   integer             :: k, k1, k2, m               
   integer             :: b                          
   integer 	       :: ounit	 
   real                :: summ                       
   real                :: total_variance             
   real                :: cumul_variance             

   real*8, allocatable   :: eval(:)                    
   real*8, allocatable   :: evec(:,:)                  
   real, allocatable   :: work(:,:)                  
   real, allocatable   :: LamInvET(:,:)              

   allocate( eval(1:nk) )
   allocate( evec(1:nk,1:nk) )
   allocate( work(1:nk,1:nk) )
   allocate( LamInvET(1:nk,1:nk) )


   do b = 1, num_bins2d  
 
      LamInvET(:,:) = 0.0
      work(1:nk,1:nk) = var(1:nk,1:nk,b)
      call da_eof_decomposition( nk, work, evec, eval )

      if ( testing_eofs ) then
         call da_eof_decomposition_test( nk, work, evec, eval )
      end if
   

      summ = 0.0
      do k = 1, nk
         summ = summ + eval(k)
      end do
      total_variance = summ
   

      cumul_variance = 0.0
      mmax = nk
      do k = 1, nk
         cumul_variance = cumul_variance + eval(k) / total_variance
         if ( cumul_variance > 1.0 - variance_threshold ) then
            mmax = k - 1
            exit
         end if
      end do
      write(6,'(2(a,i6),2(a,1pe11.5))') ' Bin = ', b, ', <psipsi> truncation = ', mmax, &
                                        ', Total Variance = ', total_variance, &
                                        ', Condition number = ', eval(1) / eval(nk-1)
   

      do k = 1, nk
         do m = 1, mmax
            LamInvET(m,k) = evec(k,m) / eval(m)
         end do
      end do


   
      do k1 = 1, nk
         do k2 = 1, k1
            summ = 0.0
            do m = 1, nk
               summ = summ + evec(k1,m) * LamInvET(m,k2)
            end do
            var_inv(k1,k2,b) = summ
         end do
      end do
   
      do k1 = 1, nk
         do k2 = k1+1, nk 
            var_inv(k1,k2,b) = var_inv(k2,k1,b)
         end do
      end do
   
   end do

   deallocate(eval)
   deallocate(evec)
   deallocate(work)
   deallocate(LamInvET)

end subroutine compute_var_inv
end program gen_mbe_stage2
