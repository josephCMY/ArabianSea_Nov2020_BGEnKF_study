












program gen_be_diags_read











   use da_control, only : do_normalize,filename_len,stderr,stdout,use_rf
   use da_tools_serial, only : da_get_unit
   use da_gen_be, only : da_print_be_stats_h_regional, &
      da_print_be_stats_h_global, da_print_be_stats_p, da_print_be_stats_v

   implicit none

   character*10        :: variable          
   character*80        :: var80c

   character*8         :: uh_method         
   integer             :: n_smth_sl                  
   integer             :: cv_options                 
   character(len=filename_len)        :: filename          
   integer             :: outunit           
   integer             :: ni, nj, nk, nk_3d 
   integer             :: bin_type          
   integer             :: num_bins          
   integer             :: num_bins2d        
   integer             :: k                 
   integer             :: kdum              
   integer             :: max_wavenumber    
   real                :: lat_min, lat_max  
   real                :: binwidth_lat      
   real                :: binwidth_hgt      
   real                :: hgt_min, hgt_max  
   real                :: scale_length_ps   
   logical             :: dummy

   integer, allocatable:: bin(:,:,:)        
   integer, allocatable:: bin2d(:,:)        
   real, allocatable   :: regcoeff_psi_chi(:)        
   real, allocatable   :: regcoeff_psi_t(:,:,:)   
   real, allocatable   :: regcoeff_psi_ps(:,:)    
   real, allocatable   :: regcoeff_psi_rh(:,:,:)  
   real, allocatable   :: regcoeff_chi_u_t(:,:,:) 
   real, allocatable   :: regcoeff_chi_u_ps(:,:)  
   real, allocatable   :: regcoeff_chi_u_rh(:,:,:)
   real, allocatable   :: regcoeff_t_u_rh(:,:,:)  
   real, allocatable   :: regcoeff_ps_u_rh(:,:)   

   real, allocatable   :: e_vec(:,:)        
   real, allocatable   :: e_val(:)          
   real, allocatable   :: e_vec_loc(:,:,:)  
   real, allocatable   :: e_val_loc(:,:)    
   real, allocatable   :: total_power(:)    
   real, allocatable   :: scale_length(:)   


   namelist / gen_be_diags_nl / cv_options, do_normalize, n_smth_sl, uh_method, use_rf


   integer :: iunit,namelist_unit

   
   
   integer, parameter :: ounit = 71

   stderr = 0
   stdout = 6

   uh_method = 'scale'
   n_smth_sl = 0
   cv_options = 5

   call da_get_unit(iunit)
   call da_get_unit(namelist_unit)

   open(unit=namelist_unit, file='gen_be_diags_nl.nl', &
        form='formatted', status='old', action='read')
   read(namelist_unit, gen_be_diags_nl)
   close(namelist_unit)

   filename = 'be.dat'
   print '("*** Unit=",i3,3X,"filename=",a40)',iunit, filename
   open (iunit, file = filename, form='unformatted')

   
   
   

   
   read(iunit)ni, nj, nk
   nk_3d = nk

   allocate( bin(1:ni,1:nj,1:nk) )
   allocate( bin2d(1:ni,1:nj) )

   

   read(iunit)bin_type
   read(iunit)lat_min, lat_max, binwidth_lat
   read(iunit)hgt_min, hgt_max, binwidth_hgt
   read(iunit)num_bins, num_bins2d
   read(iunit)bin(1:ni,1:nj,1:nk)
   read(iunit)bin2d(1:ni,1:nj)


   if ( cv_options /= 7 ) then
      
      allocate  (regcoeff_psi_chi(1:num_bins))
      allocate  (regcoeff_psi_t(1:nk,1:nk,1:num_bins2d))
      allocate  (regcoeff_psi_ps(1:nk,1:num_bins2d))
      allocate  (regcoeff_psi_rh(1:nk,1:nk,1:num_bins2d))
      allocate  (regcoeff_chi_u_t(1:nk,1:nk,1:num_bins2d))
      allocate  (regcoeff_chi_u_ps(1:nk,1:num_bins2d))
      allocate  (regcoeff_chi_u_rh(1:nk,1:nk,1:num_bins2d))
      allocate  (regcoeff_t_u_rh(1:nk,1:nk,1:num_bins2d))
      allocate  (regcoeff_ps_u_rh(1:nk,1:num_bins2d))

      regcoeff_psi_chi = 0.
      regcoeff_psi_t   = 0.
      regcoeff_psi_ps  = 0.
      regcoeff_psi_rh  = 0.
      regcoeff_chi_u_t = 0.
      regcoeff_chi_u_ps= 0.
      regcoeff_chi_u_rh= 0.
      regcoeff_t_u_rh  = 0.
      regcoeff_ps_u_rh = 0.

      if ( cv_options == 5 ) then
         read (iunit) regcoeff_psi_chi
         read (iunit) regcoeff_psi_ps
         read (iunit) regcoeff_psi_t
      else
         do k = 1 , 9
            read (iunit) var80c
            select case( trim(adjustl(var80c)) )
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
            case default;
               write(6,fmt='(A,A)')'Error in reading regression coefficients for: ',trim(adjustl(var80c))
               stop
            end select
         end do
      end if

      outunit = ounit + 100
      call da_print_be_stats_p( outunit, ni, nj, nk, num_bins, num_bins2d, &
                    bin, bin2d, regcoeff_psi_chi,regcoeff_psi_ps,regcoeff_psi_t)     
   else
      
      
      outunit = ounit + 103 
   end if

   
   
   

   read(iunit)variable
   read(iunit)nk, num_bins2d

   allocate( e_vec(1:nk,1:nk) )
   allocate( e_val(1:nk) )
   allocate( e_vec_loc(1:nk,1:nk,1:num_bins2d) )
   allocate( e_val_loc(1:nk,1:num_bins2d) )

   read(iunit)e_vec
   read(iunit)e_val
   read(iunit)e_vec_loc
   read(iunit)e_val_loc
   call da_print_be_stats_v( outunit, variable, nk, num_bins2d, &
                             e_vec, e_val, e_vec_loc, e_val_loc )
   read(iunit)variable
   read(iunit)nk, num_bins2d
   read(iunit)e_vec
   read(iunit)e_val
   read(iunit)e_vec_loc
   read(iunit)e_val_loc
   call da_print_be_stats_v( outunit, variable, nk, num_bins2d, &
                             e_vec, e_val, e_vec_loc, e_val_loc )

   read(iunit)variable
   read(iunit)nk, num_bins2d
   read(iunit)e_vec
   read(iunit)e_val
   read(iunit)e_vec_loc
   read(iunit)e_val_loc
   call da_print_be_stats_v( outunit, variable, nk, num_bins2d, &
                             e_vec, e_val, e_vec_loc, e_val_loc )

   read(iunit)variable
   read(iunit)nk, num_bins2d
   read(iunit)e_vec
   read(iunit)e_val
   read(iunit)e_vec_loc
   read(iunit)e_val_loc
   call da_print_be_stats_v( outunit, variable, nk, num_bins2d, &
                             e_vec, e_val, e_vec_loc, e_val_loc )

   deallocate( e_vec )
   deallocate( e_val )
   deallocate( e_vec_loc )
   deallocate( e_val_loc )

   if (uh_method /= 'power') then
      

      read(iunit)variable
      read(iunit)nk, num_bins2d

      allocate( e_vec(1:nk,1:nk) )
      allocate( e_val(1:nk) )
      allocate( e_vec_loc(1:nk,1:nk,1:num_bins2d) )
      allocate( e_val_loc(1:nk,1:num_bins2d) )

      read(iunit)e_vec
      read(iunit)e_val
      read(iunit)e_vec_loc
      read(iunit)e_val_loc
      call da_print_be_stats_v( outunit, variable, nk, num_bins2d, &
                                e_vec, e_val, e_vec_loc, e_val_loc )
      deallocate( e_vec )
      deallocate( e_val )
      deallocate( e_vec_loc )
      deallocate( e_val_loc )
   end if
   
   nk = nk_3d

   if (uh_method == 'power') then
       write(6,'(/a)') '[3] Gather horizontal error power spectra.'

      do k = 1, nk
         read(iunit)variable
         read(iunit)max_wavenumber, kdum
         if ( k == 1 ) allocate( total_power(0:max_wavenumber) )
         read(iunit) dummy 
         read(iunit)total_power(:)
         call da_print_be_stats_h_global( outunit, variable, k, max_wavenumber, total_power )
      end do

      do k = 1, nk
         read(iunit)variable
         read(iunit)max_wavenumber, kdum
         read(iunit) dummy 
         read(iunit)total_power(:)
         call da_print_be_stats_h_global( outunit, variable, k, max_wavenumber, total_power )
      end do

      do k = 1, nk
         read(iunit)variable
         read(iunit)max_wavenumber, kdum
         read(iunit) dummy 
         read(iunit)total_power(:)
         call da_print_be_stats_h_global( outunit, variable, k, max_wavenumber, total_power )
      end do

      do k = 1, nk
         read(iunit)variable
         read(iunit)max_wavenumber, kdum
         read(iunit) dummy 
         read(iunit)total_power(:)
         call da_print_be_stats_h_global( outunit, variable, k, max_wavenumber, total_power )
      end do

      read(iunit)variable
      read(iunit)max_wavenumber, kdum
      read(iunit) dummy 
      read(iunit)total_power(:)
      call da_print_be_stats_h_global( outunit, variable, k, max_wavenumber, total_power )

   else if (uh_method == 'scale   ') then

      allocate (scale_length(1:nk))

      
      read(iunit) variable
      read(iunit) scale_length
      call da_print_be_stats_h_regional( outunit, variable, nk, scale_length )

      
      read(iunit) variable
      read(iunit) scale_length
      call da_print_be_stats_h_regional( outunit, variable, nk, scale_length )

      
      read(iunit) variable
      read(iunit) scale_length
      call da_print_be_stats_h_regional( outunit, variable, nk, scale_length )

      
      read(iunit) variable
      read(iunit) scale_length
      call da_print_be_stats_h_regional( outunit, variable, nk, scale_length )

      
      read(iunit) variable
      read(iunit) scale_length_ps 
      write(6,'(3a,i5)')' Scale length for variable ', trim(variable), ' in unit ', outunit
      write(outunit,'(a,i4,1pe15.5)')trim(variable), 1, scale_length_ps
      outunit = outunit + 1
      write(6,*)

      deallocate (scale_length)

   else
      write(6,'(a,": uh_method=",a," not implemented.")')"gen_be_diags_read.b",uh_method
   endif

   close(iunit)

end program gen_be_diags_read
