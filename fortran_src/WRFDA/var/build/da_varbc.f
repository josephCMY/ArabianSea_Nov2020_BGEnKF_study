












module da_varbc

   
   
   

   use module_dm, only : wrf_dm_sum_real, wrf_dm_sum_reals, wrf_dm_sum_integer
   use module_radiance, only : q2ppmv, satinfo
   use da_control, only : trace_use,missing_r, qc_varbc_bad, rtm_option, &
      stdout,rtm_option_rttov,rtm_option_crtm, filename_len, cv_size_domain, &
      cv_size_domain_jp, use_varbc, freeze_varbc, varbc_factor, varbc_nobsmin, &
      rootproc, varbc_nbgerr, ierr, comm, max_ext_its
   use da_define_structures, only : iv_type, y_type, be_type, &
      varbc_info_type,varbc_type
   use da_radiance1, only : stats_rad_type
   use da_radiance1, only : da_predictor_crtm
   use da_reporting, only : da_error, message, da_warning, da_message
   use da_tools, only : da_eof_decomposition
   use da_tools_serial, only : da_free_unit, da_get_unit
   use da_tracing, only : da_trace_entry, da_trace_exit, da_trace, &
      da_trace_int_sort
   
   implicit none
   
contains

  subroutine da_varbc_direct (iv)

   !---------------------------------------------------------------------------
   !  PURPOSE: Apply bias correction to radiance innovations.
   !
   !  METHOD:  d   = y - H(x) - bc
   !          bc   = SUM( beta_i Pred_i )
   !                  i
   !          beta = VarBC parameters
   !          Pred = Bias predictors
   !
   !  Called from da_get_innov_vector_radiance
   !
   !  HISTORY: 10/26/2007 - Creation                     Tom Auligne
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout)  :: iv        ! O-B structure.

   integer              :: inst, k, n, npred, npredmax
    
   if (trace_use) call da_trace_entry("da_varbc_direct")

      do inst = 1, iv%num_inst                 ! loop for sensor
         npredmax = iv%instid(inst)%varbc_info%npredmax
         write(unit=stdout,fmt='(A,A)') 'VARBC: Applying bias correction for ', &
            trim(iv%instid(inst)%rttovid_string)
	          
         do k=1,iv%instid(inst)%nchan          ! loop for channel
            npred          = iv%instid(inst)%varbc(k)%npred
	    if (npred <= 0) cycle	               ! VarBC channel only

            do n= iv%instid(inst)%info%n1,iv%instid(inst)%info%n2  ! loop for pixel      	     
              ! apply bias correction through linear regression
              !-------------------------------------------------
              if ( iv%instid(inst)%tb_inv(k,n) > missing_r ) then
                 iv%instid(inst)%tb_inv(k,n) = iv%instid(inst)%tb_inv(k,n) - &
                       SUM( iv%instid(inst)%varbc(k)%param(1:npred) * &
                            iv%instid(inst)%varbc_info%pred( &
                            iv%instid(inst)%varbc(k)%ipred(1:npred),n) )
              end if
	    end do
         end do
      end do                                   ! end loop for sensor

   if (trace_use) call da_trace_exit("da_varbc_direct")

 end subroutine da_varbc_direct
  subroutine da_varbc_tl (cv_size, cv, iv, y)

   !---------------------------------------------------------------------------
   !  PURPOSE: Apply bias correction to radiance innovations.
   !
   !  METHOD: delta_d    = y - delta_bc
   !          y = H (delta_x)
   !          delta_bc   = SUM( delta_beta_i Pred_i )
   !                       i
   !          delta_beta = VarBC parameters
   !          Pred = Bias predictors
   ! 
   !  Called from da_transform_vtoy
   !
   !  HISTORY: 10/27/2007 - Creation                     Tom Auligne
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)           :: cv_size         ! Size of cv array.
   real, intent(in)              :: cv(1:cv_size)   ! control variables.
   type (iv_type), intent(in)    :: iv              ! O-B structure.
   type (y_type), intent(inout)  :: y               ! y = h (xa)

   integer              :: inst, i, k, num_rad, n, npred
   real, allocatable    :: varbc_param_tl(:)
   
   if (trace_use) call da_trace_entry("da_varbc_tl")

      do inst = 1, iv%num_inst                 ! loop for sensor
	 num_rad = iv%instid(inst)%num_rad
         if (num_rad < 1) cycle

         allocate(varbc_param_tl(iv%instid(inst)%varbc_info%npredmax))
  
         do k=1,iv%instid(inst)%nchan          ! loop for channel
            npred = iv%instid(inst)%varbc(k)%npred
	    if (npred <= 0) cycle              ! VarBC channels only
	    	     
           !---------------------------------------------------------------
           ! Change of variable (preconditioning) 
           !---------------------------------------------------------------
            varbc_param_tl(:) = 0.0
	    do i = 1, npred
    	       varbc_param_tl(i) = SUM( cv(iv%instid(inst)%varbc(k)%index(1:npred)) * &
	                                iv%instid(inst)%varbc(k)%vtox(i,1:npred) )
	    end do	
		
           !---------------------------------------------------------------
           ! TL of bias correction through linear regression
           !---------------------------------------------------------------
	    do n = 1, num_rad                                         ! loop for pixel      
               if (iv%instid(inst)%tb_qc(k,n) <= qc_varbc_bad) cycle  ! good obs only

  	        y%instid(inst)%tb(k,n) = y%instid(inst)%tb(k,n) +   &
		   SUM( varbc_param_tl(1:npred) *                   &
		        iv%instid(inst)%varbc_info%pred(            &
		        iv%instid(inst)%varbc(k)%ipred(1:npred),n) )
            end do
         end do

         deallocate(varbc_param_tl)
      end do                                   ! end loop for sensor

   if (trace_use) call da_trace_exit("da_varbc_tl")

 end subroutine da_varbc_tl
  subroutine da_varbc_adj (cv_size, cv, iv, y)

   !---------------------------------------------------------------------------
   !  PURPOSE: Apply bias correction to radiance innovations.
   !
   !  Called from da_transform_vtoy_adj
   !
   !  HISTORY: 10/27/2007 - Creation                     Tom Auligne
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)           :: cv_size         ! Size of cv array.
   real, intent(inout)           :: cv(1:cv_size)   ! control variables.
   type (iv_type), intent(inout) :: iv              ! O-B structure.
   type (y_type), intent(in)     :: y               ! y = h (xa)

   integer              :: inst, i, k, num_rad, n, npred
   real, allocatable    :: varbc_param_adj(:)
   real                 :: cv_local

   if (trace_use) call da_trace_entry("da_varbc_adj")

      do inst = 1, iv%num_inst                 ! loop for sensor
         num_rad = iv%instid(inst)%num_rad
   
         allocate(varbc_param_adj(iv%instid(inst)%varbc_info%npredmax))

         do k=1,iv%instid(inst)%nchan          ! loop for channel
            npred = iv%instid(inst)%varbc(k)%npred
	    if (npred <= 0) cycle              ! VarBC channels only
	   
           !---------------------------------------------------------------
           ! Adjoint to bias correction through linear regression
           !---------------------------------------------------------------
            varbc_param_adj(:) = 0.0
	    if (num_rad >0) then
	       do n = 1, num_rad                                     ! loop for pixel      
                  if (iv%instid(inst)%tb_qc(k,n) <= qc_varbc_bad) cycle  ! good obs only
                  if (.not. iv%instid(inst)%info%proc_domain(1,n))  cycle   ! do not sum up HALO data
	       
		  varbc_param_adj(1:npred) = varbc_param_adj(1:npred)      + &
	       		   		     y%instid(inst)%tb(k,n)        * &
		                            iv%instid(inst)%varbc_info%pred( &
					    iv%instid(inst)%varbc(k)%ipred(1:npred),n ) 
	       end do	        
	    end if
		       
           !---------------------------------------------------------------
           ! Change of variable (preconditioning) + sum across processors
           !---------------------------------------------------------------
	    do i = 1, npred
               cv_local = SUM(varbc_param_adj(1:npred) * iv%instid(inst)%varbc(k)%vtox(i,1:npred))
               cv(iv%instid(inst)%varbc(k)%index(i)) = cv(iv%instid(inst)%varbc(k)%index(i)) + &
	                                               wrf_dm_sum_real(cv_local)	
            end do
	 end do
	 
         deallocate(varbc_param_adj)	     
      end do                                   ! end loop for sensor

   if (trace_use) call da_trace_exit("da_varbc_adj")

 end subroutine da_varbc_adj
subroutine da_varbc_pred(iv)

   !---------------------------------------------------------------------------
   !  PURPOSE: Calculate bias predictors
   !
   ! pred(1) - 1 (Constant)
   ! pred(2) - 1000-300 thickness 43 1005.43-521.46 thickness
   ! pred(3) - 200-50 thickness   43  194.36-56.73  thickness
   ! pred(4) - T_skin
   ! pred(5) - total column precipitable water
   ! pred(6) - satellite scan position
   ! pred(7) - satellite scan position **2
   ! pred(8) - satellite scan position **3
   ! pred(9) - 1 Gamma Correction
   !
   !  Called from da_varbc_coldstart
   !
   !  HISTORY: 10/26/2007 - Creation                  Tom Auligne
   !---------------------------------------------------------------------------
   
  implicit none

  integer,parameter   :: npred_hk = 4       ! Number of H&K bias predictors

  type (iv_type), intent(inout) :: iv       ! O-B structure.

  integer   :: inst, nlevels, i, npredmax, n, num_rad
  real      :: pred_hk(npred_hk)

  if ( iv%num_inst < 1 ) RETURN

  if (trace_use) call da_trace_entry("da_varbc_pred")

  do inst = 1, iv%num_inst                                  ! loop for sensor
    npredmax = iv%instid(inst)%varbc_info%npredmax
    if (npredmax <= 0) cycle                                ! VarBC instr only
       
    num_rad = iv%instid(inst)%num_rad
    if (iv%instid(inst)%info%n2 < iv%instid(inst)%info%n1) cycle
    do n = iv%instid(inst)%info%n1, iv%instid(inst)%info%n2  ! loop for pixel
      ! get H&K predictors
      !--------------------
      if (rtm_option==rtm_option_rttov) then
      else if (rtm_option==rtm_option_crtm) then
        nlevels=iv%instid(inst)%nlevels-1
        call da_predictor_crtm( &
	     pred_hk(1:npred_hk), npred_hk, nlevels,&
             iv%instid(inst)%tm(1:nlevels,n), &
             iv%instid(inst)%qm(1:nlevels,n), &
	     iv%instid(inst)%ts(n), &
	     iv%instid(inst)%pf(0:nlevels,n))
      end if     

      ! Populate predictors
      !---------------------  
      iv%instid(inst)%varbc_info%pred(1:npredmax,n) = 0.0
  
      ! Constant predictor
      iv%instid(inst)%varbc_info%pred(1,n) = 1.0 
  
      ! H&K predictors
      if (npredmax >= 2) iv%instid(inst)%varbc_info%pred(2,n) = pred_hk(1)
      if (npredmax >= 3) iv%instid(inst)%varbc_info%pred(3,n) = pred_hk(2)
      if (npredmax >= 4) iv%instid(inst)%varbc_info%pred(4,n) = pred_hk(3)
      if (npredmax >= 5) iv%instid(inst)%varbc_info%pred(5,n) = pred_hk(4)
		
      ! Scan predictors	
      if (npredmax >= 6) iv%instid(inst)%varbc_info%pred(6,n) = iv%instid(inst)%scanpos(n)
      if (npredmax >= 7) iv%instid(inst)%varbc_info%pred(7,n) = iv%instid(inst)%scanpos(n)**2
      if (npredmax >= 8) iv%instid(inst)%varbc_info%pred(8,n) = iv%instid(inst)%scanpos(n)**3

    end do                       ! pixel loop   
  end do                         ! sensor loop
   
   if (trace_use) call da_trace_exit("da_varbc_pred")

end subroutine da_varbc_pred
  subroutine da_varbc_coldstart (iv)

   !---------------------------------------------------------------------------
   !  PURPOSE: [1]: If a cold-start is needed, calculate mode of histogram 
   !                of (uncorrected) innovations for each channel
   !
   !           [2]: Calculate statistics (mean/std) for cold-started predictors
   !
   !           [3]: Normalize predictors
   !
   !  Called from da_get_innov_vector_radiance
   !
   !  HISTORY: 10/26/2007 - Creation                     Tom Auligne
   !           11/03/2008 - Bug-fix: exclude HALO data for
   !                        computation of predictor statistics    Tom/Hui-Chuan  
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type), intent (inout)   :: iv             ! Innovation

   integer                          :: inst, k, n, i
   integer                          :: npredmax, npred, num_rad, num_rad_domain, num_rad_tot

   integer, parameter               :: nbins   = 200                       ! Number of Hist bins.
   real,    parameter               :: maxhist = 10.                       ! Maximum bin value.
   real,    parameter               :: zbin    = 2 * maxhist / real(nbins) ! Hist bin width.
   integer                          :: ibin                                ! Hist bin number.
   real                             :: modetmp(1), mode
   integer, allocatable             :: hist(:), hist_tot(:)
   
   real, allocatable                :: mean(:), rms(:), mean_tot(:), rms_tot(:)
   integer, allocatable             :: ipred(:)
   logical                          :: global_cs
   
   if ( iv%num_inst < 1 ) RETURN

   if (trace_use) call da_trace_entry("da_varbc_coldstart")   

   do inst = 1, iv%num_inst                                            !! loop for sensors

      npredmax       = iv%instid(inst)%varbc_info%npredmax
      num_rad        = iv%instid(inst)%num_rad
      if (num_rad > 0) then
         num_rad_domain = COUNT(iv%instid(inst)%info%proc_domain(1,1:num_rad)) !do not count HALO
      else
         num_rad_domain = 0
      end if
      num_rad_tot    = wrf_dm_sum_integer(num_rad_domain)
      
      if (npredmax <= 0) cycle                                         !! VarBC instr only
      
      allocate( ipred(npredmax) )
    
    !---------------------------------------------------------------------------
    ! [1]: Calculate mode of histogram of (uncorrected) innovations 
    !---------------------------------------------------------------------------
      do k = 1, iv%instid(inst)%nchan                                  !! loop for channels
         npred          = iv%instid(inst)%varbc(k)%npred
	 ipred(1:npred) = iv%instid(inst)%varbc(k)%ipred(1:npred)
	 if (npred <= 0) cycle                                         !! VarBC channels only
	 if (ALL(iv%instid(inst)%varbc(k)%pred_use /= 0)) cycle        !! Coldstart channels only
	 
	 where (iv%instid(inst)%varbc(k)%pred_use(ipred(1:npred)) == 0) &
	        iv%instid(inst)%varbc(k)%param(1:npred) = 0.0

         if (iv%instid(inst)%varbc(k)%pred_use(1) == 0) then
	    Allocate ( hist(nbins),  hist_tot(nbins))
            hist(:) = 0  
            mode    = 0.0
	 
           ! Accumulate statistics for histogram
           ! -----------------------------------
            do n = 1, num_rad      !! loop for pixel      
               if (iv%instid(inst)%info%proc_domain(1,n)) then ! do not count HALO
                  ibin = NINT( (iv%instid(inst)%tb_inv(k,n)+maxhist)/zbin )
 	          if ((ibin>0).AND.(ibin<=nbins)) &
	             hist(ibin) = hist(ibin) + 1
               end if          
            end do             ! end loop for pixels
	       
           ! Do inter-processor communication to gather statistics
           ! ------------------------------------------------------
	    do ibin = 1, nbins
	       hist_tot(ibin) = wrf_dm_sum_integer(hist(ibin))
	    end do
		   
           ! Determine mode of Histogram
           !----------------------------
            if ( SUM(hist_tot(:)) > 0 ) then
	      modetmp(1:1) = MAXLOC(hist_tot(:))*zbin - maxhist
              mode = modetmp(1)
	    end if
	 
	   ! Use mode to initialize VarBC 
	   !-----------------------------
	    if (iv%instid(inst)%varbc(k)%ipred(1) == 1) &
	        iv%instid(inst)%varbc(k)%param(1) = mode

	    Deallocate ( hist, hist_tot )
            if ( satinfo(inst)%iuse(k) == 1 ) &
               write(unit=stdout,fmt='(A,A,I5,A,F5.2)') 'VARBC: Cold-starting ', &
                  trim(adjustl(iv%instid(inst)%rttovid_string)),iv%instid(inst)%ichan(k),&
	          ' --> ',mode			       	
	 end if
      end do                                                              !  end loop for channels

    !---------------------------------------------------------------------------
    !  [2]: Calculate statistics for cold-started predictors 
    !---------------------------------------------------------------------------
      global_cs = .true.
      do k = 1, iv%instid(inst)%nchan 
	 if (iv%instid(inst)%varbc(k)%npred <= 0) cycle                   !! VarBC channels only      
	 if (ANY(iv%instid(inst)%varbc(k)%pred_use > 0)) global_cs = .false.      
      end do	 
    
      if (global_cs) then                                                 !! Instrument coldstart only

         allocate (mean(npredmax), rms(npredmax), mean_tot(npredmax), rms_tot(npredmax))

        ! Accumulate statistics for predictor mean/std
        ! ---------------------------------------------
	 if (num_rad > 0) then
            do i = 1, npredmax
               mean(i) = SUM( iv%instid(inst)%varbc_info%pred(i,1:num_rad),    &
                         MASK=iv%instid(inst)%info%proc_domain(1,1:num_rad))   ! do not count HALO  
               rms(i)  = SUM( iv%instid(inst)%varbc_info%pred(i,1:num_rad)**2, &
                         MASK=iv%instid(inst)%info%proc_domain(1,1:num_rad))   ! do not count HALO 
            end do
	 else
	    mean = 0.0
	    rms  = 0.0
	 end if
  
        ! Do inter-processor communication to gather statistics
        ! ------------------------------------------------------
	 call wrf_dm_sum_reals(mean, mean_tot)
  	 call wrf_dm_sum_reals(rms,  rms_tot)	 
         if (num_rad_tot >= varbc_nobsmin) then
	    mean_tot = mean_tot / num_rad_tot
            rms_tot  = rms_tot  / num_rad_tot
	 else
	    mean_tot = 0.0
	    rms_tot  = 1.0   
	 end if
	 
        ! Store statistics
        !------------------
	 iv%instid(inst)%varbc_info%pred_mean = mean_tot
         iv%instid(inst)%varbc_info%pred_std  = sqrt(rms_tot - mean_tot**2)
      
         deallocate(mean, rms, mean_tot, rms_tot)

      end if
      deallocate(ipred)  	           	     	

    !---------------------------------------------------------------------------
    !  [3]: Normalize predictors
    !---------------------------------------------------------------------------
      do i = 1,  npredmax
         if ( iv%instid(inst)%varbc_info%pred_std(i) <= 0.0 ) cycle
         do n = 1, num_rad      
            iv%instid(inst)%varbc_info%pred(i,n) = &
          ( iv%instid(inst)%varbc_info%pred(i,n) - &
            iv%instid(inst)%varbc_info%pred_mean(i) ) / &
            iv%instid(inst)%varbc_info%pred_std(i)
	 end do     
      end do
      
   end do                           !  end loop for sensor
   
   if (trace_use) call da_trace_exit("da_varbc_coldstart")

 end subroutine da_varbc_coldstart
  subroutine da_varbc_precond (iv)

   !---------------------------------------------------------------------------
   !  PURPOSE:  Calculate covariance matrix b/w bias predictors 
   !            for VarBC preconditioning
   !
   !  Called from da_get_innov_vector_radiance
   !
   !  HISTORY: 11/07/2007 - Creation                     Tom Auligne
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type), intent (inout)   :: iv             ! Innovation

   integer                          :: inst, n, i, j, k, ii, jj
   integer                          :: npred, npredmax, num_rad, num_rad_active
   real, allocatable 		    :: hessian(:,:)
   real*8, allocatable 		    :: eignvec(:,:), eignval(:)
   real                             :: hessian_local, bgerr_local, pred_i, pred_j
   
   if ( iv%num_inst < 1 ) RETURN

   if (trace_use) call da_trace_entry("da_varbc_precond")
   
   write(unit=stdout,fmt='(A)') 'VARBC: Estimate Hessian for preconditioning'

   do inst = 1, iv%num_inst      ! loop for sensors
      npredmax = iv%instid(inst)%varbc_info%npredmax
 
      allocate ( hessian(npredmax, npredmax) )
      allocate ( eignvec(npredmax, npredmax) )
      allocate ( eignval(npredmax)           )

      do k = 1, iv%instid(inst)%nchan         ! loop for channels
         npred = iv%instid(inst)%varbc(k)%npred
         if (npred   <= 0) cycle              ! VarBC channel only     

         num_rad        = iv%instid(inst)%num_rad
         if (num_rad > 0) then
            num_rad_active = COUNT( (iv%instid(inst)%info%proc_domain(1,1:num_rad)) &
                               .AND.(iv%instid(inst)%tb_qc(k,1:num_rad) > qc_varbc_bad) )
         else
            num_rad_active = 0
         end if
         iv%instid(inst)%varbc(k)%nobs  = wrf_dm_sum_integer(num_rad_active)
	    
         if ( satinfo(inst)%iuse(k) == 1 ) &
            write(unit=stdout,fmt='(A,I6,3A,I5)') &
	    'VARBC:',iv%instid(inst)%varbc(k)%nobs,' active observations for ', &
            trim(adjustl(iv%instid(inst)%rttovid_string)),' channel',           &
	    iv%instid(inst)%ichan(k)
	    
         if (iv%instid(inst)%varbc(k)%nobs == 0) cycle
	 
        !---------------------------------------------------------	 
	! Calculate estimation of the Hessian for preconditioning
        !---------------------------------------------------------	 
         do i = 1, npred
	    ii = iv%instid(inst)%varbc(k)%ipred(i)

	   ! Observation term
           !------------------		 
            do j = i, npred
	       jj = iv%instid(inst)%varbc(k)%ipred(j)
	       hessian_local = 0.0	
	       
               do n= 1, num_rad      ! loop for pixel      
                  if (iv%instid(inst)%tb_qc(k,n) <= qc_varbc_bad)   cycle  ! good obs only
	          if (.not. iv%instid(inst)%info%proc_domain(1,n))  cycle  ! do not sum up HALO data
		  
		  if (ii == iv%instid(inst)%varbc_info%gammapred) then
		     pred_i = iv%instid(inst)%gamma_jacobian(k,n)
		  else
   	    	     pred_i = iv%instid(inst)%varbc_info%pred(ii,n)
		  end if   

		  if (jj == iv%instid(inst)%varbc_info%gammapred) then
		     pred_j = iv%instid(inst)%gamma_jacobian(k,n)
		  else
		     pred_j = iv%instid(inst)%varbc_info%pred(jj,n)		  
		  end if   
			      
       	          hessian_local = hessian_local + pred_i * pred_j / &
	                          iv%instid(inst)%tb_error(k,n)**2				     
               end do                                               !  end loop for pixel
	       
	       ! Sum hessian preconditioning across processors
               hessian(i,j) = wrf_dm_sum_real(hessian_local)	
	       hessian(j,i) = hessian(i,j)  	       
	    end do
         
	   ! Background term
           !-----------------
            if (iv%instid(inst)%varbc_info%nbgerr(ii) <= 0) cycle
	    bgerr_local = 0.0
	    do n= 1, num_rad    
               if (iv%instid(inst)%tb_qc(k,n) <= qc_varbc_bad)   cycle  ! good obs only
               if (.not. iv%instid(inst)%info%proc_domain(1,n))  cycle  ! do not sum up HALO data

   	       bgerr_local = bgerr_local + iv%instid(inst)%tb_error(k,n)**2  / &
                                             varbc_nbgerr
	                                   ! iv%instid(inst)%varbc_info%nbgerr(ii)
	    end do

            iv%instid(inst)%varbc(k)%bgerr(i) = wrf_dm_sum_real(bgerr_local) / &
	                                        iv%instid(inst)%varbc(k)%nobs

	    if (iv%instid(inst)%varbc(k)%bgerr(i) > 0) &          
   	       hessian(i,i) = hessian(i,i) + 1/iv%instid(inst)%varbc(k)%bgerr(i)
	 end do   

        !--------------------------------------------------	 
	! Preconditioning = inverse square root of Hessian
        !--------------------------------------------------	 	 
         hessian = hessian / varbc_factor**2
	 
	 if (npred == 1) then
	    iv%instid(inst)%varbc(k)%vtox(1,1)     = 1.0 / sqrt(hessian(1,1))
         else 
	    call da_eof_decomposition(npred, hessian(1:npred,1:npred), &
	    	            eignvec(1:npred,1:npred),eignval(1:npred))

	    if (ANY(eignval(1:npred) <= 0)) then		    
	       write(unit=stdout,fmt='(3A,I4,A,10F18.5)') &
	       'VARBC: non-positive Hessian for ', iv%instid(inst)%rttovid_string, &
	              ' ,channel ',k,'--> Eigenvalues =',eignval(1:npred) 
	       do i = 1, npred
	          if (hessian(i,i) > 0) &
	             iv%instid(inst)%varbc(k)%vtox(i,i) = 1.0 / sqrt(hessian(i,i))
	       end do
	    else
 	       do i = 1, npred
	          do j = i, npred
                     iv%instid(inst)%varbc(k)%vtox(i,j) = sum(          &
		                           eignvec(i,1:npred)         * &
			                   sqrt(1.0/eignval(1:npred)) * &
			                   eignvec(j,1:npred) )
		       
		     iv%instid(inst)%varbc(k)%vtox(j,i) = &
          	     iv%instid(inst)%varbc(k)%vtox(i,j)
                  end do
	       end do
	    end if
	 end if   
      end do                     !  end loop for channels	
      deallocate(hessian, eignvec, eignval)   
   end do                        !  end loop for sensor


   if (trace_use) call da_trace_exit("da_varbc_precond")

 end subroutine da_varbc_precond
  subroutine da_varbc_init (iv, be)

   !---------------------------------------------------------------------------
   !  PURPOSE: Initialize Variational bias correction from VARBC.in file
   !
   !  Called from da_radiance_init routine
   !
   !  HISTORY: 10/26/2007 - Creation                     Tom Auligne
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout)  :: iv       ! O-B structure.
 
   type (be_type), intent(inout) :: be        ! background error.

   !
   !  local arguments
   !------------------- 
   integer   :: n, i, ii, j, k, nchan, npred, npred_ws, ipred, ichan, np, nsensor
   integer   :: iunit, iost
   integer   :: size_jp, jp_start
   integer   :: platform_id, satellite_id, sensor_id, nchanl, npredmax, num_rad
   integer, allocatable :: pred_use(:), ipred_ws(:)
   real, allocatable    :: pred_ws(:)
   character(len=filename_len) :: filename
   character(len=120) :: cline
   logical   :: lvarbc_read, limatch, lmatch
   
   if (trace_use) call da_trace_entry("da_varbc_init")

   !--------------------------------------------------------------
   ! [1] Initializations
   !--------------------------------------------------------------
   size_jp  = 0                            !! Jp Control Variable size
   jp_start = be % cv % size_jb + be % cv % size_je
 
   do n = 1, iv % num_inst
      nchan = iv%instid(n)%nchan
      iv%instid(n)%varbc_info%npredmax = 0
      allocate ( iv%instid(n)%varbc(nchan) )
      do k = 1, nchan
         iv%instid(n)%varbc(k)%npred = 0   ! by default, do not use VarBC
         iv%instid(n)%varbc(k)%nobs  = 0 
      end do
   end do   

   !--------------------------------------------------------------
   ! [2] Read VarBC.in file and match with obs data set 
   !--------------------------------------------------------------
   filename = 'VARBC.in'
   inquire(file=trim(adjustl(filename)), exist=lvarbc_read)

   if (lvarbc_read) then
    write(unit=stdout,fmt='(A)') 'VARBC: Reading VARBC.in file'
    call da_get_unit(iunit)
    open(unit=iunit,file=filename, form='formatted', status='old')
	
    read(iunit, *) 
    read(iunit, *) nsensor
   
    do j = 1, nsensor
      read(iunit, *) ; read(iunit, *) ;read(iunit, *)
      read(iunit, *) platform_id, satellite_id, sensor_id, nchanl, npredmax
      limatch = .false.
       do n = 1, iv % num_inst
         if ( (platform_id  == iv%instid(n)%platform_id) .and. &   !! found match !!
              (satellite_id == iv%instid(n)%satellite_id) .and. &  !! found match !!
	      (sensor_id    == iv%instid(n)%sensor_id) ) then      !! found match !!	      
	    limatch = .true.                                       !! found match !!

!            write(unit=stdout,fmt='(A,A)') &
!     	       'VARBC: Matching with obs for ', iv%instid(n)%rttovid_string
	       
            num_rad  = iv%instid(n)%num_rad    
            allocate ( iv%instid(n)%varbc_info%pred(npredmax, num_rad) )	 
	    
            iv%instid(n)%varbc_info%platform_id  = platform_id
            iv%instid(n)%varbc_info%satellite_id = satellite_id
            iv%instid(n)%varbc_info%sensor_id    = sensor_id
            iv%instid(n)%varbc_info%nchanl       = nchanl
            iv%instid(n)%varbc_info%npredmax     = npredmax
            iv%instid(n)%varbc_info%gammapred    = 9    ! 1 Gamma Correction
            read(iunit, *)
            allocate ( iv%instid(n)%varbc_info%pred_mean(npredmax) )
            allocate ( iv%instid(n)%varbc_info%pred_std (npredmax) )
            allocate ( iv%instid(n)%varbc_info%nbgerr   (npredmax) )
            read(iunit, *) iv%instid(n)%varbc_info%pred_mean
            read(iunit, *) iv%instid(n)%varbc_info%pred_std
            read(iunit, *) iv%instid(n)%varbc_info%nbgerr
            read(iunit, *)

            allocate ( pred_use(npredmax), ipred_ws(npredmax), pred_ws(npredmax) )
            do i = 1, nchanl
	       lmatch = .false.
	       read(iunit, '(A)',iostat=iost) cline
               read(cline, *) ii, ichan, pred_use
               npred    = COUNT (pred_use >= 0)
               if (npred <= 0) cycle                !! VarBC channels only
	       
	       chan_loop: do k = 1, iv%instid(n)%nchan
	          if (iv%instid(n)%ichan(k) == ichan) then         !! found match !!		  
		     lmatch = .true.                               !! found match !!
	             iv%instid(n)%varbc(k)%ichanl = ichan
                     iv%instid(n)%varbc(k)%npred  = npred
	             allocate ( iv%instid(n)%varbc(k)%pred_use (npredmax) )
	             allocate ( iv%instid(n)%varbc(k)%ipred (npred) )
	             allocate ( iv%instid(n)%varbc(k)%index (npred) )
                     allocate ( iv%instid(n)%varbc(k)%param (npred) )
	             allocate ( iv%instid(n)%varbc(k)%bgerr (npred) )
                     allocate ( iv%instid(n)%varbc(k)%vtox    (npred,npred) )
		     iv%instid(n)%varbc(k)%pred_use = pred_use
                     iv%instid(n)%varbc(k)%vtox(1:npred, 1:npred)     = 0.0
		     do ipred = 1, npred
		        iv%instid(n)%varbc(k)%vtox(ipred, ipred)      = 1.0
                     end do
		     iv%instid(n)%varbc(k)%param(1:npred)             = 0.0
   	             iv%instid(n)%varbc(k)%bgerr(1:npred)             = 0.0
                     iv%instid(n)%varbc(k)%ipred(1:npred)             = &
                           PACK((/ (ipred, ipred=1,npredmax) /), mask = (pred_use >= 0))
		
		    ! VarBC warm-start parameters
   	      	    !-----------------------------
                     npred_ws = COUNT (pred_use > 0)
	             if (npred_ws > 0) then
		        read(cline, *) ii, ichan, pred_use, pred_ws(1:npred_ws)
 		        ipred_ws(1:npred_ws) = PACK((/ (ipred, ipred=1,npred) /), &
  			                       mask = (pred_use(iv%instid(n)%varbc(k)%ipred) >  0))
		        iv%instid(n)%varbc(k)%param(ipred_ws(1:npred_ws))= pred_ws(1:npred_ws)
                     end if
		     
                    ! Jp Control Variable size 
		    !--------------------------
                     do ipred = 1, npred
     	                 size_jp =  size_jp + 1
	                 iv%instid(n)%varbc(k)%index(ipred) = jp_start + size_jp
                     end do
		     
		     exit chan_loop
		  end if 
	       end do chan_loop
	       
	       if (.not. lmatch) write(unit=stdout,fmt='(A,I4)') &
	       			       'VARBC: no matching for channel:',ichan
            end do
            deallocate(pred_use, ipred_ws, pred_ws)   
         end if          
       end do
       if (.not. limatch) then
          read(iunit, *) ; read(iunit, *) ; read(iunit, *) ; read(iunit, *) ; read(iunit, *)
	  do i = 1, nchanl
	     read(iunit, *)
	  end do 
	  write(unit=stdout,fmt='(A,3I4)') &
	        'VARBC: no matching for platform/satid/sensor',platform_id, satellite_id, sensor_id
       end if
    end do
    close(iunit)
    call da_free_unit(iunit)
   else
      write(unit=stdout,fmt='(A)') 'VARBC: could not find VARBC.in file --> VARBC switched off'
      use_varbc    = .false.
      freeze_varbc = .false.
   end if 
   
   !--------------------------------------------------------------
   ! [3] Define VarBC control variable size:
   !--------------------------------------------------------------
   use_varbc = use_varbc.and.(.not.freeze_varbc)   
   if (use_varbc) then
      be % cv % size_jp = size_jp
      cv_size_domain_jp = size_jp
   end if
   
   if (trace_use) call da_trace_exit("da_varbc_init")

 end subroutine da_varbc_init 
  subroutine da_varbc_update (it, cv_size, cv, iv)

   !---------------------------------------------------------------------------
   !  PURPOSE: Update VarBC parameters and write into file
   !
   ! [1] Update of VarBC parameters (including change of variable)
   !
   ! [2] Write VARBC.out file with VarBC information
   !
   !  Called from da_solve
   !
   !  HISTORY: 10/29/2007 - Creation                          Tom Auligne
   !           01/20/2010 - Update for multiple outer-loops   Tom Auligne
   !---------------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: it          !outer loop counting
   integer,        intent(in)    :: cv_size
   real,           intent(in)    :: cv(cv_size) ! Control variable structure.
   type (iv_type), intent(inout) :: iv          ! Obs. increment structure.

   !
   !  local arguments
   !------------------- 
   integer              :: inst, ichan, npred, i, npredmax, id
   integer              :: iunit, iost
   character(len=filename_len) :: filename
   character(len=120)   :: fparam
   real, allocatable    :: varbc_param_tl(:)
   logical, allocatable :: lvarbc_inst(:)

   if (trace_use) call da_trace_entry("da_varbc_update")

   write(unit=stdout,fmt='(A)') 'VARBC: Updating bias parameters'
      
   allocate( lvarbc_inst(iv%num_inst) )
   do inst = 1, iv % num_inst   
      lvarbc_inst(inst) = ANY(iv%instid(inst)%varbc(1:iv%instid(inst)%nchan)%npred>0)	
   end do   
   
   do inst = 1, iv % num_inst   
      if (.not. lvarbc_inst(inst)) cycle             !! VarBC instruments only   
      allocate(varbc_param_tl(iv%instid(inst)%varbc_info%npredmax))
      
      do ichan = 1, iv%instid(inst)%nchan
         npred    = iv%instid(inst)%varbc(ichan)%npred
         if (npred <= 0) cycle               !! VarBC channels only	 
	 if (iv%instid(inst)%varbc(ichan)%nobs >= varbc_nobsmin) then
	    where (iv%instid(inst)%varbc(ichan)%pred_use == 0) &
	           iv%instid(inst)%varbc(ichan)%pred_use = 1
	 else
            if ( satinfo(inst)%iuse(ichan) == 1) then
	       if (count(iv%instid(inst)%varbc(ichan)%pred_use == 0) > 0) &
	          write(unit=stdout,fmt='(A,A,I5)') &	   
	          'VARBC: Not enough data to keep statistics for ', &
	          trim(iv%instid(inst)%rttovid_string),             &
		  iv%instid(inst)%ichan(ichan)
	    end if
	 end if
	 
 	 if (.not. use_varbc) cycle
         
	!---------------------------------------------------------------
        ! Change of variable (preconditioning) for parameters increments
        !---------------------------------------------------------------
            
	 do i = 1, npred
            id = iv%instid(inst)%varbc(ichan)%index(i)
    	    varbc_param_tl(i) = &
	       SUM(cv(id) * iv%instid(inst)%varbc(ichan)%vtox(i,1:npred))
         end do	
	    
        !---------------------------------------------------------------
        ! Update VarBC parameters
        !---------------------------------------------------------------
	 iv%instid(inst)%varbc(ichan)%param = &
	 iv%instid(inst)%varbc(ichan)%param + varbc_param_tl(1:npred)
      end do
      
      deallocate(varbc_param_tl)
   end do
   
   if (.not. rootproc) then
      if (trace_use) call da_trace_exit("da_varbc_update")
      return
   end if

   !---------------------------------------------------------------
   ! Write VARBC.out file
   !---------------------------------------------------------------
   
   write(unit=stdout,fmt='(A)') 'VARBC: Writing information in VARBC.out file'

     
   call da_get_unit(iunit)
   if ( it == max_ext_its ) then
      filename = 'VARBC.out'
   else
      write(unit=filename, fmt='(a,i2.2)') 'VARBC.out_',it
   end if
   open(unit=iunit,file=filename,form='formatted',iostat = iost,status='replace')
	
   if (iost /= 0) then
      message(1)="Cannot open bias correction file "//adjustl(filename)
      call da_error("da_varbc_update.inc",106,message(1:1))
   end if

   write(iunit, *) ' VARBC version 1.0 - Number of instruments: '
   write(iunit, *) COUNT(lvarbc_inst)
   
   do inst = 1, iv % num_inst   
      if (.not. lvarbc_inst(inst)) cycle             !! VarBC instruments only   
      npredmax = iv%instid(inst)%varbc_info%npredmax

      write(iunit, *) '------------------------------------------------'
      write(iunit, *) 'Platform_id  Sat_id  Sensor_id  Nchanl  Npredmax'
      write(iunit, *) '------------------------------------------------'
      write(iunit, *) iv%instid(inst)%varbc_info%platform_id, &
                      iv%instid(inst)%varbc_info%satellite_id, &
		      iv%instid(inst)%varbc_info%sensor_id, &
      		      iv%instid(inst)%varbc_info%nchanl,&
		      iv%instid(inst)%varbc_info%npredmax

      write(iunit, *) ' -----> Bias predictor statistics:  Mean & Std & Nbgerr'
      write(iunit,'(30F10.1)') iv%instid(inst)%varbc_info%pred_mean(1:npredmax)
      write(iunit,'(30F10.1)') iv%instid(inst)%varbc_info%pred_std (1:npredmax)
      write(iunit,'(30I10)')   iv%instid(inst)%varbc_info%nbgerr   (1:npredmax)
            
      write(iunit, *) ' -----> Chanl_id Chanl_nb  Pred_use(-1/0/1)  Param'
      
      do ichan = 1, iv%instid(inst)%nchan
         npred    = iv%instid(inst)%varbc(ichan)%npred
         if (npred <= 0) cycle               !! VarBC channels only	 
	 write(fparam,*) '(I4,I6,',npredmax,'I3,',npred,'F8.3)'
         write(iunit, fmt=trim(adjustl(fparam))) &
	     ichan, iv%instid(inst)%varbc(ichan)%ichanl,   &
	            iv%instid(inst)%varbc(ichan)%pred_use, &
  	            iv%instid(inst)%varbc(ichan)%param
      end do
   end do

   deallocate(lvarbc_inst)
   close(iunit)
   call da_free_unit(iunit)

   if (trace_use) call da_trace_exit("da_varbc_update")


 end subroutine da_varbc_update


end module da_varbc
