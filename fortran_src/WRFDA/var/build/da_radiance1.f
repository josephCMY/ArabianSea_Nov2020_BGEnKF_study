












module da_radiance1

   
   
   

   use module_radiance, only : satinfo,q2ppmv,rttov_inst_name
   use module_radiance, only : CRTM_Planck_Radiance, CRTM_Planck_Temperature

   use da_control, only : trace_use,missing_r, rootproc, &
      stdout,myproc,qc_good,num_fgat_time,qc_bad, &
      use_error_factor_rad,biasprep_unit,obs_qc_pointer, filename_len, &
      print_detail_rad, rtm_option, trace_use_dull, &
      rtm_option_rttov,rtm_option_crtm, radiance, only_sea_rad, &
      global, gas_constant, gravity, monitor_on,kts,kte,use_rttov_kmatrix, &
      use_pseudo_rad, pi, t_triple, crtm_cloud, DT_cloud_model,write_jacobian, &
      use_crtm_kmatrix,use_clddet_mmr, use_satcv, cv_size_domain, &
      cv_size_domain_js, calc_weightfunc, use_clddet_ecmwf, deg_to_rad, rad_to_deg
   use da_define_structures, only : info_type,model_loc_type,maxmin_type, &
      iv_type, y_type, jo_type,bad_data_type,bad_data_type,number_type, &
      be_type
   use module_dm, only : wrf_dm_sum_real, wrf_dm_sum_integer
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int,da_proc_sum_ints
   use da_reporting, only : da_error, message
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_residual_new, da_eof_decomposition
   use da_tools_serial, only : da_free_unit, da_get_unit
   use da_tracing, only : da_trace_entry, da_trace_exit, da_trace_int_sort

   use da_control, only : rtminit_sensor,write_profile,num_procs,tovs_min_transfer
   use da_reporting, only : da_warning, da_message
   use da_tracing, only : da_trace

   implicit none
   
   type datalink_type

      type (info_type)        :: info
      type (model_loc_type)   :: loc

      integer   ::  ifgat, landsea_mask, rain_flag
      integer   ::  scanline, scanpos
      real      ::  satzen, satazi, solzen, solazi  
      
      real,    pointer   ::   emiss(:)
      
      integer, pointer   ::   cloud_flag(:)
      real,    pointer   ::   t(:), mr(:), zk(:)
      real,    pointer   ::   pm(:), tm(:), qm(:), qrn(:), qcw(:),qci(:),qsn(:),qgr(:)
      real               ::   ps,ts,t2m,mr2m,u10,v10, clwp
      real               ::   smois, tslb, snowh, elevation,soiltyp,vegtyp,vegfra
      real               ::   clw
      integer            ::   isflg

      real, pointer             :: tb_ob(:)
      real, pointer             :: tb_inv(:)
      real, pointer             :: tb_qc(:)
      real, pointer             :: tb_error(:)
      integer                   :: sensor_index
      type (datalink_type), pointer  :: next 
   end type datalink_type

   type con_vars_type
      integer            ::  nlevels
      real   ,  pointer  ::  t(:)
      real   ,  pointer  ::  q(:)
      real               ::  ps
      real   ,  pointer  ::  t_jac(:,:) => null()
      real   ,  pointer  ::  q_jac(:,:) => null()
      real   ,  pointer  ::  ps_jac(:)  => null()
   end type con_vars_type

   type con_cld_vars_type
      integer            ::  nwp_levels
      real   ,  pointer  ::  p(:)
      real   ,  pointer  ::  ph(:)
      real   ,  pointer  ::  t(:)
      real   ,  pointer  ::  cc(:)
      real   ,  pointer  ::  clw(:)   
      real   ,  pointer  ::  ciw(:)   
      real   ,  pointer  ::  rain(:)  
      real   ,  pointer  ::  sp(:)    
   end type con_cld_vars_type

   type aux_vars_type
      integer            ::  surftype
      real               ::  surft, t2m, q2m, u10, v10
      real               ::  satzen, satazi  
      real               ::  solzen, solazi
      real               ::  elevation ,rlat
   end type aux_vars_type

   type maxmin_rad_stats_type
      type (maxmin_type)         :: maximum, minimum
      real                       :: ave, rms
      integer                    :: num
   end type maxmin_rad_stats_type

   type stats_rad_type
      type (maxmin_rad_stats_type), pointer  :: ichan(:)
   end type stats_rad_type

   type rad_header_type                       
      character (LEN = 19) :: date_char       
      integer              :: assim_win       
      character(LEN=20)    :: rttovid_string  
      integer              :: platform_id     
      integer              :: satellite_id    
      integer              :: sensor_id       
      integer              :: num_rad         
      integer              :: nchan           
      integer ,  pointer   :: ichan(:)        
      integer              :: nemis           
                                              
                                              
      integer              :: nlevel_fix      
                                              
      real   ,   pointer   :: pres(:)         
      integer              :: nlevel_cld      
   end type rad_header_type

   type rad_data_type                       

      
      integer            :: landmask      
      integer            :: scanline      
      integer            :: scanpos       
      real               :: lat           
      real               :: lon           
      real               :: elv           
      real               :: satzen        
      real               :: satazi        
      real               :: solzen        
      real               :: solazi        
      real,    pointer   :: tb(:)         
      real,    pointer   :: inv(:)        
      real,    pointer   :: bias(:)       
      real,    pointer   :: err(:)        
      real,    pointer   :: qc(:)         
                                          
      real,    pointer   :: emiss(:)      

      
      integer            :: surftype      
                                          
                                          
      integer            :: terrain       
      integer            :: soiltyp       
      integer            :: vegtyp        
      real               :: vegfra        
      real               :: soilm         
      real               :: soilt         
      real               :: snowh         
      real               :: ps            
      real               :: ts            
      real               :: t2m           
      real               :: mr2m          
      real               :: u10,v10       
      real,    pointer   :: t(:)          
      real,    pointer   :: mr(:)         
      real,    pointer   :: zk(:)         
      real,    pointer   :: pm(:)         
      real,    pointer   :: phm(:)        
      real,    pointer   :: tm(:)         
      real,    pointer   :: cc(:)         
      real,    pointer   :: rain(:)       
      real,    pointer   :: solidp(:)     
      real,    pointer   :: clw(:)        
      real,    pointer   :: ciw(:)        

   end type rad_data_type

   type bias_type
      integer :: nchan     
      integer :: npred     
      integer :: platform_id,satellite_id,sensor_id
      integer :: year, month, day, hour, min, sec
      integer :: scanline,scanpos
      integer :: landmask
      integer, pointer :: qc_flag(:) 
      integer, pointer :: cloud_flag(:) 
      integer :: surf_flag  
      real    :: elevation,lat,lon,ps, t2m, q2m, tsk, clwp
      real, pointer  :: tb(:), omb(:), bias(:)
      real, pointer  :: pred(:)
   end type bias_type

   integer, allocatable :: num_tovs_before(:,:)
   integer, allocatable :: num_tovs_after(:,:)
   integer, allocatable :: tovs_send_pe(:,:)
   integer, allocatable :: tovs_send_start(:,:)
   integer, allocatable :: tovs_send_count(:,:)
   integer, allocatable :: tovs_recv_pe(:,:)
   integer, allocatable :: tovs_recv_start(:,:)
   integer, allocatable :: tovs_copy_count(:)

contains

subroutine da_jo_and_grady_rad(iv, re, jo, jo_grad_y) 

   !---------------------------------------------------------------------------
   ! Purpose: Calculate Gradient_y i and cost function Jo for radiance data.
   !
   ! Method:  grad_y = -R^-1 (d - H delta_x)
   !              Jo = -(d - H delta_x) grad_y
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type) , intent(in)    :: re          ! Residual vector.
   type (y_type) , intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo          ! Obs cost function.

   integer                       :: n, k, i

   if (trace_use) call da_trace_entry("da_jo_and_grady_rad")

   do i =1, iv%num_inst

      jo % rad(i)%jo_ichan(:) = 0.0
      jo % rad(i)%num_ichan(:) = 0

      if (iv%instid(i)%num_rad < 1 .or. iv%instid(i)%rad_monitoring == monitor_on ) cycle

      do n=1, iv%instid(i)%num_rad
         do k=1, iv%instid(i)%nchan
            jo_grad_y%instid(i)%tb(k,n) = -re%instid(i)%tb(k,n) / &
               (iv%instid(i)%tb_error(k,n) * iv%instid(i)%tb_error(k,n))
         end do
         if (iv%instid(i)%info%proc_domain(1,n)) then
            do k=1, iv%instid(i)%nchan
               if (iv%instid(i)%tb_qc(k,n) >= obs_qc_pointer) then
                  jo % rad(i) % jo_ichan(k) = jo % rad(i) % jo_ichan(k) - &
                     re%instid(i)%tb(k,n) * jo_grad_y%instid(i)%tb(k,n)
                  jo % rad(i) % num_ichan(k) = jo % rad(i) % num_ichan(k) + 1
               end if
            end do
         end if
      end do
      jo % rad(i)%jo_ichan(:) = 0.5 * jo % rad(i)%jo_ichan(:)
   end do

   if (trace_use) call da_trace_exit("da_jo_and_grady_rad")

end subroutine da_jo_and_grady_rad


subroutine da_residual_rad(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !---------------------------------------------------------------------------
   ! Purpose: Calculate Obs Residual and counting obs number.
   !
   ! Method:  re = (d - H delta_x)
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv     ! Innovation vector (O-B).
   type (y_type) , intent(in)    :: y      ! y = H (xa)
   type (y_type) , intent(inout) :: re     ! Residual vector (O-A).

   integer       , intent(inout) :: np_available
   integer       , intent(inout) :: np_obs_used
   integer       , intent(inout) :: np_missing
   integer       , intent(inout) :: np_bad_data

   type (bad_data_type)              :: n_obs_bad
   integer                           :: i

   if (trace_use) call da_trace_entry("da_residual_rad")

   do i = 1, iv%num_inst
      if (iv%instid(i)%num_rad < 1) cycle

      n_obs_bad % rad % num = number_type(0, 0, 0)
      call da_residual_new(y%instid(i)%tb(:,:), iv%instid(i)%tb_qc(:,:), &
         iv%instid(i)%tb_inv(:,:), re%instid(i)%tb(:,:))
      np_available = iv%instid(i)%nchan*iv%instid(i)%num_rad
   end do

   np_missing  = np_missing  + n_obs_bad % rad % num % miss
   np_bad_data = np_bad_data + n_obs_bad % rad % num % bad
   np_obs_used = np_obs_used + n_obs_bad % rad % num % use 

   if (trace_use) call da_trace_exit("da_residual_rad")
   
end subroutine da_residual_rad


subroutine da_biascorr ( i, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform bias correction for radiance data.
   !
   ! METHOD:  omb(corrected)=omb-scanbias-airbias
   !---------------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: i       ! sensor index.
   type (y_type),  intent(in)    :: ob      ! Observation structure.
   type (iv_type), intent(inout) :: iv      ! O-B structure.

   ! Local variables
   integer   :: k,iband,iscan, n,j,npred,nlevels, num_rad
   real      :: pred(4),airbias
   real,allocatable      :: q(:), temp(:), hum(:), pf(:)

   num_rad = iv%instid(i)%info%n2-iv%instid(i)%info%n1+1
   if (num_rad < 1) return

   if (trace_use) call da_trace_entry("da_biascorr")

   npred=4
   nlevels=iv%instid(i)%nlevels-1
   allocate(temp(nlevels))
   allocate(hum(nlevels))
   allocate(pf(0:nlevels))

   do n=iv%instid(i)%info%n1,iv%instid(i)%info%n2
      ! get airmass predictors
      !-------------------------
      if (rtm_option==rtm_option_rttov) then
      else if (rtm_option==rtm_option_crtm) then
! FIX? problems with IBM AIX COMPILER
         temp(1:nlevels) = iv%instid(i)%tm(1:nlevels,n)
         hum(1:nlevels) = iv%instid(i)%qm(1:nlevels,n)
         pf(0:nlevels) = iv%instid(i)%pf(0:nlevels,n)
         call da_predictor_crtm(pred(1:npred), npred, nlevels,temp, &
            hum, iv%instid(i)%ts(n), pf)
      end if
      iscan = iv%instid(i)%scanpos(n)
      iband = floor(iv%instid(i)%info%lat(1,n)/10.0001) + 10
      do k=1,iv%instid(i)%nchan
         ! only do bias correction for used channels
         if ( (satinfo(i)%iuse(k) == 1) .and. (iv%instid(i)%tb_inv(k,n) > missing_r) ) then
            ! scan bias correction
            !-----------------------
            if (global) then
               iv%instid(i)%tb_inv(k,n) = iv%instid(i)%tb_inv(k,n) - satinfo(i)%scanbias_b(k,iscan,iband)
            else
               iv%instid(i)%tb_inv(k,n) = iv%instid(i)%tb_inv(k,n) - satinfo(i)%scanbias(k,iscan) 
            end if
            ! airmass bias correction
            !----------------------------
            airbias = satinfo(i)%bcoef0(k)
            do j=1,npred
               airbias= airbias + satinfo(i)%bcoef(k,j)*pred(j)
            end do
            iv%instid(i)%tb_inv(k,n) = iv%instid(i)%tb_inv(k,n)-airbias
         end if
      end do
   end do

   deallocate(q)
   deallocate(temp)
   deallocate(hum)
   deallocate(pf)

   if (trace_use) call da_trace_exit("da_biascorr")

end subroutine da_biascorr

subroutine da_read_biascoef(string,nchan,nscan,nband,npred,global, &
      scanbias,scanbias_b,coef,coef0, vstd_dep)

   integer,           intent(in)    :: nchan,nscan,nband,npred
   logical,           intent(in)    :: global
   character(len=20), intent(in)    :: string
   real,              intent(inout) :: scanbias(nchan,nscan)
   real,              intent(inout) :: scanbias_b(nchan,nscan,nband)
   real,              intent(inout) :: coef(nchan,npred),coef0(nchan)
   real,              intent(inout) :: vstd_dep(nchan)

   integer :: iunit,j,i, ii,jj, iost
   integer :: nobs(nchan)
   real    :: vmean_abs(nchan)
   real    :: vstd_abs(nchan)
   real    :: vmean_dep(nchan)

   character(len=80) :: filename

   call da_trace_entry("da_read_biascoef")

   call da_get_unit(iunit)
   filename='biascorr/'//trim(adjustl(string))//'.bcor'
   open(unit=iunit,file=filename, form='formatted',iostat = iost, status='old')
   if (iost /= 0) then
      message(1)="Cannot open radiance biascorr file "//adjustl(filename)
      call da_error("da_read_biascoef.inc",27,message(1:1))
   end if

   ! read (iunit,'(4i6)') nchan,nscan,nband,npred
   read (iunit,'(4i6)')
   do i=1, nchan
      read (iunit,'(i5,i10,4F8.2)') ii,nobs(i),vmean_abs(i),vstd_abs(i), vmean_dep(i),vstd_dep(i)
   end do

   do i=1, nchan
      read (iunit,'(i5,5F12.5)') ii,(coef(i,j),j=1,npred),coef0(i)
   end do

   read (iunit,*)
   read (iunit,*)
   if (global) then   ! global coefs not available now, use regional one
      do j=1, nchan
         read(iunit,'(i5,90F7.2)') jj, scanbias(j,1:nscan)
      end do
   else
      do j=1, nchan
         read(iunit,'(i5,90F7.2)') jj, scanbias(j,1:nscan)
      end do
   end if

   close(iunit)
   call da_free_unit(iunit)

   call da_trace_exit("da_read_biascoef")

end subroutine da_read_biascoef


subroutine da_biasprep(inst,ob,iv)

   !-----------------------------------------------------------------------
   ! Purpose: Output information files for bias correction progs
   !-----------------------------------------------------------------------

   implicit none

   integer       ,  intent(in)      :: inst
   type (y_type) ,  intent(in)      :: ob         ! O structure.
   type (iv_type),  intent(in)      :: iv         ! O-B structure.

   integer  :: n,jx,npred,nchan,num_rad,nlevels
   character(len=80)  :: filename
   character(len=1)   :: s1
   real               :: pred(6)
   type (bias_type)   :: radbias
   real,allocatable      :: q(:), temp(:), hum(:), pf(:)

   num_rad = iv%instid(inst)%info%n2-iv%instid(inst)%info%n1+1

   if (num_rad < 1) return

   if (trace_use) call da_trace_entry("da_biasprep")

   write(filename, '(a,i4.4)') 'biasprep_'//trim(iv%instid(inst)%rttovid_string)//'.', myproc

   call da_get_unit(biasprep_unit)
   open(unit=biasprep_unit,FILE=filename,FORM='unformatted')

   !---------------------------------------------------------------------------
   npred = 4
   nchan = iv%instid(inst)%nchan 
   nlevels = iv%instid(inst)%nlevels-1
   allocate(q(nlevels))
   allocate(temp(nlevels))
   allocate(hum(nlevels))
   allocate(pf(0:nlevels))

   allocate(radbias%tb(nchan))
   allocate(radbias%omb(nchan))
   allocate(radbias%bias(nchan))
   allocate(radbias%qc_flag(nchan))
   allocate(radbias%cloud_flag(nchan))
   allocate(radbias%pred(npred))

   do n=iv%instid(inst)%info%n1,iv%instid(inst)%info%n2
      if (iv%instid(inst)%info%proc_domain(1,n)) then 

        if (rtm_option==rtm_option_rttov) then
        else if (rtm_option==rtm_option_crtm) then
! FIX? problems with IBM AIX COMPILER
         temp(1:nlevels) = iv%instid(inst)%tm(1:nlevels,n)
         hum(1:nlevels) = iv%instid(inst)%qm(1:nlevels,n)
         pf(0:nlevels) = iv%instid(inst)%pf(0:nlevels,n)
         call da_predictor_crtm(pred(1:npred), npred, nlevels,temp, &
            hum, iv%instid(inst)%ts(n), pf)
        end if

         ! transfer information to bias structure
         radbias%platform_id  = iv%instid(inst)%platform_id
         radbias%satellite_id = iv%instid(inst)%satellite_id
         radbias%sensor_id    = iv%instid(inst)%sensor_id

         read(iv%instid(inst)%info%date_char(n),'(i4,5(a1,i2))') &
                                   radbias%year,s1, radbias%month,s1, radbias%day, &
                                   s1,radbias%hour, s1,radbias%min, s1,radbias%sec

         radbias%scanline     = iv%instid(inst)%scanline(n)    ! not available
         radbias%scanpos      = iv%instid(inst)%scanpos(n)
         radbias%landmask     = iv%instid(inst)%landsea_mask(n)
         radbias%elevation    = iv%instid(inst)%info%elv(n)
         radbias%lat          = iv%instid(inst)%info%lat(1,n)
         radbias%lon          = iv%instid(inst)%info%lon(1,n)
         radbias%surf_flag    = iv%instid(inst)%isflg(n)
         radbias%ps           = iv%instid(inst)%ps(n)
         radbias%t2m          = iv%instid(inst)%t2m(n)
         radbias%q2m          = iv%instid(inst)%mr2m(n)/q2ppmv
         radbias%tsk          = iv%instid(inst)%ts(n)
         radbias%clwp         = iv%instid(inst)%clwp(n)  ! in mm

         radbias%nchan        = nchan 
         radbias%tb(1:nchan)  = ob%instid(inst)%tb(1:nchan,n)
         radbias%omb(1:nchan) = ob%instid(inst)%tb(1:nchan,n)-iv%instid(inst)%tb_xb(1:nchan,n)
         radbias%bias(1:nchan) = 0.0

         radbias%npred         = npred
         radbias%pred(1:npred) = pred(1:npred)

         radbias%qc_flag(1:nchan)= iv%instid(inst)%tb_qc(1:nchan,n)
         radbias%cloud_flag(1:nchan)= iv%instid(inst)%cloud_flag(1:nchan,n)

         ! set missing data and bad data to missing
         do jx=1,nchan   
            if (radbias%tb(jx) < 150.0 .or. radbias%tb(jx) > 400.0 ) then
               radbias%tb(jx)   = missing_r
               radbias%omb(jx)  = missing_r 
            end if
         end do

         !write(unit=biasprep_unit) radbias ! can not compiled with pointer

         call da_write_biasprep(radbias)

      end if
   end do
  
   close(unit=biasprep_unit)
   call da_free_unit(biasprep_unit)

   deallocate(q)
   deallocate(temp)
   deallocate(hum)
   deallocate(pf)

   deallocate(radbias%tb)
   deallocate(radbias%omb)
   deallocate(radbias%bias)
   deallocate(radbias%qc_flag)
   deallocate(radbias%cloud_flag)
   deallocate(radbias%pred)

   if (trace_use) call da_trace_exit("da_biasprep")

end subroutine da_biasprep


subroutine da_write_biasprep(radbias)

   implicit none

   type (bias_type), intent(in)  :: radbias

   if (trace_use) call da_trace_entry("da_write_biasprep")

      write(unit=biasprep_unit) radbias%nchan,radbias%npred
      write(unit=biasprep_unit) radbias%platform_id , &
                                radbias%satellite_id, &
                                radbias%sensor_id,    &
                                radbias%year,radbias%month,&
                                radbias%day, radbias%hour, &
                                radbias%min, radbias%sec,  &
                                radbias%scanline, &
                                radbias%scanpos,  &
                                radbias%landmask, &
                                radbias%elevation,&
                                radbias%lat,radbias%lon, &
                                radbias%ps,radbias%t2m, &
                                radbias%q2m,radbias%tsk, &
                                radbias%tb(1:radbias%nchan), &
                                radbias%omb(1:radbias%nchan), &
                                radbias%bias(1:radbias%nchan), &
                                radbias%pred(1:radbias%npred), &
                                radbias%qc_flag(1:radbias%nchan), &
                                radbias%cloud_flag(1:radbias%nchan), &
                                radbias%surf_flag, radbias%clwp

   if (trace_use) call da_trace_exit("da_write_biasprep")

end subroutine da_write_biasprep

subroutine da_predictor_crtm(pred,npred,nlevels,temp,hum,t_skin,pf)

   implicit none

   ! temp - model level temperatures (k)
   ! hum  - model level moistures    (g/kg)
   ! t-skin - model skin temperature (k)
   ! nlevels - number of model levels (0:model top)
   ! pm   - model level pressure (hpa)
   ! pf   - full level pressure  (hpa)

   ! pred(1) - 1000-300 thickness
   ! pred(2) - 200-50 thickness
   ! pred(3) - t_skin
   ! pred(4) - total column precipitable water

   integer, intent(in)  :: npred,nlevels
   real,    intent(in)  :: temp(nlevels), hum(nlevels), t_skin
   real,    intent(in)  :: pf(0:nlevels)
   real,    intent(out) :: pred(npred)

   real, parameter :: kth = gas_constant/gravity
   real, parameter :: kpc = 100.0/gravity
   real :: tv(nlevels), qm(nlevels)
   real :: dlp(nlevels), dp(nlevels)
   real    :: add_thk
   integer :: itmp(1)
   integer :: index1000, index300, index200, index50

   if (trace_use) call da_trace_entry("da_predictor_crtm")

   qm=hum*0.001  ! g/kg to kg/kg

   dlp(1:nlevels) = log(pf(1:nlevels)) - log(pf(0:nlevels-1))
    dp(1:nlevels) =     pf(1:nlevels)  -     pf(0:nlevels-1)

   ! 0.0 find the pressure level index that is closest to 
   ! 1000hPa, 300hPa, 200hPa, and 50hPa respectively

   ! note: pf levels are 0:nlevels, 
   ! return values of minloc are 1:nlevels+1
   itmp(1:1) = minloc(abs(pf(:)-1000.0))
   index1000 = itmp(1)-1
   itmp(1:1) = minloc(abs(pf(:)-300.0))
   index300  = itmp(1)-1
   itmp(1:1) = minloc(abs(pf(:)-200.0))
   index200  = itmp(1)-1
   itmp(1:1) = minloc(abs(pf(:)-50.0))
   index50   = itmp(1)-1

   ! 1.0 convert all temperatures to virtual temperatures
   ! ----------------------------------------------------
   tv = temp*(1.0+0.608*qm)

   if ( (1000.0 > pf(nlevels)) ) then
      add_thk = kth*( tv(nlevels)*(log(1000.0)-log(pf(nlevels))) )  ! approximation
   else
      add_thk = 0.0
   end if

   ! 2.0 construct averages for nesdis thick layers
   ! ----------------------------------------------

   pred(1) = kth*sum( tv(index300+1:index1000)*dlp(index300+1:index1000) ) + add_thk
   pred(2) = kth*sum( tv(index50+1:index200)*dlp(index50+1:index200) )
   pred(3) = t_skin
   pred(4) = kpc*sum( qm(1:nlevels)*dp(1:nlevels) )

   if (trace_use) call da_trace_exit("da_predictor_crtm")

end subroutine da_predictor_crtm

subroutine da_qc_crtm (it, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for radiance data.
   !
   ! METHOD:  seperated QC for each sensor
   !---------------------------------------------------------------------------

   implicit none

   integer      ,  intent(in)      :: it
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.

   integer :: i, nchan
   logical   :: amsua, amsub, hirs, msu,airs,hsb, ssmis, mhs

   if (trace_use) call da_trace_entry("da_qc_crtm")

   do i = 1, iv%num_inst

      !if (iv%instid(i)%info%n2 < iv%instid(i)%info%n1) cycle

      nchan    = iv%instid(i)%nchan

      amsua = trim(rttov_inst_name(rtminit_sensor(i))) == 'amsua'
      amsub = trim(rttov_inst_name(rtminit_sensor(i))) == 'amsub'
      hirs  = trim(rttov_inst_name(rtminit_sensor(i))) == 'hirs'
      msu   = trim(rttov_inst_name(rtminit_sensor(i))) == 'msu'
      airs  = trim(rttov_inst_name(rtminit_sensor(i))) == 'airs'
      hsb   = trim(rttov_inst_name(rtminit_sensor(i))) == 'hsb'
      ssmis = trim(rttov_inst_name(rtminit_sensor(i))) == 'ssmis'
      mhs   = trim(rttov_inst_name(rtminit_sensor(i))) == 'mhs'

      if (hirs) then
         call da_qc_hirs(it, i,nchan,ob,iv)
      else if (airs) then
         call da_qc_airs(it, i,nchan,ob,iv)
      else if ( hsb ) then
         ! call da_qc_hsb(it, i,nchan,ob,iv)
         call da_warning("da_qc_crtm.inc",42,(/'QC Not implemented for HSB'/))
      else if (amsua) then
         call da_qc_amsua(it,i,nchan,ob,iv)
      else if ( amsub ) then
         call da_qc_amsub(it,i,nchan,ob,iv)
      else if (msu) then
         ! call da_qc_msu(it, i,nchan, ob,iv)
         call da_warning("da_qc_crtm.inc",49,(/'QC Not implemented for MSU'/))
      else if (ssmis) then
         call da_qc_ssmis(it, i,nchan,ob,iv)
      else if (mhs) then
         call da_qc_mhs(it,i,nchan,ob,iv)
      else
         write(unit=message(1),fmt='(A,A)') &
            "Unrecognized instrument",trim(rttov_inst_name(rtminit_sensor(i)))
         call da_error("da_qc_crtm.inc",57,message(1:1))
      end if

   end do

   if (trace_use) call da_trace_exit("da_qc_crtm")

end subroutine da_qc_crtm


subroutine da_cloud_sim(KINDIC,KDIM,PX,PF,PG,IZS,RZS,DZS)

! Purpose :
! -------
! Simulate the cloud as a linear combination of grey clouds at model levels

! Interface :
! ---------
! KINDIC
! KDIM           : Dimension of cloud fraction variable
! PX             : Cloud fraction variable             -> Input
! PF             : Fitness Error
! PG             : Gradient of cloud fraction variable -> Output

! Externals :
! ---------

! Method :
! ------

! Reference :
! ---------

! Author :
! ------
! 01/08/2005  Thomas Auligne         *ECMWF*

! Modifications :
! -------------

! ---------------------------------------------------------------------------------------------

IMPLICIT NONE

!! Parameters !!
INTEGER,INTENT(IN)      :: KINDIC 
INTEGER,INTENT(IN)      :: KDIM 
INTEGER,INTENT(IN)      :: IZS(2)
double precision   ,INTENT(INOUT)   :: PX(KDIM) !izs(2)
double precision   ,INTENT(OUT)     :: PF 
double precision   ,INTENT(OUT)     :: PG(KDIM)
real               ,INTENT(IN)      :: RZS(kdim*izs(2))      ! Eigenvectors
DOUBLE PRECISION   ,INTENT(IN)      :: DZS(IZS(1)*KDIM) ! AMAT
!! Local arrays !!
INTEGER                 :: JCH, ilev, JLEV, nchan, neignvec
REAL                    :: ZNORM_PG, ZCLR, ZDCLR, eignvec(kdim,izs(2)), eignval(izs(2))
double precision        :: AMAT(IZS(1),KDIM)
double precision        :: alpha, beta
double precision        :: zx(KDIM), zgx(KDIM, KDIM), zx_eof(KDIM)

!IF (KINDIC == 1) RETURN
 PF       = 0.0
 PG       = 0.0
 nchan    = izs(1)
 neignvec = izs(2)
 !eignvec  = RESHAPE(rzs(1:KDIM*neignvec),(/KDIM,neignvec/))
 !eignval  = rzs(KDIM*neignvec+1:(KDIM+1)*neignvec)

 AMAT     = RESHAPE(DZS(1:NCHAN*KDIM),(/NCHAN,KDIM/))
 PX(KDIM) = 1.0 - SUM(PX(1:kdim-1))
! where (PX < 0.0) PX = 0.0
! where (PX > 1.0) PX = 1.0

 !ZX_EOF   = MATMUL(eignvec,eignval*PX)
!!! ZX_EOF   = MATMUL(eignvec,MATMUL(TRANSPOSE(eignvec),PX))
 zx_eof = PX
 zx     = zx_eof

 ! Softmax (= multiple-logistic) variable transform
 !beta = 100.0
  
 !zx = exp(beta*zx_eof) / SUM(exp(beta*zx_eof))
 !do ilev = 1, kdim
 !   do jlev = 1, kdim
 !      zgx(ilev,jlev) = - beta * zx(ilev) * zx(jlev)
 !      if (ilev == jlev) zgx(ilev,jlev) =  zgx(ilev,jlev) + zx(ilev) * beta
 !   end do
 !end do
 
 DO JCH=1,NCHAN 
   PF = PF + 0.5 * (SUM(ZX*AMAT(JCH,:)) - 1.0)**2
   DO JLEV=1,KDIM
      PG(JLEV) = PG(JLEV) + (AMAT(JCH,JLEV)-AMAT(JCH,KDIM)) * (SUM(ZX*AMAT(JCH,:)) - 1.0)
   ENDDO
 ENDDO
 !PG = MATMUL(PG, zgx)

 alpha = float(nchan)*100.0
 PF = PF + 0.5*alpha*SUM(ZX**2, MASK=ZX<0.0)
 WHERE (ZX<0.0) PG = PG + alpha*ZX

!write(*,'(a,2f10.2,50f6.1)') 'ACD_PX',PF,sqrt(sum(pg**2)),sum(px(1:kdim-1))*100.,PX*100.
!write(*,'(a,2f10.5,f10.2,50f7.2)') '888888 ',PF,sqrt(sum(pg**2)),sum(zx(1:kdim-1))*100.,ZX*100.

end subroutine da_cloud_sim
subroutine da_cloud_detect_airs(isensor,nchannels,ndim,kts,kte,n,iv)

!** *CLOUD_DETECT_AIRS* - CLOUD FLAGGING FOR AIRS AND IASI

!    AUTHOR: THOMAS AULIGNE      DATE : 01/08/2005
!
!    PURPOSE.
!    -------
!    FLAG THE PRESENCE OF CLOUD CONTAMINATION IN AIRS AND IASI CHANNELS
!
!**  INTERFACE.
!    ---------
!    WHERE nchannels    : Number of channels
!          kts          : model level corresponding to 100hPa (top of initial cloud search)
!          kte          : model level corresponding to surface (lower extent of cloud)            
!          rad_obs      : Potentially cloudy observations
!          rad_clr      : Clear radiance from Model
!          rad_ovc      : Model overcast radiance estimates
!          cloud_flag   : Cloud flag by channel; 1=clear, -1=cloudy
!
!**  EXTERNALS
!    ---------
!    N2QN1  - Minimization algorithm (double-precision constrained version of M1QN3)
!
!    MODIFICATIONS
!    -------------
!    NONE
!**  -----------------------------------------

IMPLICIT NONE

!* 0.1 Global arrays
INTEGER,INTENT(IN)    :: isensor             ! sensor index. 
INTEGER,INTENT(IN)    :: nchannels           ! number of channels 
INTEGER,INTENT(IN)    :: ndim                ! model levels between surface (lower extent of cloud) and 100hPa (top of cloud search)     
INTEGER,INTENT(IN)    :: kts                 ! model level corresponding to 100hPa (top of initial cloud search)
INTEGER,INTENT(IN)    :: kte                 ! model level corresponding to surface (lower extent of cloud)
INTEGER,INTENT(IN)    :: n                   ! pixel index 
type (iv_type), intent(inout)  :: iv         ! O-B structure.

INTEGER,PARAMETER     :: NITER             = 100
INTEGER,PARAMETER     :: NBAND             = 1
LOGICAL,PARAMETER     :: LPRECON           = .false.
INTEGER,PARAMETER     :: NEIGNVEC          = 4
INTEGER,PARAMETER     :: AIRS_Max_Channels = 2378

!! local declarations 
INTEGER               :: ichan(nchannels)    ! AIRS and IASI channel IDs
REAL                  :: rad_obs(nchannels)  ! Observed radiance
REAL                  :: rad_clr(nchannels)  ! Model clear radiance estimates
REAL                  :: rad_ovc(nchannels,ndim-1) ! RT overcast radiance estimates
double precision      :: px(ndim) !neignvec)        ! Cloud fractions
REAL                  :: rad_cld(nchannels)
INTEGER               :: ich,ilev,jlev,i,j,JBAND
double precision      :: ZF, ZF_CLR
double precision      :: ZG(ndim)
double precision      :: binf(ndim), bsup(ndim)
REAL                  :: AMAT(nchannels,ndim)
INTEGER               :: NCHAN
LOGICAL               :: LMATCH
INTEGER               :: Band_Size(5)
INTEGER               :: Bands(AIRS_Max_Channels,5) 
integer               :: cldtoplevel

! Hessian evaluation
 REAL                 :: hessian(ndim,ndim), eignvec(ndim,ndim), eignval(ndim)

!! local declarations for N2QN1 !!
INTEGER               :: NRZ, impres, io, IMODE, NSIM, nit, izs(2)
double precision      :: ZDF1, ZDXMIN, ZEPSG
double precision ,ALLOCATABLE :: ZRZ(:)
real, allocatable     :: RZS(:)
INTEGER, ALLOCATABLE  :: IZ(:)
DOUBLE PRECISION, ALLOCATABLE :: DZS(:)

REAL :: ZHOOK_HANDLE

! Initializations
      Band_Size(:)   = 0
      Bands(:,:)     = 0 
      Band_Size(1:5) = (/86, 0, 0, 16, 0 /)
 
      Bands(1:Band_Size(1),1) = &
&    (/                                                 &              !&      1,   6,   7,  10,  11,  15,  16,  17,  20,  21, &
&                                                       &              !&     22,  24,  27,  28,  30,  36,  39,  40,  42,  51, &
&                                                       &              !&     52,  54,  55,  56,  59,  62,  63,  68,  69,  71, &
&                                                       &              !&     72,  73,  74,  75,  76,  77,  78,  79,  80,  82, &
&                     92,  93,  98,  99, 101, 104, 105, &              !&     83,  84,  86,  92,  93,  98,  99, 101, 104, 105, &
&     108, 110, 111, 113, 116, 117, 123, 124, 128, 129, &
&     138, 139, 144, 145, 150, 151, 156, 157, 159, 162, &
&     165, 168, 169, 170, 172, 173, 174, 175, 177, 179, &
&     180, 182, 185, 186, 190, 192,      198, 201, 204, &              !&     180, 182, 185, 186, 190, 192, 193, 198, 201, 204, &
&     207, 210,      215, 216,      221,      226, 227, &              !&     207, 210, 213, 215, 216, 218, 221, 224, 226, 227, &
&     232,                     252, 253, 256, 257, 261, &              !&     232, 239, 248, 250, 251, 252, 253, 256, 257, 261, &
&     262, 267, 272, 295, 299,      305,           310, &              !&     262, 267, 272, 295, 299, 300, 305, 308, 309, 310, &
&          321, 325, 333, 338, 355, 362, 375, 453, 475, &              !&     318, 321, 325, 333, 338, 355, 362, 375, 453, 475, &
&     484, 497, 528, 587, 672, 787, 791, 843, 870, 914, &
&     950 /)

!      Bands(1:Band_Size(2),2) = &
!&    (/ 1003, 1012, 1019, 1024, 1030, 1038, 1048, 1069, 1079, 1082,  &
!&       1083, 1088, 1090, 1092, 1095, 1104, 1111, 1115, 1116, 1119,  &
!&       1120, 1123, 1130, 1138, 1142, 1178, 1199, 1206, 1221, 1237,  &
!&       1252, 1260, 1263, 1266, 1278, 1285 /)

!      Bands(1:Band_Size(3),3) = &
!&    (/       1301, 1304, 1329, 1371, 1382, 1415, 1424, 1449, 1455, &  !&    1290, 1301, 1304, 1329, 1371, 1382, 1415, 1424, 1449, 1455, &  
!&       1466,       1477,             1500, 1519,       1538, 1545, &  !&    1466, 1471, 1477, 1479, 1488, 1500, 1519, 1520, 1538, 1545, &  
!&       1565, 1574, 1583, 1593,       1627, 1636,       1652, 1669, &  !&    1565, 1574, 1583, 1593, 1614, 1627, 1636, 1644, 1652, 1669, & 
!&                   1694, 1708,       1723, 1740, 1748,       1756, &  !&    1674, 1681, 1694, 1708, 1717, 1723, 1740, 1748, 1751, 1756, &
!&             1766, 1771, 1777,       1783, 1794, 1800,       1806, &  !&    1763, 1766, 1771, 1777, 1780, 1783, 1794, 1800, 1803, 1806, &
!&             1826, 1843  /)                                           !&    1812, 1826, 1843  /)

      Bands(1:Band_Size(4),4) = &
&    (/ 1852, 1865, 1866,       1868, 1869, 1872, 1873,       1876, &  !&    1852, 1865, 1866, 1867, 1868, 1869, 1872, 1873, 1875, 1876, 
&             1881, 1882, 1883,                   1911, 1917, 1918, &  !&    1877, 1881, 1882, 1883, 1884, 1897, 1901, 1911, 1917, 1918, &
&                   1924, 1928        /)                               !&    1921, 1923, 1924, 1928, 1937  /)   

!      Bands(1:Band_Size(5),5) = &
!&    (/ 1938, 1939, 1941, 1946, 1947, 1948, 1958, 1971, 1973, 1988, &
!&       1995, 2084, 2085, 2097, 2098, 2099, 2100, 2101, 2103, 2104, &
!&       2106, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, 2115, &
!&       2116, 2117, 2118, 2119, 2120, 2121, 2122, 2123, 2128, 2134, &
!&       2141, 2145, 2149, 2153, 2164, 2189, 2197, 2209, 2226, 2234, &
!&       2280, 2318, 2321, 2325, 2328, 2333, 2339, 2348, 2353, 2355, &
!&       2363, 2370, 2371, 2377  /)  

    ichan        = iv%instid(isensor)%ichan(1:nchannels)
    rad_clr      = iv%instid(isensor)%rad_xb(1:nchannels,n)              !iv%instid(isensor)%tb_xb(1:nchan,n)
    rad_obs      = iv%instid(isensor)%rad_obs(1:nchannels,n)             !iv%instid(isensor)%tb_inv(1:nchan,n) + rad_clr
    rad_ovc      = iv%instid(isensor)%rad_ovc(1:nchannels,kts+1:kte,n)

    nchan        = 0
    AMAT(:,:)    = 0.0
    px(1:ndim-1) = 0.0
    px(ndim)     = 1.0
    ZF_CLR       = 0.0
    nit          = niter

!do ich=1,nchannels  
!   CALL CRTM_Planck_Radiance(11,ichan(ich),tb_obs(ich),rad_obs(ich))
!   CALL CRTM_Planck_Radiance(11,ichan(ich),tb_clr(ich),rad_clr(ich))
!end do                

    !--------------------!
    !   Loop over band   ! 
    !--------------------!
    BAND_LOOP: DO JBAND = 1, NBAND
      DO i = 1, Band_Size(JBAND)
        LMATCH = .FALSE.
        DO ich=1,nchannels        
          IF (ichan(ich)/= Bands(i,JBAND)) CYCLE 
          IF ((rad_obs(ich)<=0.0).OR.(rad_obs(ich)>1000.0)) CYCLE
          IF ((rad_clr(ich)<=0.0).OR.(rad_clr(ich)>1000.0)) CYCLE
          IF (ANY(rad_ovc(ich,1:NDIM-1)<=0.0)) CYCLE
          IF (ANY(rad_ovc(ich,1:NDIM-1)>1000.0)) CYCLE

          LMATCH               = .TRUE.    !! Found match for channel
          nchan                = nchan +1
          AMAT(nchan,1:ndim-1) = rad_ovc(ich,1:NDIM-1) / rad_obs(ich)
	  AMAT(nchan,ndim)     = rad_clr(ich) / rad_obs(ich)
          ZF_CLR               = ZF_CLR + 0.5*(AMAT(nchan,ndim)-1.0)**2   
        ENDDO
        IF (.NOT. LMATCH) WRITE(*,*) &
           'CLOUD_DETECT_AIRS: No matching for channel:',i,Bands(i,JBAND) 
      ENDDO	
    ENDDO BAND_LOOP                      ! Loop over band
    
    !--------------------!
    ! Hessian evaluation !
    !--------------------!
    IF (LPRECON) THEN
!      write(76,*) '**** HESSIAN ****'
      hessian(:,:)= 0.0      
      DO ilev=1, NDIM
        DO jlev=ilev, NDIM
          DO J=1,NCHAN
              hessian(ilev,jlev) = hessian(ilev,jlev)  + &
                                  (AMAT(J,ilev)-AMAT(J,NDIM)) * &
				  (AMAT(J,jlev)-AMAT(J,NDIM))
          ENDDO
          hessian(jlev,ilev) = hessian(ilev,jlev)   
        ENDDO
!	write(76,*) hessian(ilev,1:NDIM)
      ENDDO  
!!!      call da_eof_decomposition(ndim, hessian, eignvec, eignval)
    ENDIF
       
     !-----------------!
     ! n2qn1 minimizer !
     !-----------------!
      impres = 2
      io     = 66
      NSIM   = NITER+5
      ZDXMIN = 1.e-6
      ZEPSG  = 1.e-3 !e-9
      IMODE  = 1
      NRZ    = NDIM*(NDIM+9)/2 ! N2QN1
      ALLOCATE(IZ(2*NDIM +1))
      ALLOCATE(ZRZ(NRZ))
      ALLOCATE(DZS(NCHAN*NDIM))
      allocate(rzs(ndim*neignvec))
      binf   = -1000.0
      bsup   = 1000.0
      izs(1) = nchan
      izs(2) = neignvec
      rzs    = 0.0
      ZRZ    = 0.0
      dzs(1:NCHAN*NDIM)=RESHAPE(AMAT(1:NCHAN,1:NDIM),(/NCHAN*NDIM/))
	  
      IF (LPRECON) THEN
        IMODE = 2
        i     = 0
        DO ilev=1, NDIM
          DO jlev=ilev, NDIM
            i = i + 1
            ZRZ(i) = hessian(jlev,ilev)
          ENDDO
        ENDDO    
      ENDIF
!      rzs(1:ndim*neignvec)                  = RESHAPE(eignvec(1:ndim,1:neignvec),(/ndim*neignvec/))
!      rzs(ndim*neignvec+1:(ndim+1)*neignvec)= eignval(1:neignvec)

      call da_cloud_sim(0,NDIM,px,ZF,ZG,izs,RZS,DZS)
      ZDF1      = 1.e-1*ZF 


      call da_error("da_cloud_detect_airs.inc",228, &
             (/"inria_n2qn1 is not implemented here, please contact the author of this subroutine."/))
!     call inria_n2qn1(da_cloud_sim,NDIM,px,ZF,ZG,(/(ZDXMIN,jlev=1,NDIM)/),ZDF1, &
!                ZEPSG,impres,io,IMODE,nit,NSIM,binf,bsup,IZ,ZRZ,izs,RZS,DZS)
  
      IF (ALLOCATED(IZ))  DEALLOCATE(IZ)       
      IF (ALLOCATED(ZRZ)) DEALLOCATE(ZRZ)       
      IF (ALLOCATED(DZS)) DEALLOCATE(DZS)       
      if (allocated(rzs)) deallocate(rzs)                 
         
      !-----------------!
      ! Cloudy radiance !
      !-----------------!
      DO ich=1,nchannels
        rad_cld(ich) = SUM(px(1:ndim-1) * rad_ovc(ich,1:ndim-1)) + px(ndim) * rad_clr(ich) 
	
	if (ABS(rad_cld(ich)-rad_clr(ich)) < 0.01*rad_clr(ich)) then
	   iv%instid(isensor)%cloud_flag(ich,n) = qc_good
	else
   	   iv%instid(isensor)%cloud_flag(ich,n) = qc_bad
	end if   
      ENDDO    
      
    ! Dump cloud top pressure
    do ilev = kte, kts+2, -1
      if (px(ilev-kts+1) > 0.01) cldtoplevel = ilev
    end do   
    
    if (rtm_option == rtm_option_rttov) then
    elseif (rtm_option == rtm_option_crtm) then
       iv%instid(isensor)%clwp(n) = iv%instid(isensor)%pm(cldtoplevel,n)
    end if  	    
    
 
end subroutine da_cloud_detect_airs
subroutine da_cloud_detect_iasi(isensor,nchannels,ndim,kts,kte,n,iv)

!** *CLOUD_DETECT_AIRS* - CLOUD FLAGGING FOR AIRS AND IASI

!    AUTHOR: THOMAS AULIGNE      DATE : 01/08/2005
!
!    PURPOSE.
!    -------
!    FLAG THE PRESENCE OF CLOUD CONTAMINATION IN AIRS AND IASI CHANNELS
!
!**  INTERFACE.
!    ---------
!    WHERE nchannels    : Number of channels
!          kts          : model level corresponding to 100hPa (top of initial cloud search)
!          kte          : model level corresponding to surface (lower extent of cloud)            
!          rad_obs      : Potentially cloudy observations
!          rad_clr      : Clear radiance from Model
!          rad_ovc      : Model overcast radiance estimates
!          cloud_flag   : Cloud flag by channel; 1=clear, -1=cloudy
!
!**  EXTERNALS
!    ---------
!    N2QN1  - Minimization algorithm (double-precision constrained version of M1QN3)
!
!    MODIFICATIONS
!    -------------
!    NONE
!**  -----------------------------------------

IMPLICIT NONE

!* 0.1 Global arrays

INTEGER,INTENT(IN)    :: isensor             ! sensor index. 
INTEGER,INTENT(IN)    :: nchannels           ! number of channels 
INTEGER,INTENT(IN)    :: ndim                ! model levels between surface (lower extent of cloud) and 100hPa (top of cloud search)     
INTEGER,INTENT(IN)    :: kts                 ! model level corresponding to 100hPa (top of initial cloud search)
INTEGER,INTENT(IN)    :: kte                 ! model level corresponding to surface (lower extent of cloud)
INTEGER,INTENT(IN)    :: n                   ! pixel index 
type (iv_type), intent(inout)  :: iv         ! O-B structure.

INTEGER,PARAMETER     :: NITER             = 100
INTEGER,PARAMETER     :: NBAND             = 1
LOGICAL,PARAMETER     :: LPRECON           = .false.
INTEGER,PARAMETER     :: NEIGNVEC          = 4
INTEGER,PARAMETER     :: AIRS_Max_Channels = 2378
INTEGER,PARAMETER     :: IASI_Max_Channels = 8079
!! local declarations 
INTEGER               :: ichan(nchannels)    ! AIRS and IASI channel IDs
REAL                  :: rad_obs(nchannels)  ! Observed radiance
REAL                  :: rad_clr(nchannels)  ! Model clear radiance estimates
REAL                  :: rad_ovc(nchannels,ndim-1) ! RT overcast radiance estimates
double precision      :: px(ndim) !neignvec)        ! Cloud fractions
REAL                  :: rad_cld(nchannels)
INTEGER               :: ich,ilev,jlev,i,j,JBAND
double precision      :: ZF, ZF_CLR
double precision      :: ZG(ndim)
double precision      :: binf(ndim), bsup(ndim)
REAL                  :: AMAT(nchannels,ndim)
INTEGER               :: NCHAN,k
LOGICAL               :: LMATCH
INTEGER               :: Band_Size(5)
INTEGER               :: Bands(IASI_Max_Channels,5)
integer               :: cldtoplevel

! Hessian evaluation
 REAL                 :: hessian(ndim,ndim), eignvec(ndim,ndim), eignval(ndim)

!! local declarations for N2QN1 !!
INTEGER               :: NRZ, impres, io, IMODE, NSIM, nit, izs(2)
double precision      :: ZDF1, ZDXMIN, ZEPSG
double precision ,ALLOCATABLE :: ZRZ(:)
real, allocatable     :: RZS(:)
INTEGER, ALLOCATABLE  :: IZ(:)
DOUBLE PRECISION, ALLOCATABLE :: DZS(:)
INTEGER               :: gn
REAL :: ZHOOK_HANDLE

! Initializations
      Band_Size(:)   = 0
      Bands(:,:)     = 0 
      Band_Size(1:5) = (/ 193, 15, 116, 4, 15 /)
 
      Bands(1:Band_Size(1),1) = &
&      (/    16,   38,   49,   51,   55,   57,   59,   61,   63,   66, &
&            70,   72,   74,   79,   81,   83,   85,   87,   89,   92, &
&            95,   97,   99,  101,  104,  106,  109,  111,  113,  116, &
&           119,  122,  125,  128,  131,  133,  135,  138,  141,  144, &
&           146,  148,  151,  154,  157,  159,  161,  163,  165,  167, &
&           170,  173,  176,  178,  179,  180,  183,  185,  187,  189, &
&           191,  193,  195,  197,  199,  201,  203,  205,  207,  210, &
&           212,  214,  217,  219,  222,  224,  226,  228,  230,  232, &
&           234,  236,  239,  241,  242,  243,  246,  249,  252,  254, &
&           256,  258,  260,  262,  265,  267,  269,  271,  272,  273, &
&           275,  278,  280,  282,  284,  286,  288,  290,  292,  294, &
&           296,  299,  301,  303,  306,  308,  310,  312,  314,  316, &
&           318,  320,  323,  325,  327,  329,  331,  333,  335,  337, &
&           339,  341,  343,  345,  347,  350,  352,  354,  356,  358, &
&           360,  362,  364,  366,  369,  371,  373,  375,  377,  379, &
&           381,  383,  386,  389,  398,  401,  404,  407,  410,  414, &
&           416,  426,  428,  432,  434,  439,  445,  457,  515,  546, &
&           552,  559,  566,  571,  573,  646,  662,  668,  756,  867, &
&           921, 1027, 1090, 1133, 1191, 1194, 1271, 1805, 1884, 1946, &
&          1991, 2094, 2239 /)


      Bands(1:Band_Size(2),2) = &
&      (/ 1479, 1509, 1513, 1521, 1536, 1574, 1579, 1585, 1587, 1626, &
&         1639, 1643, 1652, 1658, 1671  /)

      Bands(1:Band_Size(3),3) = &
&      (/ 2119, 2213, 2271, 2321, 2398, 2701, 2741, 2819, 2889, 2907, 2910, &
&         2919, 2939, 2944, 2948, 2951, 2958, 2977, 2985, 2988, 2991, &
&         2993, 3002, 3008, 3014, 3027, 3029, 3036, 3047, 3049, 3053, &
&         3058, 3064, 3069, 3087, 3093, 3098, 3105, 3107, 3110, 3127, &
&         3136, 3151, 3160, 3165, 3168, 3175, 3178, 3207, 3228, 3244, &
&         3248, 3252, 3256, 3263, 3281, 3303, 3309, 3312, 3322, 3375, &
&         3378, 3411, 3438, 3440, 3442, 3444, 3446, 3448, 3450, 3452, &
&         3454, 3458, 3467, 3476, 3484, 3491, 3497, 3499, 3504, 3506, &
&         3509, 3518, 3527, 3555, 3575, 3577, 3580, 3582, 3586, 3589, &
&         3599, 3653, 3658, 3661, 4032, 5368, 5371, 5379, 5381, 5383, &
&         5397, 5399, 5401, 5403, 5405, 5455, 5480, 5483, 5485, 5492, &
&         5502, 5507, 5509, 5517, 5558  /)                                           !&    1812, 1826, 1843  /)

      Bands(1:Band_Size(4),4) = &
&      (/   5988, 5992, 5994, 6003  /)                              !&    1921, 1923, 1924, 1928, 1937  /)   

      Bands(1:Band_Size(5),5) = &
&      (/  6982, 6985, 6987, 6989, 6991, 6993, 6995, 6997, 7267, 7269, &
&          7424, 7426, 7428, 7885, 8007 /)   

    ichan        = iv%instid(isensor)%ichan(1:nchannels)
    rad_clr      = iv%instid(isensor)%rad_xb(1:nchannels,n)              !iv%instid(isensor)%tb_xb(1:nchan,n)
    rad_obs      = iv%instid(isensor)%rad_obs(1:nchannels,n)             !iv%instid(isensor)%tb_inv(1:nchan,n) + rad_clr
    rad_ovc      = iv%instid(isensor)%rad_ovc(1:nchannels,kts+1:kte,n)

    nchan        = 0
    AMAT(:,:)    = 0.0
    px(1:ndim-1) = 0.0
    px(ndim)     = 1.0
    ZF_CLR       = 0.0
    nit          = niter

!do ich=1,nchannels  
!   CALL CRTM_Planck_Radiance(11,ichan(ich),tb_obs(ich),rad_obs(ich))
!   CALL CRTM_Planck_Radiance(11,ichan(ich),tb_clr(ich),rad_clr(ich))
!end do                

    !--------------------!
    !   Loop over band   ! 
    !--------------------!
    BAND_LOOP: DO JBAND = 1, NBAND
      DO i = 1, Band_Size(JBAND)
        LMATCH = .FALSE.
        DO ich=1,nchannels        
          IF (ichan(ich)/= Bands(i,JBAND)) CYCLE 
          IF ((rad_obs(ich)<=0.0).OR.(rad_obs(ich)>1000.0)) CYCLE
          IF ((rad_clr(ich)<=0.0).OR.(rad_clr(ich)>1000.0)) CYCLE
          IF (ANY(rad_ovc(ich,1:NDIM-1)<=0.0)) CYCLE
          IF (ANY(rad_ovc(ich,1:NDIM-1)>1000.0)) CYCLE

          LMATCH               = .TRUE.    !! Found match for channel
          nchan                = nchan +1
          AMAT(nchan,1:ndim-1) = rad_ovc(ich,1:NDIM-1) / rad_obs(ich)
	  AMAT(nchan,ndim)     = rad_clr(ich) / rad_obs(ich)
          ZF_CLR               = ZF_CLR + 0.5*(AMAT(nchan,ndim)-1.0)**2   
        ENDDO
        IF (.NOT. LMATCH) then
           if (print_detail_rad) then
              write(unit=message(1),fmt='(A,2I8)') 'CLOUD_DETECT_IASI: No match for channel:',i,Bands(i,JBAND)
              call da_message(message(1:1))
           endif
        ENDIF
      ENDDO
    ENDDO BAND_LOOP                      ! Loop over band
    
    !--------------------!
    ! Hessian evaluation !
    !--------------------!
    IF (LPRECON) THEN
      hessian(:,:)= 0.0      
      DO ilev=1, NDIM
        DO jlev=ilev, NDIM
          DO J=1,NCHAN
              hessian(ilev,jlev) = hessian(ilev,jlev)  + &
                                  (AMAT(J,ilev)-AMAT(J,NDIM)) * &
				  (AMAT(J,jlev)-AMAT(J,NDIM))
          ENDDO
          hessian(jlev,ilev) = hessian(ilev,jlev)   
        ENDDO
      ENDDO  
    ENDIF
       
     !-----------------!
     ! n2qn1 minimizer !
     !-----------------!
      impres = 2
      io     = 66
      NSIM   = NITER+5
      ZDXMIN = 1.e-6
      ZEPSG  = 1.e-3 !e-9
      IMODE  = 1
      NRZ    = NDIM*(NDIM+9)/2 ! N2QN1
      ALLOCATE(IZ(2*NDIM +1))
      ALLOCATE(ZRZ(NRZ))
      ALLOCATE(DZS(NCHAN*NDIM))
      allocate(rzs(ndim*neignvec))
      binf   = -1000.0
      bsup   = 1000.0
      izs(1) = nchan
      izs(2) = neignvec
      rzs    = 0.0
      ZRZ    = 0.0
      dzs(1:NCHAN*NDIM)=RESHAPE(AMAT(1:NCHAN,1:NDIM),(/NCHAN*NDIM/))
	  
      IF (LPRECON) THEN
        IMODE = 2
        i     = 0
        DO ilev=1, NDIM
          DO jlev=ilev, NDIM
            i = i + 1
            ZRZ(i) = hessian(jlev,ilev)
          ENDDO
        ENDDO    
      ENDIF
!      rzs(1:ndim*neignvec)                  = RESHAPE(eignvec(1:ndim,1:neignvec),(/ndim*neignvec/))
!      rzs(ndim*neignvec+1:(ndim+1)*neignvec)= eignval(1:neignvec)

      call da_cloud_sim(0,NDIM,px,ZF,ZG,izs,RZS,DZS)
      ZDF1      = 1.e-1*ZF 


!      call da_error("da_cloud_detect_iasi.inc",233, &
!             (/"inria_n2qn1 is not implemented here, please contact the author of this subroutine."/))
!    call inria_n2qn1(da_cloud_sim,NDIM,px,ZF,ZG,(/(ZDXMIN,jlev=1,NDIM)/),ZDF1, &
!                ZEPSG,impres,io,IMODE,nit,NSIM,binf,bsup,IZ,ZRZ,izs,RZS,DZS)
  
      IF (ALLOCATED(IZ))  DEALLOCATE(IZ)       
      IF (ALLOCATED(ZRZ)) DEALLOCATE(ZRZ)       
      IF (ALLOCATED(DZS)) DEALLOCATE(DZS)       
      if (allocated(rzs)) deallocate(rzs)                 
         
      !-----------------!
      ! Cloudy radiance !
      !-----------------!
      DO ich=1,nchannels
        rad_cld(ich) = SUM(px(1:ndim-1) * rad_ovc(ich,1:ndim-1)) + px(ndim) * rad_clr(ich) 
	
	if (ABS(rad_cld(ich)-rad_clr(ich)) < 0.01*rad_clr(ich)) then
	   iv%instid(isensor)%cloud_flag(ich,n) = qc_good
	else
   	   iv%instid(isensor)%cloud_flag(ich,n) = qc_bad
	end if   
      ENDDO 
  
    ! Dump cloud top pressure
    do ilev = kte, kts+2, -1
      if (px(ilev-kts+1) > 0.01) cldtoplevel = ilev
    end do   
    
    if (rtm_option == rtm_option_rttov) then
    elseif (rtm_option == rtm_option_crtm) then
       iv%instid(isensor)%clwp(n) = iv%instid(isensor)%pm(cldtoplevel,n)
    end if  	    
    
end subroutine da_cloud_detect_iasi
subroutine da_qc_airs (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for AQUA/EOS-2-AIRS data.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.


   ! local variables
   integer   :: n,k,isflg,ios,fgat_rad_unit
   integer :: scanpos
   ! logical   :: lmix
   ! real    :: satzen

   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),nrej_limb,     &
                nrej_landsurface,nrej_windowchshort,nrej_windowchlong,    &
                nrej_clw,nrej_sst,nrej_topo, num_proc_domain,sensor_id,nrej_eccloud

   real      :: SST_model, SST_airs, SST_pred, diffSST, diffSST2
   real      :: inv_grosscheck

   character(len=30)  :: filename

   ! AIRS Cloud Detection Variables
   integer              :: kts_100hPa(1), kte_surf, ndim
   integer              :: numrad_local(nchan), numrad_global(nchan)
   real                 :: tstore
   real                 :: bias_local(nchan), bias_global(nchan)
   integer              :: kmin, kmax
   integer, allocatable  :: k_cloud_flag(:) ! cloud flags   
   if (trace_use_dull) call da_trace_entry("da_qc_airs")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_landsurface = 0
   nrej_windowchshort= 0
   nrej_windowchlong= 0
   nrej_sst= 0
   nrej_clw        = 0
   nrej_topo       = 0
   nrej_eccloud       = 0   
   sensor_id = 11
!   nrej_mixsurface = 0

   nrej_limb       = 0
   num_proc_domain = 0
   numrad_local    = 0
   bias_local      = 0.0
    
      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

         if (iv%instid(i)%info%proc_domain(1,n)) &
               num_proc_domain = num_proc_domain + 1

         !  0.0  initialise QC by flags assuming good obs
         !---------------------------------------------
         iv%instid(i)%tb_qc(:,n) = qc_good

         !  a.  reject channels over land/sea-ice/snow and mixture 
         !------------------------------------------------------------
         isflg = iv%instid(i)%isflg(n) 
         if (isflg > 0) then
            ! Check surface emissivity Jacobian 
            !-----------------------------------
            if (rtm_option == rtm_option_crtm .and. use_crtm_kmatrix) then
               do k = 1, nchan
                  if ( abs(iv%instid(i)%emiss_jacobian(k,n)) > 0.1 ) &
                  iv%instid(i)%tb_qc(k,n)  =  qc_bad
	       end do  	      
            end if
        
            if (use_rttov_kmatrix) then
               do k = 1, nchan
                  if ( abs(iv%instid(i)%emiss_jacobian(k,n)) > 0.1 ) &
                  iv%instid(i)%tb_qc(k,n)  =  qc_bad
               end do
            end if


            ! reject all channels 
            !--------------------
            if (only_sea_rad) then
	       iv%instid(i)%tb_qc(:,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_landsurface = nrej_landsurface + 1
	    end if	     
         end if

         ! a (bis) Check T and Q Jacobians for sensitivity to model top  
         !-----------------------------------------------------------
         if (rtm_option == rtm_option_crtm .and. use_crtm_kmatrix) then
            do k = 1, nchan
              if ( abs(iv%instid(i)%t_jacobian(k,1,n)) > 0.1 * SUM( &
	           abs(iv%instid(i)%t_jacobian(k,1:41,n))) ) &
                 iv%instid(i)%tb_qc(k,n)  =  qc_bad
              if ( abs(iv%instid(i)%q_jacobian(k,1,n)) > 0.1 * SUM( &
	           abs(iv%instid(i)%q_jacobian(k,1:41,n))) ) &
                 iv%instid(i)%tb_qc(k,n)  =  qc_bad
     	    end do
	 end if     
         
         if (use_rttov_kmatrix) then
            do k = 1, nchan
              if ( abs(iv%instid(i)%t_jacobian(k,1,n)) > 0.1 * SUM( &
                   abs(iv%instid(i)%t_jacobian(k,1:41,n))) ) &
                 iv%instid(i)%tb_qc(k,n)  =  qc_bad
              if ( abs(iv%instid(i)%q_jacobian(k,1,n)) > 0.1 * SUM( &
                   abs(iv%instid(i)%q_jacobian(k,1:41,n))) ) &
                 iv%instid(i)%tb_qc(k,n)  =  qc_bad
            end do
         end if  

!         !  a.  reject all channels over mixture surface type
!         !------------------------------------------------------
!         isflg = iv%instid(i)%isflg(n)
!         lmix  = (isflg==4) .or. (isflg==5) .or. (isflg==6) .or. (isflg==7)
!         if (lmix) then
!            iv%instid(i)%tb_qc(:,n)  =  qc_bad
!            if (iv%instid(i)%info%proc_domain(1,n)) &
!               nrej_mixsurface = nrej_mixsurface + 1
!         end if

         !  b.  reject limb obs 
         !------------------------------------------------------
         scanpos = iv%instid(i)%scanpos(n)
         if (scanpos <= 3 .or. scanpos >= 88) then
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_limb = nrej_limb + 1
         end if

	 
         !  c. Check for model clouds
         !-----------------------------------------------------------
         if (iv%instid(i)%clwp(n) > 0.05) then
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_clw = nrej_clw + 1
         end if

         !  d. Crude check for clouds in obs (assuming obs are used over ocean only)
         !   Use long wave window channel #914 - 10.662 nm (965.43 cm^-1)
         !   should be warmer than freezing temperature of the sea  
         !-----------------------------------------------------------
         !
         if(ob%instid(i)%tb(129,n) < 271.0) then
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_windowchlong = nrej_windowchlong + 1
         end if

         !  e. Check for contaminated obs in warmest near-infrared: Sun contamination during day 
         !-----------------------------------------------------------
         !
!         SST_airs=ob%instid(i)%tb(272,n)   !! short wave window channel #2333 - 3.882 nm (2616.38 cm^-1)
!         if(SST_airs > 307.0) then
!           iv%instid(i)%tb_qc(257:281,n)  = qc_bad
!            if (iv%instid(i)%info%proc_domain(1,n)) &
!               nrej_windowchshort = nrej_windowchshort + 1
!         end if


         !  f. Check for cloud free in obs (assuming obs are used over ocean only)
         !  Criterion: model SST within 2 K of transparent (hottest) short wave window channel
         !             includes check for contaminated near-infrared
         !-----------------------------------------------------------
         !
	 SST_model=iv%instid(i)%ts(n)       !! SST in the model
!         diffSST=abs(SST_model-SST_airs)
!         if(iv%instid(i)%solzen(n)>85.0 .and. diffSST > 2.0) then !! night-time only
!            iv%instid(i)%tb_qc(:,n)  = qc_bad
!            if (iv%instid(i)%info%proc_domain(1,n)) &
!               nrej_sst = nrej_sst + 1
!         end if
	    
	 !  g. Test on SST versus AIRS predicted SST from shortwave and longwave
	 !  Use channels #791, 914, 1285 and 1301.
	 !----------------------------------------------------------
         !
         SST_pred=8.28206-0.97957*ob%instid(i)%tb(126,n)+0.60529*ob%instid(i)%tb(129,n) + &
                  1.7444*ob%instid(i)%tb(165,n)-0.40379*ob%instid(i)%tb(166,n)
         diffSST2=SST_model-SST_pred
	 
         if ((diffSST2<-0.6).or.(diffSST2>3.5)) then
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_sst = nrej_sst + 1
         end if
	    
	    
         !  h. Test AIRS/VISNIR cloud fraction (in %) 
         !  Criterion : 5% cloud coverage within AIRS pixel
         !----------------------------------------------------------
         !
         if ((iv%instid(i)%rain_flag(n)>5).and.(iv%instid(i)%rain_flag(n)<100)) then
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_sst = nrej_sst + 1
         end if
 
 	    
         !  i. check surface height/pressure
         !-----------------------------------------------------------
         ! sfchgt = ivrad%info(n)%elv
         ! if (sfchgt >=) then
         ! else 
         ! end if

         !if ((isflg .ne. 0) .and. (iv%instid(i)%ps(n) < 850.0)) then
         !   iv%instid(i)%tb_qc(5,n)  = qc_bad
         !   if (iv%instid(i)%info%proc_domain(1,n)) &
         !      nrej_topo = nrej_topo + 1
         !end if

         !  j. check innovation
         !-----------------------------------------------------------
         do k = 1, nchan

           !  j.1. check absolute value of innovation
           !------------------------------------------------
	    inv_grosscheck = 15.0
 	    if (use_satcv(2)) inv_grosscheck = 100.0
            if (abs(iv%instid(i)%tb_inv(k,n)) > inv_grosscheck) then 
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_abs(k) = nrej_omb_abs(k) + 1
            end if
	    if (use_satcv(2)) cycle

           !  j.2. check relative value of innovation
           !      and assign of the observation error (standard deviation)
           !------------------------------------------------------------------------
            if (use_error_factor_rad) then         ! if use error tuning factor
               iv%instid(i)%tb_error(k,n) = &
                   satinfo(i)%error(k)*satinfo(i)%error_factor(k)
            else
               iv%instid(i)%tb_error(k,n) = satinfo(i)%error(k)
            end if
	    
	    ! M-estimator using Huber function (with k=sigmaO)
!	    if (abs(iv%instid(i)%tb_inv(k,n)) > iv%instid(i)%tb_error(k,n)) &
!	       iv%instid(i)%tb_error(k,n) = sqrt( &
!	       iv%instid(i)%tb_error(k,n) * abs(iv%instid(i)%tb_inv(k,n)) )

            if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_std(k) = nrej_omb_std(k) + 1
            end if
	    
	    if ( (iv%instid(i)%tb_qc     (k,n) == qc_good) .and. &
	         (iv%instid(i)%cloud_flag(k,n) == qc_good) ) then
	       bias_local(k)   = bias_local(k) + ob%instid(i)%tb(k,n) - iv%instid(i)%tb_xb(k,n)
	       numrad_local(k) = numrad_local(k) + 1
	    end if   
	    
         end do ! chan
      end do ! end loop pixel      
      			  
    ! Do inter-processor communication to gather statistics.
      
      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

           ! 1. Cloud detection scheme (NESDIS or MMR)
           !---------------------------------------------
         if (use_clddet_mmr .and. (.not.use_satcv(2))) then
            iv%instid(i)%cloud_flag(:,n) = qc_good
	    
	    if (rtm_option == rtm_option_rttov) then
            elseif (rtm_option == rtm_option_crtm) then
  	       kte_surf   = kte
   	       kts_100hPa = MAXLOC(iv%instid(i)%pm(kts:kte,n), &
	                    MASK = iv%instid(i)%pm(kts:kte,n) < 100.0)
            end if
	    
	    ndim = kte_surf - kts_100hPa(1) + 1
	 
            call da_cloud_detect_airs(i,nchan,ndim,kts_100hPa(1),kte_surf,n,iv)	       
         end if	          
	 
	 do k = 1, nchan
	    if (iv%instid(i)%cloud_flag(k,n) == qc_bad) iv%instid(i)%tb_qc(k,n) = qc_bad
         end do
 	 
           !  2. Check iuse from information file (channel selection)
           !-----------------------------------------------------------
         do k = 1, nchan
            if (satinfo(i)%iuse(k) .eq. -1) &
               iv%instid(i)%tb_qc(k,n)  = qc_bad
         end do
	      
           ! 3. Final QC decision
           !---------------------------------------------
         do k = 1, nchan
            if (iv%instid(i)%tb_qc(k,n) == qc_bad) then  ! bad obs
               iv%instid(i)%tb_error(k,n) = 500.0
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej(k) = nrej(k) + 1
            else                                         ! good obs
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  ngood(k) = ngood(k) + 1
            end if
         end do ! chan
	 
      end do ! end loop pixel 
 
   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int  (num_proc_domain)
   call da_proc_sum_int  (nrej_landsurface )
   call da_proc_sum_int  (nrej_windowchlong)
   call da_proc_sum_int  (nrej_windowchshort)
   call da_proc_sum_int  (nrej_sst)
   call da_proc_sum_int  (nrej_clw  )
   call da_proc_sum_int  (nrej_eccloud )   
   call da_proc_sum_int  (nrej_topo )
   call da_proc_sum_int  (nrej_limb)
   call da_proc_sum_ints (nrej_omb_abs(:))
   call da_proc_sum_ints (nrej_omb_std(:))
   call da_proc_sum_ints (nrej(:))
   call da_proc_sum_ints (ngood(:))

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it, '_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if

      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_airs.inc",410,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') ' Quality Control Statistics for '//iv%instid(i)%rttovid_string
      write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_landsurface  = ', nrej_landsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchlong = ', nrej_windowchlong
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchshort = ', nrej_windowchshort
      write(fgat_rad_unit,'(a20,i7)') ' nrej_sst = ', nrej_sst
      !      write(fgat_rad_unit,'(a20,i7)') ' nrej_mixsurface  = ', nrej_mixsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw         = ', nrej_clw
      write(fgat_rad_unit,'(a20,i7)') ' nrej_eccloud        = ', nrej_eccloud		  
      write(fgat_rad_unit,'(a20,i7)') ' nrej_topo        = ', nrej_topo
      write(fgat_rad_unit,'(a20,i7)') ' nrej_limb        = ', nrej_limb
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if

   if (trace_use_dull) call da_trace_exit("da_qc_airs")

end subroutine da_qc_airs

subroutine da_qc_amsua (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for amsua data.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.


   ! local variables
   integer   :: n,scanpos,k,isflg,ios,fgat_rad_unit
   logical   :: lmix
   real      :: si
   ! real    :: satzen
   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),      &
                nrej_mixsurface,nrej_windowchanl, nrej_si,    &
                nrej_clw,nrej_topo, num_proc_domain,  &
                nrej_limb

   character(len=30)  :: filename

   if (trace_use) call da_trace_entry("da_qc_amsua")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_mixsurface = 0
   nrej_windowchanl= 0
   nrej_si         = 0
   nrej_clw        = 0
   nrej_topo       = 0
   nrej_limb       = 0
   num_proc_domain = 0


      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

         if (iv%instid(i)%info%proc_domain(1,n)) &
               num_proc_domain = num_proc_domain + 1

         !  0.0  initialise QC by flags assuming good obs
         !---------------------------------------------
         iv%instid(i)%tb_qc(:,n) = qc_good

         !  a.  reject all channels over mixture surface type
         !------------------------------------------------------
         isflg = iv%instid(i)%isflg(n)
         lmix  = (isflg==4) .or. (isflg==5) .or. (isflg==6) .or. (isflg==7)
         if (lmix) then
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_mixsurface = nrej_mixsurface + 1
         end if
         !  b.  reject channels 1~4 over land/sea-ice/snow
         !------------------------------------------------------
         if (isflg > 0) then 
            iv%instid(i)%tb_qc(1:4,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_windowchanl = nrej_windowchanl + 1
           ! reject whole pixel if not over sea for global case
            if (global) iv%instid(i)%tb_qc(:,n)  = qc_bad
            if (only_sea_rad) iv%instid(i)%tb_qc(:,n)  = qc_bad
         end if

         !  c.  reject channels 13,14(above top model 10mb),15 
         !------------------------------------------------------
         iv%instid(i)%tb_qc(13:15,n)  = qc_bad

         !    reject limb obs 
         !------------------------------------------------------
         scanpos = iv%instid(i)%scanpos(n)
         if (scanpos <= 3 .or. scanpos >= 28) then
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_limb = nrej_limb + 1
         end if

         ! satzen  = rad%satzen
         ! if (abs(satzen) > 45.0) iv%instid(i)%tb_qc(:,n)  =  qc_bad

         !  d. check precipitation 
         !-----------------------------------------------------------
         if (ob%instid(i)%tb(1,n) > 0.0 .and. &
             ob%instid(i)%tb(15,n) > 0.0) then
            si = ob%instid(i)%tb(1,n) - ob%instid(i)%tb(15,n)
            if (si >= 3.0) then
               iv%instid(i)%tb_qc(:,n) = qc_bad
               iv%instid(i)%cloud_flag(:,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            end if
         end if

         if (iv%instid(i)%clwp(n) >= 0.2) then
            iv%instid(i)%tb_qc(:,n) = qc_bad
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_clw = nrej_clw + 1
         end if

         !   3.1 Estimate Cloud Liquid Water (CLW) in mm over sea
         !       (Grody etal. 2001, JGR, Equation 5b,7c,7d,9)
         !---------------------------------------------------------
         ! if (isflg == 0) then
         !    coszen =  cos(iv%instid(i)%satzen(n))
         !    d0     =  8.24-(2.622-1.846*coszen)*coszen
         !    d1     =  0.754
         !    d2     =  -2.265
         !    ts     =  iv%instid(i)%ts(n)
         !    tb1    =  ob%instid(i)%tb(1,n)
         !    tb2    =  ob%instid(i)%tb(2,n)
         !    clw    =  coszen*(d0+d1*log(ts-tb1)+d2*log(ts-tb2))
         !    clw    =  clw - 0.03
         ! end if


         !  e. check surface height/pressure
         !-----------------------------------------------------------
         ! sfchgt = ivrad%info(n)%elv
         ! if (sfchgt >=) then
         ! else 
         ! end if

         if ((isflg .ne. 0) .and. (iv%instid(i)%ps(n) < 850.0)) then
            iv%instid(i)%tb_qc(5,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_topo = nrej_topo + 1
         end if

         !  g. check iuse
         !-----------------------------------------------------------
         do k = 1, nchan
            if (satinfo(i)%iuse(k) .eq. -1) &
               iv%instid(i)%tb_qc(k,n)  = qc_bad
         end do

         !  f. check innovation
         !-----------------------------------------------------------
         do k = 1, nchan

         ! absolute departure check
            if (abs(iv%instid(i)%tb_inv(k,n)) > 15.0) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_abs(k) = nrej_omb_abs(k) + 1
            end if

         ! relative departure check
            if (use_error_factor_rad) then
               iv%instid(i)%tb_error(k,n) = &
                   satinfo(i)%error_std(k)*satinfo(i)%error_factor(k)
            else
               iv%instid(i)%tb_error(k,n) = satinfo(i)%error_std(k)
            end if

            if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
                iv%instid(i)%tb_qc(k,n)  = qc_bad
                if (iv%instid(i)%info%proc_domain(1,n)) &
                   nrej_omb_std(k) = nrej_omb_std(k) + 1
            end if

         ! final QC decsion
            if (iv%instid(i)%tb_qc(k,n) == qc_bad) then
               iv%instid(i)%tb_error(k,n) = 500.0
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej(k) = nrej(k) + 1
            else
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  ngood(k) = ngood(k) + 1
            end if

         end do ! chan
      end do ! end loop pixel
 
   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (num_proc_domain)
   call da_proc_sum_int (nrej_mixsurface)
   call da_proc_sum_int (nrej_windowchanl)
   call da_proc_sum_int (nrej_si )
   call da_proc_sum_int (nrej_clw)
   call da_proc_sum_int (nrej_topo)
   call da_proc_sum_int (nrej_limb)
   call da_proc_sum_ints (nrej_omb_abs(:))
   call da_proc_sum_ints (nrej_omb_std(:))
   call da_proc_sum_ints (nrej(:))
   call da_proc_sum_ints (ngood(:))

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if

      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_amsua.inc",207,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') ' Quality Control Statistics for '//iv%instid(i)%rttovid_string
      write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_mixsurface  = ', nrej_mixsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchanl = ', nrej_windowchanl
      write(fgat_rad_unit,'(a20,i7)') ' nrej_si          = ', nrej_si
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw         = ', nrej_clw
      write(fgat_rad_unit,'(a20,i7)') ' nrej_topo        = ', nrej_topo
      write(fgat_rad_unit,'(a20,i7)') ' nrej_limb        = ', nrej_limb
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if
   if (trace_use) call da_trace_exit("da_qc_amsua")

end subroutine da_qc_amsua


subroutine da_qc_amsub (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for amsub data.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.

   ! local variables
   integer   :: n,scanpos,k,isflg,ios,fgat_rad_unit
   logical   :: lmix
   real      :: si
   ! real     :: satzen
   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),      &
                nrej_mixsurface,nrej_windowchanl, nrej_si,    &
                nrej_clw,nrej_topo, num_proc_domain,  &
                nrej_limb

   character(len=30)  :: filename

   if (trace_use) call da_trace_entry("da_qc_amsub")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_mixsurface = 0
   nrej_windowchanl= 0
   nrej_si         = 0
   nrej_clw        = 0
   nrej_topo       = 0
   nrej_limb       = 0

   num_proc_domain = 0

      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

         if (iv%instid(i)%info%proc_domain(1,n)) &
               num_proc_domain = num_proc_domain + 1

         ! 0.0  initialise QC flags by assuming good obs
         !---------------------------------------------
         iv%instid(i)%tb_qc(:,n) = qc_good
         if (crtm_cloud) go to 2508

         ! a.  reject all channels over mixture surface type
         !------------------------------------------------------
         isflg = iv%instid(i)%isflg(n)
         lmix  = (isflg==4) .or. (isflg==5) .or. (isflg==6) .or. (isflg==7)
         if (lmix) then 
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_mixsurface = nrej_mixsurface + 1
         end if

         !  b.  reject channels 1~2 over land/sea-ice/snow
         !------------------------------------------------------
         if (isflg > 0) then
            iv%instid(i)%tb_qc(1:2,n)    = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_windowchanl = nrej_windowchanl + 1
            if (only_sea_rad) iv%instid(i)%tb_qc(:,n)    = qc_bad
         end if

         ! reject limb obs 
         !------------------------------------------------------
         scanpos = iv%instid(i)%scanpos(n)
         if (scanpos <= 8 .or. scanpos >= 83) then
             iv%instid(i)%tb_qc(:,n)  =  qc_bad
             if (iv%instid(i)%info%proc_domain(1,n)) &
                 nrej_limb = nrej_limb + 1
         end if

         ! satzen  = rad%satzen
         ! if (abs(satzen) > 45.0) &
         !    iv%instid(i)%tb_qc(:,n)  =  qc_bad

         !  d. check cloud/precipitation 
         !-----------------------------------------------------------
         if (ob%instid(i)%tb(1,n) > 0.0 .and. &
              ob%instid(i)%tb(2,n) > 0.0) then
            si = ob%instid(i)%tb(1,n) - ob%instid(i)%tb(2,n)
            if (si >= 3.0) then
               iv%instid(i)%tb_qc(:,n) = qc_bad
               iv%instid(i)%cloud_flag(:,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            end if
         end if

         if (iv%instid(i)%clwp(n) >= 0.2) then
            iv%instid(i)%tb_qc(:,n) = qc_bad
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_clw = nrej_clw + 1
         end if

         !  e. check surface height/pressure
         !------------------------------------------------------
         ! sfchgt = ivrad%info(n)%elv
         ! if (sfchgt >=) then
         !
         ! else 
         ! end if

         if ((isflg .ne. 0) .and. (iv%instid(i)%ps(n) < 800.0)) then
            iv%instid(i)%tb_qc(5,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_topo = nrej_topo + 1
         end if

         !  g. check iuse (pre-rejected channels by .info files)
         !------------------------------------------------------
         do k = 1, nchan
           if (satinfo(i)%iuse(k) .eq. -1) &
                iv%instid(i)%tb_qc(k,n) = qc_bad
         end do
 
2508      continue

         !  f. check innovation
         !-----------------------------------------------------------
         do k = 1, nchan
          ! absolute departure check
          if (.not. crtm_cloud) then
            if (abs(iv%instid(i)%tb_inv(k,n)) > 15.0) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_abs(k) = nrej_omb_abs(k) + 1
            end if
          end if 

          ! relative departure check
            if (use_error_factor_rad) then
                 iv%instid(i)%tb_error(k,n) =  &
                    satinfo(i)%error_std(k)*satinfo(i)%error_factor(k)
            else
                 iv%instid(i)%tb_error(k,n) = satinfo(i)%error_std(k)
            end if

          if (.not. crtm_cloud) then
            if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                     nrej_omb_std(k) = nrej_omb_std(k) + 1
            end if

           ! final QC decision
            if (iv%instid(i)%tb_qc(k,n) == qc_bad) then
               iv%instid(i)%tb_error(k,n) = 500.0
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej(k) = nrej(k) + 1
            else
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  ngood(k) = ngood(k) + 1
            end if
          end if

         end do ! chan
      end do ! end loop pixel

   ! Do inter-processor communication to gather statistics.

   call da_proc_sum_int (num_proc_domain)
   call da_proc_sum_int (nrej_mixsurface)
   call da_proc_sum_int (nrej_windowchanl)
   call da_proc_sum_int (nrej_si)
   call da_proc_sum_int (nrej_clw)
   call da_proc_sum_int (nrej_topo)
   call da_proc_sum_int (nrej_limb)
   call da_proc_sum_ints (nrej_omb_abs(:))
   call da_proc_sum_ints (nrej_omb_std(:))
   call da_proc_sum_ints (nrej(:))
   call da_proc_sum_ints (ngood(:))

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it, '_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it, '_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if
      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_amsub.inc",193,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') &
         'Quality Control Statistics for '//iv%instid(i)%rttovid_string
      write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_mixsurface  = ', nrej_mixsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchanl = ', nrej_windowchanl
      write(fgat_rad_unit,'(a20,i7)') ' nrej_si          = ', nrej_si
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw         = ', nrej_clw
      write(fgat_rad_unit,'(a20,i7)') ' nrej_topo        = ', nrej_topo
      write(fgat_rad_unit,'(a20,i7)') ' nrej_limb        = ', nrej_limb
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if

   if (trace_use) call da_trace_exit("da_qc_amsub")

end subroutine da_qc_amsub


subroutine da_qc_hirs (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for HIRS data.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.


   ! local variables
   integer   :: n,scanpos,k,isflg,ios,fgat_rad_unit
   logical   :: lmix
   ! real    :: satzen
   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),      &
                nrej_mixsurface,nrej_windowchanl, nrej_si,    &
                nrej_clw,nrej_topo, num_proc_domain,  &
                nrej_limb

   character(len=30)  :: filename

   if (trace_use) call da_trace_entry("da_qc_hirs")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_mixsurface = 0
   nrej_windowchanl= 0
   nrej_si         = 0
   nrej_clw        = 0
   nrej_topo       = 0
   nrej_limb       = 0
   num_proc_domain = 0

      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

         if (iv%instid(i)%info%proc_domain(1,n)) num_proc_domain = num_proc_domain + 1

         !  0.0  initialise QC by flags assuming good obs
         !---------------------------------------------
         iv%instid(i)%tb_qc(:,n) = qc_good

         !  a.  reject all channels over mixture surface type
         !------------------------------------------------------
         isflg = iv%instid(i)%isflg(n)
         lmix  = (isflg==4) .or. (isflg==5) .or. (isflg==6) .or. (isflg==7)
         if (lmix) then
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_mixsurface = nrej_mixsurface + 1
         end if

         !  b.  reject all channels over land/sea-ice/snow
         !------------------------------------------------------
         if (isflg > 0) then 
            iv%instid(i)%tb_qc(:,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_windowchanl = nrej_windowchanl + 1
            if (only_sea_rad) iv%instid(i)%tb_qc(:,n)    = qc_bad
         end if

         !  c.  reject channels 13,14(above top model 10mb),15 
         !------------------------------------------------------
         !iv%instid(i)%tb_qc(13:15,n)  = qc_bad

         !    reject limb obs 
         !------------------------------------------------------
         scanpos = iv%instid(i)%scanpos(n)
         if (scanpos <= 3 .or. scanpos >= 54) then
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_limb = nrej_limb + 1
         end if

         !  d. cloud detection to be added
         !-----------------------------------------------------------
         if (iv%instid(i)%clwp(n) >= 0.2) then
            iv%instid(i)%tb_qc(:,n) = qc_bad
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_clw = nrej_clw + 1
         end if

         !  e. check surface height/pressure
         !-----------------------------------------------------------
         ! sfchgt = ivrad%info(n)%elv
         ! if (sfchgt >=) then
         ! else 
         ! end if

         !if ((isflg .ne. 0) .and. (iv%instid(i)%ps(n) < 850.0)) then
         !   iv%instid(i)%tb_qc(5,n)  = qc_bad
         !   if (iv%instid(i)%info%proc_domain(1,n)) &
         !      nrej_topo = nrej_topo + 1
         !end if

         !  g. check iuse from information file (channel selection)
         !-----------------------------------------------------------
         do k = 1, nchan
            if (satinfo(i)%iuse(k) .eq. -1) &
               iv%instid(i)%tb_qc(k,n)  = qc_bad
         end do

         !  f. check innovation
         !-----------------------------------------------------------
         do k = 1, nchan

         !  1. check absolute value of innovation
         !------------------------------------------------
            if (abs(iv%instid(i)%tb_inv(k,n)) > 15.0) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_abs(k) = nrej_omb_abs(k) + 1
            end if

         !  2. check relative value of innovation
         !      and assign of the observation error (standard deviation)
         !------------------------------------------------------------------------
            if (use_error_factor_rad) then         ! if use error tuning factor
               iv%instid(i)%tb_error(k,n) = &
                  satinfo(i)%error(k)*satinfo(i)%error_factor(k)
            else
                iv%instid(i)%tb_error(k,n) = satinfo(i)%error(k)
            end if

            if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_std(k) = nrej_omb_std(k) + 1
            end if

           ! 3. Final QC decision
           !---------------------------------------------
            if (iv%instid(i)%tb_qc(k,n) == qc_bad) then  ! bad obs
               iv%instid(i)%tb_error(k,n) = 500.0
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej(k) = nrej(k) + 1
            else                                         ! good obs
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  ngood(k) = ngood(k) + 1
            end if
         end do ! chan
      end do ! end loop pixel
 
   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (num_proc_domain)
   call da_proc_sum_int (nrej_mixsurface)
   call da_proc_sum_int (nrej_windowchanl)
   call da_proc_sum_int (nrej_si )
   call da_proc_sum_int (nrej_clw)
   call da_proc_sum_int (nrej_topo)
   call da_proc_sum_int (nrej_limb)
   call da_proc_sum_ints (nrej_omb_abs(:))
   call da_proc_sum_ints (nrej_omb_std(:))
   call da_proc_sum_ints (nrej(:))
   call da_proc_sum_ints (ngood(:))

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it, '_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if

      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_hirs.inc",176,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') ' Quality Control Statistics for '//iv%instid(i)%rttovid_string
      write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_mixsurface  = ', nrej_mixsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchanl = ', nrej_windowchanl
      write(fgat_rad_unit,'(a20,i7)') ' nrej_si          = ', nrej_si
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw         = ', nrej_clw
      write(fgat_rad_unit,'(a20,i7)') ' nrej_topo        = ', nrej_topo
      write(fgat_rad_unit,'(a20,i7)') ' nrej_limb        = ', nrej_limb
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if

   if (trace_use) call da_trace_exit("da_qc_hirs")

end subroutine da_qc_hirs

subroutine da_qc_ssmis (it,i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for ssmis data.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.


   ! local variables
   integer   :: n,k,isflg,ios,fgat_rad_unit
   logical   :: lmix
   real      :: si37, q19, q37, term22v
   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),      &
                nrej_mixsurface,nrej_windowchanl, nrej_rain, nrej_si,    &
                nrej_clw,nrej_clw_rv,nrej_topo,num_proc_domain

   character(len=30)  :: filename

   if (trace_use) call da_trace_entry("da_qc_ssmis")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_mixsurface = 0
   nrej_windowchanl= 0
   nrej_rain       = 0
   nrej_si         = 0
   nrej_clw        = 0
   nrej_clw_rv     = 0
   nrej_topo       = 0
   num_proc_domain = 0

   do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

      if (iv%instid(i)%info%proc_domain(1,n)) &
            num_proc_domain = num_proc_domain + 1

      !  0.0  initialise QC by flags assuming good obs
      !---------------------------------------------
      iv%instid(i)%tb_qc(:,n) = qc_good

      !  a.  reject all channels over mixture surface type
      !------------------------------------------------------
      isflg = iv%instid(i)%isflg(n)
      lmix  = (isflg==4) .or. (isflg==5) .or. (isflg==6) .or. (isflg==7)
      if (lmix) then
         iv%instid(i)%tb_qc(:,n)  =  qc_bad
         if (iv%instid(i)%info%proc_domain(1,n)) &
            nrej_mixsurface = nrej_mixsurface + 1
      end if
      !  b.  reject channels 1~2,8 over land/sea-ice/snow
      !------------------------------------------------------
      if (isflg > 0) then 
         iv%instid(i)%tb_qc(1:2,n)  = qc_bad
         iv%instid(i)%tb_qc(8,n)    = qc_bad
         if (iv%instid(i)%info%proc_domain(1,n)) &
            nrej_windowchanl = nrej_windowchanl + 1
         if (only_sea_rad) iv%instid(i)%tb_qc(:,n)  = qc_bad
      end if

      !  c. reject rain_flagged data
      !------------------------------------------------------
      if (iv%instid(i)%rain_flag(n) == 1) then
         iv%instid(i)%tb_qc(:,n) = qc_bad
         iv%instid(i)%cloud_flag(:,n) = qc_bad
         if (iv%instid(i)%info%proc_domain(1,n)) &
            nrej_rain = nrej_rain + 1
      end if

      !  d. check precipitation 
      !  Ferraro, 1997: Journal of Geophysical Research Vol 102, 16715-16735
      !  SI37 = 62.18 + 0.773 * TB19v - TB37v
      !  Q19 = -2.70 * ( ln(290-TB19v) - 2.84 - 0.40 * ln(290-TB22v) )
      !  Q37 = -1.15 * ( ln(290-TB37v) - 2.99 - 0.32 * ln(290-TB22v) )
      !-----------------------------------------------------------

      if ( isflg >= 2 ) then
         si37 = 62.18 + 0.773 * ob%instid(i)%tb(13,n) - ob%instid(i)%tb(16,n)
         if ( si37 >= 5.0 ) then
            if ( ob%instid(i)%tb(14,n) >= 264.0 ) then  ! snow check
               if ( ob%instid(i)%tb(13,n)-ob%instid(i)%tb(12,n) <= 20.0 ) then  ! desert check
                  if ( ob%instid(i)%tb(16,n) <= 253.0 .and.   &
                       ob%instid(i)%tb(13,n)-ob%instid(i)%tb(12,n) <= 7.0 ) then  ! arid soil check
                     iv%instid(i)%tb_qc(:,n) = qc_bad
                     iv%instid(i)%cloud_flag(:,n) = qc_bad
                     if (iv%instid(i)%info%proc_domain(1,n)) &
                        nrej_si = nrej_si + 1
                   end if
               end if
            end if
          end if
      end if

      if ( isflg <= 1 ) then
         if ( ob%instid(i)%tb(14,n) < 290.0 ) then  ! assure positive for log
            term22v = log(290.0-ob%instid(i)%tb(14,n))
            if ( ob%instid(i)%tb(13,n) < 290.0 ) then  ! assure positive for log
               q19 = -2.70*(log(290.0-ob%instid(i)%tb(13,n))-2.84-0.40*term22v)
            end if
            if ( ob%instid(i)%tb(16,n) < 290.0 ) then  ! assure positive for log
               q37 = -1.15*(log(290.0-ob%instid(i)%tb(16,n))-2.99-0.32*term22v)
            end if
            if ( q19 >= 0.60 .or. q37 >= 0.20 ) then
               if ( ob%instid(i)%tb(14,n) >= 44.0+0.85*ob%instid(i)%tb(13,n) .or.  &
                    (ob%instid(i)%tb(14,n) <= 264.0 .and.    &
                     (ob%instid(i)%tb(14,n)-ob%instid(i)%tb(13,n)) >= 2.0) ) then
                  iv%instid(i)%tb_qc(:,n) = qc_bad
                  iv%instid(i)%cloud_flag(:,n) = qc_bad
                  if (iv%instid(i)%info%proc_domain(1,n)) &
                     nrej_clw_rv = nrej_clw_rv + 1      ! clw retrieval
               end if
            end if
         end if
      end if
       
      if (iv%instid(i)%clwp(n) >= 0.2) then
         iv%instid(i)%tb_qc(:,n) = qc_bad
         iv%instid(i)%cloud_flag(:,n) = qc_bad
         if (iv%instid(i)%info%proc_domain(1,n)) &
            nrej_clw = nrej_clw + 1               ! clw model
      end if

!      if ((isflg .ne. 0) .and. (iv%instid(i)%ps(n) < 850.0)) then
!         iv%instid(i)%tb_qc(5,n)  = qc_bad
!         if (iv%instid(i)%info%proc_domain(1,n)) &
!            nrej_topo = nrej_topo + 1
!      end if

      !  g. check iuse
      !-----------------------------------------------------------
      do k = 1, nchan
         if (satinfo(i)%iuse(k) .eq. -1) &
            iv%instid(i)%tb_qc(k,n)  = qc_bad
      end do

      !  f. check innovation
      !-----------------------------------------------------------
      do k = 1, nchan

         ! absolute departure check
         if (abs(iv%instid(i)%tb_inv(k,n)) > 15.0) then
            iv%instid(i)%tb_qc(k,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_omb_abs(k) = nrej_omb_abs(k) + 1
         end if

         ! relative departure check
         if (use_error_factor_rad) then
            iv%instid(i)%tb_error(k,n) = &
                satinfo(i)%error_std(k)*satinfo(i)%error_factor(k)
         else
            iv%instid(i)%tb_error(k,n) = satinfo(i)%error_std(k)
         end if

         if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
             iv%instid(i)%tb_qc(k,n)  = qc_bad
             if (iv%instid(i)%info%proc_domain(1,n)) &
                nrej_omb_std(k) = nrej_omb_std(k) + 1
         end if

         ! final QC decsion
         if (iv%instid(i)%tb_qc(k,n) == qc_bad) then
            iv%instid(i)%tb_error(k,n) = 500.0
            if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej(k) = nrej(k) + 1
            else
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  ngood(k) = ngood(k) + 1
            end if

      end do ! chan
   end do ! end loop pixel
 
   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (num_proc_domain)
   call da_proc_sum_int (nrej_mixsurface)
   call da_proc_sum_int (nrej_windowchanl)
   call da_proc_sum_int (nrej_rain )
   call da_proc_sum_int (nrej_si )
   call da_proc_sum_int (nrej_clw_rv)
   call da_proc_sum_int (nrej_clw)
!   call da_proc_sum_int (nrej_topo)
   call da_proc_sum_ints (nrej_omb_abs(:))
   call da_proc_sum_ints (nrej_omb_std(:))
   call da_proc_sum_ints (nrej(:))
   call da_proc_sum_ints (ngood(:))

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if

      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_ssmis.inc",208,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') ' Quality Control Statistics for '//iv%instid(i)%rttovid_string
      write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_mixsurface  = ', nrej_mixsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchanl = ', nrej_windowchanl
      write(fgat_rad_unit,'(a20,i7)') ' nrej_rain        = ', nrej_rain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_si          = ', nrej_si
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw_rv      = ', nrej_clw_rv
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw         = ', nrej_clw
!      write(fgat_rad_unit,'(a20,i7)') ' nrej_topo        = ', nrej_topo
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if

   if (trace_use) call da_trace_exit("da_qc_ssmis")

end subroutine da_qc_ssmis


subroutine da_qc_iasi (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for metop-a IASI data.
   !---------------------------------------------------------------------------

   implicit none
   
   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.


   ! local variables
   logical   :: lmix,lcould_read
   real    :: satzen
   integer   :: n,k,isflg,ios,fgat_rad_unit,sensor_id
   integer :: scanpos
   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),nrej_limb,     &
                nrej_landsurface,nrej_windowchshort,nrej_windowchlong,    &
                nrej_clw,nrej_sst,nrej_eccloud, num_proc_domain, nrej_mixsurface

   real      :: SST_model, SST_airs, SST_pred, diffSST, diffSST2
   real      :: inv_grosscheck

   character(len=30)  :: filename

   ! iasi Cloud Detection Variables
   integer              :: kmin, kmax, ndim, iunit
   integer, allocatable  :: k_cloud_flag(:) ! cloud flags   
   ! mmr Cloud Detection Variables
   integer              :: kts_100hPa(1), kte_surf,kts_200hPa(1)
   integer              :: numrad_local(nchan), numrad_global(nchan)
   real                 :: tstore
   real                 :: bias_local(nchan), bias_global(nchan)
   if (trace_use_dull) call da_trace_entry("da_qc_iasi")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_landsurface = 0
   nrej_windowchshort= 0
   nrej_windowchlong= 0
   nrej_sst= 0
   nrej_clw        = 0
   nrej_eccloud       = 0

!   nrej_mixsurface = 0

   nrej_limb       = 0
   num_proc_domain = 0
   numrad_local    = 0
   bias_local      = 0.0
   sensor_id = 16 


 
      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

         if (iv%instid(i)%info%proc_domain(1,n)) &
               num_proc_domain = num_proc_domain + 1

         !  0.0  initialise QC by flags assuming good obs
         !---------------------------------------------
         iv%instid(i)%tb_qc(:,n) = qc_good
		 
         !  a.  reject all channels over mixture surface type
         !------------------------------------------------------
         isflg = iv%instid(i)%isflg(n)
         lmix  = (isflg==4) .or. (isflg==5) .or. (isflg==6) .or. (isflg==7)
         if (lmix) then
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_mixsurface = nrej_mixsurface + 1
         end if	

        !  b.  reject limb obs 
         !------------------------------------------------------
         scanpos = iv%instid(i)%scanpos(n)
         if (scanpos <= 5 .or. scanpos >= 56) then
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_limb = nrej_limb + 1
         end if		 
         ! c.  reject channels from 565(Reject wavenumber > 2400 )
         !------------------------------------------------------
         iv%instid(i)%tb_qc(565:nchan,n)  = qc_bad
         ! c. cloud detection 
         !-----------------------------------------------------------
         if (iv%instid(i)%clwp(n) >= 0.2) then
            iv%instid(i)%tb_qc(:,n) = qc_bad
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_clw = nrej_clw + 1
         end if
	 
         !  d. check innovation
         !-----------------------------------------------------------
         do k = 1, nchan

           !  d.1. check absolute value of innovation
           !------------------------------------------------
	       inv_grosscheck = 15.0
 	       if (use_satcv(2)) inv_grosscheck = 100.0
              if (abs(iv%instid(i)%tb_inv(k,n)) > inv_grosscheck) then 
                 iv%instid(i)%tb_qc(k,n)  = qc_bad
                 if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_abs(k) = nrej_omb_abs(k) + 1
              end if


           !  d.2. check relative value of innovation
           !      and assign of the observation error (standard deviation)
           !------------------------------------------------------------------------
           if (use_error_factor_rad) then         ! if use error tuning factor
               iv%instid(i)%tb_error(k,n) = &
                   satinfo(i)%error(k)*satinfo(i)%error_factor(k)
           else
               iv%instid(i)%tb_error(k,n) = satinfo(i)%error(k)
           end if

           if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_std(k) = nrej_omb_std(k) + 1
           end if
	    
	       if ( (iv%instid(i)%tb_qc     (k,n) == qc_good) .and. &
	         (iv%instid(i)%cloud_flag(k,n) == qc_good) ) then
	        bias_local(k)   = bias_local(k) + ob%instid(i)%tb(k,n) - iv%instid(i)%tb_xb(k,n)
	        numrad_local(k) = numrad_local(k) + 1
	       end if 		    
         end do ! chan
      end do ! end loop pixel		 

      			  
    ! Do inter-processor communication to gather statistics.
      
      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

           ! 1. Cloud detection scheme (NESDIS or MMR)
           !---------------------------------------------
         if (use_clddet_mmr .and. (.not.use_satcv(2))) then
            iv%instid(i)%cloud_flag(:,n) = qc_good
	    
	    if (rtm_option == rtm_option_rttov) then
            else if (rtm_option == rtm_option_crtm) then
  	       kte_surf   = kte
   	       kts_100hPa = MAXLOC(iv%instid(i)%pm(kts:kte,n), &
	                    MASK = iv%instid(i)%pm(kts:kte,n) < 100.0)
            end if
	    
	    ndim = kte_surf - kts_100hPa(1) + 1
	 
            call da_cloud_detect_iasi(i,nchan,ndim,kts_100hPa(1),kte_surf,n,iv)	       
         end if	          
	 
	 do k = 1, nchan
	    if (iv%instid(i)%cloud_flag(k,n) == qc_bad) iv%instid(i)%tb_qc(k,n) = qc_bad
         end do
 	 
           !  2. Check iuse from information file (channel selection)
           !-----------------------------------------------------------
         do k = 1, nchan
            if (satinfo(i)%iuse(k) .eq. -1) &
               iv%instid(i)%tb_qc(k,n)  = qc_bad
         end do
	      
           ! 3. Final QC decision
           !---------------------------------------------
         do k = 1, nchan
            if (iv%instid(i)%tb_qc(k,n) == qc_bad) then  ! bad obs
               iv%instid(i)%tb_error(k,n) = 500.0
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej(k) = nrej(k) + 1
            else                                         ! good obs
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  ngood(k) = ngood(k) + 1
            end if
         end do ! chan
	 
      end do ! end loop pixel       

 
   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int  (num_proc_domain)
   call da_proc_sum_int  (nrej_landsurface )
   call da_proc_sum_int  (nrej_windowchlong)
   call da_proc_sum_int  (nrej_windowchshort)
   call da_proc_sum_int  (nrej_sst)
   call da_proc_sum_int  (nrej_clw  )
   call da_proc_sum_int  (nrej_eccloud )
   call da_proc_sum_int  (nrej_limb)
   call da_proc_sum_ints (nrej_omb_abs(:))
   call da_proc_sum_ints (nrej_omb_std(:))
   call da_proc_sum_ints (nrej(:))
   call da_proc_sum_ints (ngood(:))

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it, '_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if

      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_iasi.inc",282,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') ' Quality Control Statistics for '//iv%instid(i)%rttovid_string
      write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_landsurface  = ', nrej_landsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchlong = ', nrej_windowchlong
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchshort = ', nrej_windowchshort
      write(fgat_rad_unit,'(a20,i7)') ' nrej_sst = ', nrej_sst
      !      write(fgat_rad_unit,'(a20,i7)') ' nrej_mixsurface  = ', nrej_mixsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw         = ', nrej_clw
      write(fgat_rad_unit,'(a20,i7)') ' nrej_eccloud        = ', nrej_eccloud
      write(fgat_rad_unit,'(a20,i7)') ' nrej_limb        = ', nrej_limb
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if
    close(iunit)
   if (trace_use_dull) call da_trace_exit("da_qc_iasi")

end subroutine da_qc_iasi

subroutine da_qc_mhs (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for mhs data.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.

   ! local variables
   integer   :: n,scanpos,k,isflg,ios,fgat_rad_unit
   logical   :: lmix
   real      :: si
   ! real     :: satzen
   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),      &
                nrej_mixsurface,nrej_windowchanl, nrej_si,    &
                nrej_clw,nrej_topo, num_proc_domain,  &
                nrej_limb

   character(len=30)  :: filename

   if (trace_use) call da_trace_entry("da_qc_mhs")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_mixsurface = 0
   nrej_windowchanl= 0
   nrej_si         = 0
   nrej_clw        = 0
   nrej_topo       = 0
   nrej_limb       = 0

   num_proc_domain = 0

      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

         if (iv%instid(i)%info%proc_domain(1,n)) &
               num_proc_domain = num_proc_domain + 1

         ! 0.0  initialise QC flags by assuming good obs
         !---------------------------------------------
         iv%instid(i)%tb_qc(:,n) = qc_good

         ! a.  reject all channels over mixture surface type
         !------------------------------------------------------
         isflg = iv%instid(i)%isflg(n)
         lmix  = (isflg==4) .or. (isflg==5) .or. (isflg==6) .or. (isflg==7)
         if (lmix) then 
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_mixsurface = nrej_mixsurface + 1
         end if

         !  b.  reject channels 1~2 over land/sea-ice/snow
         !------------------------------------------------------
         if (isflg > 0) then
            iv%instid(i)%tb_qc(1:2,n)    = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_windowchanl = nrej_windowchanl + 1
            if (only_sea_rad) iv%instid(i)%tb_qc(:,n)    = qc_bad
         end if

         ! reject limb obs 
         !------------------------------------------------------
         scanpos = iv%instid(i)%scanpos(n)
         if (scanpos <= 8 .or. scanpos >= 83) then
             iv%instid(i)%tb_qc(:,n)  =  qc_bad
             if (iv%instid(i)%info%proc_domain(1,n)) &
                 nrej_limb = nrej_limb + 1
         end if

         ! satzen  = rad%satzen
         ! if (abs(satzen) > 45.0) &
         !    iv%instid(i)%tb_qc(:,n)  =  qc_bad

         !  d. check cloud/precipitation 
         !-----------------------------------------------------------
         if (ob%instid(i)%tb(1,n) > 0.0 .and. &
              ob%instid(i)%tb(2,n) > 0.0) then
            si = ob%instid(i)%tb(1,n) - ob%instid(i)%tb(2,n)
            if (si >= 3.0) then
               iv%instid(i)%tb_qc(:,n) = qc_bad
               iv%instid(i)%cloud_flag(:,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            end if
         end if

         if (iv%instid(i)%clwp(n) >= 0.2) then
            iv%instid(i)%tb_qc(:,n) = qc_bad
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_clw = nrej_clw + 1
         end if

         !  e. check surface height/pressure
         !------------------------------------------------------
         ! sfchgt = ivrad%info(n)%elv
         ! if (sfchgt >=) then
         !
         ! else 
         ! end if

         if ((isflg .ne. 0) .and. (iv%instid(i)%ps(n) < 800.0)) then
            iv%instid(i)%tb_qc(5,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_topo = nrej_topo + 1
         end if

         !  g. check iuse (pre-rejected channels by .info files)
         !------------------------------------------------------
         do k = 1, nchan
           if (satinfo(i)%iuse(k) .eq. -1) &
                iv%instid(i)%tb_qc(k,n) = qc_bad
         end do

         !  f. check innovation
         !-----------------------------------------------------------
         do k = 1, nchan
          ! absolute departure check
            if (abs(iv%instid(i)%tb_inv(k,n)) > 15.0) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_abs(k) = nrej_omb_abs(k) + 1
            end if

          ! relative departure check
            if (use_error_factor_rad) then
                 iv%instid(i)%tb_error(k,n) =  &
                    satinfo(i)%error_std(k)*satinfo(i)%error_factor(k)
            else
                 iv%instid(i)%tb_error(k,n) = satinfo(i)%error_std(k)
            end if

            if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                     nrej_omb_std(k) = nrej_omb_std(k) + 1
            end if

           ! final QC decision
            if (iv%instid(i)%tb_qc(k,n) == qc_bad) then
               iv%instid(i)%tb_error(k,n) = 500.0
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej(k) = nrej(k) + 1
            else
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  ngood(k) = ngood(k) + 1
            end if
         end do ! chan
      end do ! end loop pixel

   ! Do inter-processor communication to gather statistics.

   call da_proc_sum_int (num_proc_domain)
   call da_proc_sum_int (nrej_mixsurface)
   call da_proc_sum_int (nrej_windowchanl)
   call da_proc_sum_int (nrej_si)
   call da_proc_sum_int (nrej_clw)
   call da_proc_sum_int (nrej_topo)
   call da_proc_sum_int (nrej_limb)
   call da_proc_sum_ints (nrej_omb_abs(:))
   call da_proc_sum_ints (nrej_omb_std(:))
   call da_proc_sum_ints (nrej(:))
   call da_proc_sum_ints (ngood(:))

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it, '_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it, '_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if
      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_mhs.inc",185,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') &
         'Quality Control Statistics for '//iv%instid(i)%rttovid_string
      write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_mixsurface  = ', nrej_mixsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchanl = ', nrej_windowchanl
      write(fgat_rad_unit,'(a20,i7)') ' nrej_si          = ', nrej_si
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw         = ', nrej_clw
      write(fgat_rad_unit,'(a20,i7)') ' nrej_topo        = ', nrej_topo
      write(fgat_rad_unit,'(a20,i7)') ' nrej_limb        = ', nrej_limb
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if

   if (trace_use) call da_trace_exit("da_qc_mhs")

end subroutine da_qc_mhs


subroutine da_qc_mwts (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for fy3 mwts data.
   ! Dong Peiming modified from da_qc_amsua 2012/03/15
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.


   ! local variables
   integer   :: n,scanpos,k,isflg,ios,fgat_rad_unit
   logical   :: lmix
   real      :: si
   ! real    :: satzen
   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),      &
                nrej_mixsurface,nrej_windowchanl, nrej_si,    &
                nrej_clw,nrej_topo, num_proc_domain,  &
                nrej_limb

   character(len=30)  :: filename

   if (trace_use) call da_trace_entry("da_qc_mwts")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_mixsurface = 0
   nrej_windowchanl= 0
   nrej_si         = 0
   nrej_clw        = 0
   nrej_topo       = 0
   nrej_limb       = 0
   num_proc_domain = 0


      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

         if (iv%instid(i)%info%proc_domain(1,n)) &
               num_proc_domain = num_proc_domain + 1

         !  0.0  initialise QC by flags assuming good obs
         !---------------------------------------------
         iv%instid(i)%tb_qc(:,n) = qc_good

         !  a.  reject all channels over mixture surface type
         !------------------------------------------------------
         isflg = iv%instid(i)%isflg(n)
         lmix  = (isflg==4) .or. (isflg==5) .or. (isflg==6) .or. (isflg==7)
         if (lmix) then
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_mixsurface = nrej_mixsurface + 1
         end if
         !  b.  reject channels 1 over land/sea-ice/snow
         !------------------------------------------------------
         if (isflg > 0) then 
            iv%instid(i)%tb_qc(1:1,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_windowchanl = nrej_windowchanl + 1
           ! reject whole pixel if not over sea for global case
            if (global) iv%instid(i)%tb_qc(:,n)  = qc_bad
            if (only_sea_rad) iv%instid(i)%tb_qc(:,n)  = qc_bad
         end if

         !  c.  reject channels 13,14(above top model 10mb),15 
         !------------------------------------------------------
!Dongpm         iv%instid(i)%tb_qc(13:15,n)  = qc_bad

         !    reject limb obs 
         !------------------------------------------------------
         scanpos = iv%instid(i)%scanpos(n)
         if (scanpos <= 2 .or. scanpos >= 14) then
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_limb = nrej_limb + 1
         end if

         ! satzen  = rad%satzen
         ! if (abs(satzen) > 45.0) iv%instid(i)%tb_qc(:,n)  =  qc_bad

         !  d. check precipitation 
         !-----------------------------------------------------------
!Dongpm         if (ob%instid(i)%tb(1,n) > 0.0 .and. &
!Dongpm             ob%instid(i)%tb(15,n) > 0.0) then
!Dongpm            si = ob%instid(i)%tb(1,n) - ob%instid(i)%tb(15,n)
!Dongpm            if (si >= 3.0) then
!Dongpm               iv%instid(i)%tb_qc(:,n) = qc_bad
!Dongpm               iv%instid(i)%cloud_flag(:,n) = qc_bad
!Dongpm               if (iv%instid(i)%info%proc_domain(1,n)) &
!Dongpm                  nrej_si = nrej_si + 1
!Dongpm            end if
!Dongpm         end if
!Dongpm
         if (ob%instid(i)%tb(1,n) > 0.0) then
            si = iv%instid(i)%tb_inv(1,n)
            if (isflg .eq.0 .AND. si >= 3.0) then
               iv%instid(i)%tb_qc(1:3,n) = qc_bad
               iv%instid(i)%cloud_flag(1:3,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            elseif (isflg .gt.0 .AND. si >= 1.5) then
               iv%instid(i)%tb_qc(1:3,n) = qc_bad
               iv%instid(i)%cloud_flag(1:3,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            end if
         end if

         if (iv%instid(i)%clwp(n) >= 0.2) then
            iv%instid(i)%tb_qc(:,n) = qc_bad
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_clw = nrej_clw + 1
         end if

         !   3.1 Estimate Cloud Liquid Water (CLW) in mm over sea
         !       (Grody etal. 2001, JGR, Equation 5b,7c,7d,9)
         !---------------------------------------------------------
         ! if (isflg == 0) then
         !    coszen =  cos(iv%instid(i)%satzen(n))
         !    d0     =  8.24-(2.622-1.846*coszen)*coszen
         !    d1     =  0.754
         !    d2     =  -2.265
         !    ts     =  iv%instid(i)%ts(n)
         !    tb1    =  ob%instid(i)%tb(1,n)
         !    tb2    =  ob%instid(i)%tb(2,n)
         !    clw    =  coszen*(d0+d1*log(ts-tb1)+d2*log(ts-tb2))
         !    clw    =  clw - 0.03
         ! end if


         !  e. check surface height/pressure
         !-----------------------------------------------------------
         ! sfchgt = ivrad%info(n)%elv
         ! if (sfchgt >=) then
         ! else 
         ! end if

         if ((isflg .ne. 0) .and. (iv%instid(i)%ps(n) < 850.0)) then
            iv%instid(i)%tb_qc(2,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_topo = nrej_topo + 1
         end if

         !  g. check iuse
         !-----------------------------------------------------------
         do k = 1, nchan
            if (satinfo(i)%iuse(k) .eq. -1) &
               iv%instid(i)%tb_qc(k,n)  = qc_bad
         end do

         !  f. check innovation
         !-----------------------------------------------------------
         do k = 1, nchan

         ! absolute departure check
            if (abs(iv%instid(i)%tb_inv(k,n)) > 15.0) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_abs(k) = nrej_omb_abs(k) + 1
            end if

         ! relative departure check
            if (use_error_factor_rad) then
               iv%instid(i)%tb_error(k,n) = &
                   satinfo(i)%error_std(k)*satinfo(i)%error_factor(k)
            else
               iv%instid(i)%tb_error(k,n) = satinfo(i)%error_std(k)
            end if

            if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
                iv%instid(i)%tb_qc(k,n)  = qc_bad
                if (iv%instid(i)%info%proc_domain(1,n)) &
                   nrej_omb_std(k) = nrej_omb_std(k) + 1
            end if

         ! final QC decsion
            if (iv%instid(i)%tb_qc(k,n) == qc_bad) then
               iv%instid(i)%tb_error(k,n) = 500.0
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej(k) = nrej(k) + 1
            else
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  ngood(k) = ngood(k) + 1
            end if

         end do ! chan
      end do ! end loop pixel
 
   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (num_proc_domain)
   call da_proc_sum_int (nrej_mixsurface)
   call da_proc_sum_int (nrej_windowchanl)
   call da_proc_sum_int (nrej_si )
   call da_proc_sum_int (nrej_clw)
   call da_proc_sum_int (nrej_topo)
   call da_proc_sum_int (nrej_limb)
   call da_proc_sum_ints (nrej_omb_abs(:))
   call da_proc_sum_ints (nrej_omb_std(:))
   call da_proc_sum_ints (nrej(:))
   call da_proc_sum_ints (ngood(:))

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if

      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_mwts.inc",223,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') ' Quality Control Statistics for '//iv%instid(i)%rttovid_string
      write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_mixsurface  = ', nrej_mixsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchanl = ', nrej_windowchanl
      write(fgat_rad_unit,'(a20,i7)') ' nrej_si          = ', nrej_si
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw         = ', nrej_clw
      write(fgat_rad_unit,'(a20,i7)') ' nrej_topo        = ', nrej_topo
      write(fgat_rad_unit,'(a20,i7)') ' nrej_limb        = ', nrej_limb
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if
   if (trace_use) call da_trace_exit("da_qc_mwts")

end subroutine da_qc_mwts


subroutine da_qc_mwhs (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for fy3 mwhs data.
   ! Dong Peiming modified from da_qc_mhs 2012/3/15
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.

   ! local variables
   integer   :: n,scanpos,k,isflg,ios,fgat_rad_unit
   logical   :: lmix
   real      :: si
   ! real     :: satzen
   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),      &
                nrej_mixsurface,nrej_windowchanl, nrej_si,    &
                nrej_clw,nrej_topo, num_proc_domain,  &
                nrej_limb

   character(len=30)  :: filename

   if (trace_use) call da_trace_entry("da_qc_mwhs")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_mixsurface = 0
   nrej_windowchanl= 0
   nrej_si         = 0
   nrej_clw        = 0
   nrej_topo       = 0
   nrej_limb       = 0

   num_proc_domain = 0

      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

         if (iv%instid(i)%info%proc_domain(1,n)) &
               num_proc_domain = num_proc_domain + 1

         ! 0.0  initialise QC flags by assuming good obs
         !---------------------------------------------
         iv%instid(i)%tb_qc(:,n) = qc_good

         ! a.  reject all channels over mixture surface type
         !------------------------------------------------------
         isflg = iv%instid(i)%isflg(n)
         lmix  = (isflg==4) .or. (isflg==5) .or. (isflg==6) .or. (isflg==7)
         if (lmix) then 
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_mixsurface = nrej_mixsurface + 1
         end if

         !  b.  reject channels 1~2 over land/sea-ice/snow
         !------------------------------------------------------
         if (isflg > 0) then
            iv%instid(i)%tb_qc(1:2,n)    = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_windowchanl = nrej_windowchanl + 1
            if (only_sea_rad) iv%instid(i)%tb_qc(:,n)    = qc_bad
         end if

         ! reject limb obs 
         !------------------------------------------------------
         scanpos = iv%instid(i)%scanpos(n)
         if (scanpos <= 8 .or. scanpos >= 90) then
             iv%instid(i)%tb_qc(:,n)  =  qc_bad
             if (iv%instid(i)%info%proc_domain(1,n)) &
                 nrej_limb = nrej_limb + 1
         end if

         ! satzen  = rad%satzen
         ! if (abs(satzen) > 45.0) &
         !    iv%instid(i)%tb_qc(:,n)  =  qc_bad

         !  d. check cloud/precipitation 
         !-----------------------------------------------------------
         if (ob%instid(i)%tb(1,n) > 0.0 ) then
            si = abs(iv%instid(i)%tb_inv(1,n))
            if (isflg .eq.0 .AND. si >= 5.0) then
               iv%instid(i)%tb_qc(:,n) = qc_bad
               iv%instid(i)%cloud_flag(:,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            elseif (isflg .gt.0 .AND. si >= 3.0) then
               iv%instid(i)%tb_qc(:,n) = qc_bad
               iv%instid(i)%cloud_flag(:,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            end if
         end if

         if (iv%instid(i)%clwp(n) >= 0.2) then
            iv%instid(i)%tb_qc(:,n) = qc_bad
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_clw = nrej_clw + 1
         end if

         !  e. check surface height/pressure
         !------------------------------------------------------
         ! sfchgt = ivrad%info(n)%elv
         ! if (sfchgt >=) then
         !
         ! else 
         ! end if

         if ((isflg .ne. 0) .and. (iv%instid(i)%ps(n) < 800.0)) then
            iv%instid(i)%tb_qc(5,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_topo = nrej_topo + 1
         end if

         !  g. check iuse (pre-rejected channels by .info files)
         !------------------------------------------------------
         do k = 1, nchan
           if (satinfo(i)%iuse(k) .eq. -1) &
                iv%instid(i)%tb_qc(k,n) = qc_bad
         end do

         !  f. check innovation
         !-----------------------------------------------------------
         do k = 1, nchan
          ! absolute departure check
            if (abs(iv%instid(i)%tb_inv(k,n)) > 15.0) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_abs(k) = nrej_omb_abs(k) + 1
            end if

          ! relative departure check
            if (use_error_factor_rad) then
                 iv%instid(i)%tb_error(k,n) =  &
                    satinfo(i)%error_std(k)*satinfo(i)%error_factor(k)
            else
                 iv%instid(i)%tb_error(k,n) = satinfo(i)%error_std(k)
            end if

            if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                     nrej_omb_std(k) = nrej_omb_std(k) + 1
            end if

           ! final QC decision
            if (iv%instid(i)%tb_qc(k,n) == qc_bad) then
               iv%instid(i)%tb_error(k,n) = 500.0
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej(k) = nrej(k) + 1
            else
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  ngood(k) = ngood(k) + 1
            end if
         end do ! chan
      end do ! end loop pixel

   ! Do inter-processor communication to gather statistics.

   call da_proc_sum_int (num_proc_domain)
   call da_proc_sum_int (nrej_mixsurface)
   call da_proc_sum_int (nrej_windowchanl)
   call da_proc_sum_int (nrej_si)
   call da_proc_sum_int (nrej_clw)
   call da_proc_sum_int (nrej_topo)
   call da_proc_sum_int (nrej_limb)
   call da_proc_sum_ints (nrej_omb_abs(:))
   call da_proc_sum_ints (nrej_omb_std(:))
   call da_proc_sum_ints (nrej(:))
   call da_proc_sum_ints (ngood(:))

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it, '_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it, '_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if
      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_mwhs.inc",190,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') &
         'Quality Control Statistics for '//iv%instid(i)%rttovid_string
      write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_mixsurface  = ', nrej_mixsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchanl = ', nrej_windowchanl
      write(fgat_rad_unit,'(a20,i7)') ' nrej_si          = ', nrej_si
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw         = ', nrej_clw
      write(fgat_rad_unit,'(a20,i7)') ' nrej_topo        = ', nrej_topo
      write(fgat_rad_unit,'(a20,i7)') ' nrej_limb        = ', nrej_limb
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if

   if (trace_use) call da_trace_exit("da_qc_mwhs")

end subroutine da_qc_mwhs


subroutine da_qc_atms (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for atms data.
   ! Dongpm modified from atms 20120424
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.


   ! local variables
   integer   :: n,scanpos,k,isflg,ios,fgat_rad_unit
   logical   :: lmix
   real      :: si
   ! real    :: satzen
   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),      &
                nrej_mixsurface,nrej_windowchanl, nrej_si,    &
                nrej_clw,nrej_topo, num_proc_domain,  &
                nrej_limb

   character(len=30)  :: filename

   if (trace_use) call da_trace_entry("da_qc_atms")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_mixsurface = 0
   nrej_windowchanl= 0
   nrej_si         = 0
   nrej_clw        = 0
   nrej_topo       = 0
   nrej_limb       = 0
   num_proc_domain = 0


      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

         if (iv%instid(i)%info%proc_domain(1,n)) &
               num_proc_domain = num_proc_domain + 1

         !  0.0  initialise QC by flags assuming good obs
         !---------------------------------------------
         iv%instid(i)%tb_qc(:,n) = qc_good

         !  a.  reject all channels over mixture surface type
         !------------------------------------------------------
         isflg = iv%instid(i)%isflg(n)
         lmix  = (isflg==4) .or. (isflg==5) .or. (isflg==6) .or. (isflg==7)
         if (lmix) then
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_mixsurface = nrej_mixsurface + 1
         end if
         !  b.  reject channels 1~5 and 16~17 over land/sea-ice/snow
         !------------------------------------------------------
         if (isflg > 0) then 
            iv%instid(i)%tb_qc(1:5,n)  = qc_bad
            iv%instid(i)%tb_qc(16:17,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_windowchanl = nrej_windowchanl + 1
           ! reject whole pixel if not over sea for global case
            if (global) iv%instid(i)%tb_qc(:,n)  = qc_bad
            if (only_sea_rad) iv%instid(i)%tb_qc(:,n)  = qc_bad
         end if

         !  c.  reject channels 13,14(above top model 10mb)
         !------------------------------------------------------
         iv%instid(i)%tb_qc(13:14,n)  = qc_bad

         !    reject limb obs 
         !------------------------------------------------------
         scanpos = iv%instid(i)%scanpos(n)
         if (scanpos <= 0 .or. scanpos >= 97) then
            iv%instid(i)%tb_qc(:,n)  =  qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_limb = nrej_limb + 1
         end if

         ! satzen  = rad%satzen
         ! if (abs(satzen) > 45.0) iv%instid(i)%tb_qc(:,n)  =  qc_bad

         !  d. check precipitation 
         !-----------------------------------------------------------
         if (ob%instid(i)%tb(1,n) > 0.0) then
            si = iv%instid(i)%tb_inv(1,n)
            if (isflg .eq.0 .AND. si >= 3.0) then
               iv%instid(i)%tb_qc(1:8,n) = qc_bad
               iv%instid(i)%cloud_flag(1:8,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            elseif (isflg .gt.0 .AND. si >= 1.5) then
               iv%instid(i)%tb_qc(1:8,n) = qc_bad
               iv%instid(i)%cloud_flag(1:8,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            end if
         end if
!ECMWF for atms
         if (ob%instid(i)%tb(3,n) > 0.0) then
            si = abs(iv%instid(i)%tb_inv(3,n))
            if (si >= 5.0) then
               iv%instid(i)%tb_qc(1:8,n) = qc_bad
               iv%instid(i)%cloud_flag(1:8,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            endif
         endif

         if (ob%instid(i)%tb(16,n) > 0.0 .and. &
              ob%instid(i)%tb(17,n) > 0.0) then
            si = ob%instid(i)%tb(16,n) - ob%instid(i)%tb(17,n)
            if (si >= 3.0) then
               iv%instid(i)%tb_qc(16:22,n) = qc_bad
               iv%instid(i)%cloud_flag(16:22,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            end if
         end if
!ECMWF for atms
         if (ob%instid(i)%tb(3,n) > 0.0) then
            si = abs(iv%instid(i)%tb_inv(3,n))
            if (si >= 5.0) then
               iv%instid(i)%tb_qc(16:22,n) = qc_bad
               iv%instid(i)%cloud_flag(16:22,n) = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_si = nrej_si + 1
            endif
         endif


         if (iv%instid(i)%clwp(n) >= 0.2) then
            iv%instid(i)%tb_qc(:,n) = qc_bad
            iv%instid(i)%cloud_flag(:,n) = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_clw = nrej_clw + 1
         end if

         !   3.1 Estimate Cloud Liquid Water (CLW) in mm over sea
         !       (Grody etal. 2001, JGR, Equation 5b,7c,7d,9)
         !---------------------------------------------------------
         ! if (isflg == 0) then
         !    coszen =  cos(iv%instid(i)%satzen(n))
         !    d0     =  8.24-(2.622-1.846*coszen)*coszen
         !    d1     =  0.754
         !    d2     =  -2.265
         !    ts     =  iv%instid(i)%ts(n)
         !    tb1    =  ob%instid(i)%tb(1,n)
         !    tb2    =  ob%instid(i)%tb(2,n)
         !    clw    =  coszen*(d0+d1*log(ts-tb1)+d2*log(ts-tb2))
         !    clw    =  clw - 0.03
         ! end if


         !  e. check surface height/pressure
         !-----------------------------------------------------------
         ! sfchgt = ivrad%info(n)%elv
         ! if (sfchgt >=) then
         ! else 
         ! end if

         if ((isflg .ne. 0) .and. (iv%instid(i)%ps(n) < 850.0)) then
            iv%instid(i)%tb_qc(6,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_topo = nrej_topo + 1
         end if
         if ((isflg .ne. 0) .and. (iv%instid(i)%ps(n) < 800.0)) then
            iv%instid(i)%tb_qc(18,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
               nrej_topo = nrej_topo + 1
         end if

         !  g. check iuse
         !-----------------------------------------------------------
         do k = 1, nchan
            if (satinfo(i)%iuse(k) .eq. -1) &
               iv%instid(i)%tb_qc(k,n)  = qc_bad
         end do

         !  f. check innovation
         !-----------------------------------------------------------
         do k = 1, nchan

         ! absolute departure check
            if (abs(iv%instid(i)%tb_inv(k,n)) > 15.0) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_abs(k) = nrej_omb_abs(k) + 1
            end if

         ! relative departure check
            if (use_error_factor_rad) then
               iv%instid(i)%tb_error(k,n) = &
                   satinfo(i)%error_std(k)*satinfo(i)%error_factor(k)
            else
               iv%instid(i)%tb_error(k,n) = satinfo(i)%error_std(k)
            end if

            if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
                iv%instid(i)%tb_qc(k,n)  = qc_bad
                if (iv%instid(i)%info%proc_domain(1,n)) &
                   nrej_omb_std(k) = nrej_omb_std(k) + 1
            end if

         ! final QC decsion
            if (iv%instid(i)%tb_qc(k,n) == qc_bad) then
               iv%instid(i)%tb_error(k,n) = 500.0
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej(k) = nrej(k) + 1
            else
               if (iv%instid(i)%info%proc_domain(1,n)) &
                  ngood(k) = ngood(k) + 1
            end if

         end do ! chan
      end do ! end loop pixel
 
   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (num_proc_domain)
   call da_proc_sum_int (nrej_mixsurface)
   call da_proc_sum_int (nrej_windowchanl)
   call da_proc_sum_int (nrej_si )
   call da_proc_sum_int (nrej_clw)
   call da_proc_sum_int (nrej_topo)
   call da_proc_sum_int (nrej_limb)
   call da_proc_sum_ints (nrej_omb_abs(:))
   call da_proc_sum_ints (nrej_omb_std(:))
   call da_proc_sum_ints (nrej(:))
   call da_proc_sum_ints (ngood(:))

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if

      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_atms.inc",250,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') ' Quality Control Statistics for '//iv%instid(i)%rttovid_string
      write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_mixsurface  = ', nrej_mixsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchanl = ', nrej_windowchanl
      write(fgat_rad_unit,'(a20,i7)') ' nrej_si          = ', nrej_si
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw         = ', nrej_clw
      write(fgat_rad_unit,'(a20,i7)') ' nrej_topo        = ', nrej_topo
      write(fgat_rad_unit,'(a20,i7)') ' nrej_limb        = ', nrej_limb
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if
   if (trace_use) call da_trace_exit("da_qc_atms")

end subroutine da_qc_atms


subroutine da_qc_seviri (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for seviri data.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.


   ! local variables
   integer   :: n,scanpos,k,isflg,ios,fgat_rad_unit
   logical   :: lmix
   real      :: si
   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),      &
                nrej_mixsurface,nrej_windowchanl, nrej_si,    &
                nrej_clw,nrej_topo, num_proc_domain,  &
                nrej_limb

   character(len=30)  :: filename

   if (trace_use) call da_trace_entry("da_qc_seviri")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_mixsurface = 0
   nrej_windowchanl= 0
   nrej_si         = 0
   nrej_clw        = 0
   nrej_topo       = 0
   nrej_limb       = 0
   num_proc_domain = 0


      do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2

         if (iv%instid(i)%info%proc_domain(1,n)) &
               num_proc_domain = num_proc_domain + 1

         !  0.0  initialise QC by flags assuming good obs
         !---------------------------------------------
         iv%instid(i)%tb_qc(:,n) = qc_good

         ! observation errors
         ! assigned to satinfo%error_std instead of satinfo%error
         ! this is to make the function (read_biascoef) in da_radiance_init.inc work for seviri  
         do k = 1, nchan
            if (use_error_factor_rad) then
               iv%instid(i)%tb_error(k,n) = &
                   satinfo(i)%error_std(k)*satinfo(i)%error_factor(k)
            else
               iv%instid(i)%tb_error(k,n) = satinfo(i)%error_std(k)
            end if
        end do

         !  1.0 check iuse
         !-----------------------------------------------------------
         do k = 1, nchan
            if (satinfo(i)%iuse(k) .eq. -1) &
               iv%instid(i)%tb_qc(k,n)  = qc_bad
         end do

         !  2.0 check innovation
         !-----------------------------------------------------------
         do k = 1, nchan
            if ( iv%instid(i)%tb_qc(k,n) .eq. qc_good ) then

               ! absolute departure check
               if (abs(iv%instid(i)%tb_inv(k,n)) > 15.0) then
                  iv%instid(i)%tb_qc(k,n)  = qc_bad
                  if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_abs(k) = nrej_omb_abs(k) + 1
               end if

               ! relative departure check
               if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
                  iv%instid(i)%tb_qc(k,n)  = qc_bad
                  !iv%instid(i)%tb_error(k,n) = 500.0
                  if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_std(k) = nrej_omb_std(k) + 1
               end if

               ! final QC decsion
               if (iv%instid(i)%tb_qc(k,n) == qc_bad) then
                  iv%instid(i)%tb_error(k,n) = 500.0
                  if (iv%instid(i)%info%proc_domain(1,n)) &
                     nrej(k) = nrej(k) + 1
               else
                  if (iv%instid(i)%info%proc_domain(1,n)) &
                     ngood(k) = ngood(k) + 1
               end if

            end if   ! qc_good
         end do      ! nchan
      end do ! end loop pixel
 
   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (num_proc_domain)
   call da_proc_sum_int (nrej_mixsurface)
   call da_proc_sum_int (nrej_windowchanl)
   call da_proc_sum_int (nrej_si )
   call da_proc_sum_int (nrej_clw)
   call da_proc_sum_int (nrej_topo)
   call da_proc_sum_int (nrej_limb)
   call da_proc_sum_ints (nrej_omb_abs(:))
   call da_proc_sum_ints (nrej_omb_std(:))
   call da_proc_sum_ints (nrej(:))
   call da_proc_sum_ints (ngood(:))

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if

      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_seviri.inc",129,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') ' Quality Control Statistics for '//iv%instid(i)%rttovid_string
      write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      write(fgat_rad_unit,'(a20,i7)') ' nrej_mixsurface  = ', nrej_mixsurface
      write(fgat_rad_unit,'(a20,i7)') ' nrej_windowchanl = ', nrej_windowchanl
      write(fgat_rad_unit,'(a20,i7)') ' nrej_si          = ', nrej_si
      write(fgat_rad_unit,'(a20,i7)') ' nrej_clw         = ', nrej_clw
      write(fgat_rad_unit,'(a20,i7)') ' nrej_topo        = ', nrej_topo
      write(fgat_rad_unit,'(a20,i7)') ' nrej_limb        = ', nrej_limb
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if
   if (trace_use) call da_trace_exit("da_qc_seviri")

end subroutine da_qc_seviri


subroutine da_qc_amsr2 (it, i, nchan, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for amsr2  data.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)             :: it         ! outer loop count
   integer, intent(in)             :: i          ! sensor index.
   integer, intent(in)             :: nchan      ! number of channel
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.

   real, parameter  :: clwcutofx(1:14) =  &
          (/0.350, 0.350, 0.350, 0.350, 0.350, 0.350, 0.300, &
            0.300, 0.250, 0.250, 0.100, 0.100, 0.020, 0.020  /)

   ! local variables
   integer   :: n,k,isflg,ios,fgat_rad_unit
   real      :: si, coscon, sincon, bearaz, sun_zenith, sun_azimuth, sun_glint
   integer   :: ngood(nchan),nrej(nchan),nrej_omb_abs(nchan), &
                nrej_omb_std(nchan),      &
                nrej_clw(nchan), num_proc_domain,  &
                nrej_land, nrej_glint

   character(len=30)  :: filename

   if (trace_use) call da_trace_entry("da_qc_amsr2")

   ngood(:)        = 0
   nrej(:)         = 0
   nrej_omb_abs(:) = 0
   nrej_omb_std(:) = 0
   nrej_clw(:)     = 0
   nrej_land       = 0
   nrej_glint      = 0
   num_proc_domain = 0

   coscon = cos( (90.0 - 55.0)*deg_to_rad)
   sincon = sin( (90.0 - 55.0)*deg_to_rad)

   do n= iv%instid(i)%info%n1,iv%instid(i)%info%n2
      if (iv%instid(i)%info%proc_domain(1,n)) &
            num_proc_domain = num_proc_domain + 1

      !  0.0  initialise QC by flags assuming good obs
      !-----------------------------------------------------------------
      iv%instid(i)%tb_qc(:,n) = qc_good

      !  1.0 reject all channels over land
      !-----------------------------------------------------------------
      isflg = iv%instid(i)%isflg(n) !model surface type at ob location
      if (isflg > 0 .or. iv%instid(i)%landsea_mask(n) == 0) then
         iv%instid(i)%tb_qc(:,n)  = qc_bad
         nrej_land = nrej_land + 1
      end if

      !  2.0 check sun_glint angle for 6.9, 7.3 and 10.6 GHz
      !      It should be >=  25 degrees
      !-----------------------------------------------------------------
      bearaz      = 270.0 - iv%instid(i)%satazi(n)
      sun_zenith  = iv%instid(i)%solzen(n)
      sun_azimuth = 90.0 - iv%instid(i)%solazi(n)

      sun_glint = acos(coscon * cos( bearaz ) * cos( sun_zenith ) * cos( sun_azimuth ) + &
                  coscon * sin( bearaz ) * cos( sun_zenith ) * sin( sun_azimuth ) +  &
                  sincon * sin( sun_zenith )) * rad_to_deg

      if (sun_glint < 25.0) then
         ! apply only to 6.9,7.3 and 10.6 GHz for both V & H polarizations
         iv%instid(i)%tb_qc(1:6,n)  = qc_bad
         nrej_glint = nrej_glint + 1
      end if

      !  3.0 check iuse
      !-----------------------------------------------------------------
      do k = 1, nchan
         if (satinfo(i)%iuse(k) .eq. -1) &
               iv%instid(i)%tb_qc(k,n)  = qc_bad
      end do

      !  4.0 check cloud
      !-----------------------------------------------------------------
      if (.not. crtm_cloud ) then
         do k = 1, nchan
            ! clw check
            ! channel dependent criteria
            if( iv%instid(i)%clw(n) < 0 .or.             &  !bad obs retrieved clw
                iv%instid(i)%clw(n)  > clwcutofx(k) .or. &  !obs retrieved clw
                iv%instid(i)%clwp(n) > clwcutofx(k) ) then  !model/guess clwp
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                     nrej_clw(k) = nrej_clw(k) + 1
            end if
         end do
      end if

      !  5.0 check innovation
      !-----------------------------------------------------------------
      do k = 1, nchan

         ! absolute departure check
         if ( k <= 7 .or. k == 11 .or. k == 12) then
            if (abs(iv%instid(i)%tb_inv(k,n)) > 6.0) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                     nrej_omb_abs(k) = nrej_omb_abs(k) + 1
            end if
         else if ( k == 8 .or. k == 9 ) then
            if (abs(iv%instid(i)%tb_inv(k,n)) > 8.0) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                     nrej_omb_abs(k) = nrej_omb_abs(k) + 1
            end if
         else
            if (abs(iv%instid(i)%tb_inv(k,n)) > 10.0) then
               iv%instid(i)%tb_qc(k,n)  = qc_bad
               if (iv%instid(i)%info%proc_domain(1,n)) &
                     nrej_omb_abs(k) = nrej_omb_abs(k) + 1
            end if
         end if

         ! relative departure check
         if (use_error_factor_rad) then
            iv%instid(i)%tb_error(k,n) = &
                satinfo(i)%error_std(k)*satinfo(i)%error_factor(k)
         else
            iv%instid(i)%tb_error(k,n) = satinfo(i)%error_std(k)
         end if
         if (abs(iv%instid(i)%tb_inv(k,n)) > 3.0*iv%instid(i)%tb_error(k,n)) then
            iv%instid(i)%tb_qc(k,n)  = qc_bad
            if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej_omb_std(k) = nrej_omb_std(k) + 1
         end if

         ! final QC decsion
         if (iv%instid(i)%tb_qc(k,n) == qc_bad) then
            iv%instid(i)%tb_error(k,n) = 500.0
            if (iv%instid(i)%info%proc_domain(1,n)) &
                  nrej(k) = nrej(k) + 1
         else
            if (iv%instid(i)%info%proc_domain(1,n)) &
                  ngood(k) = ngood(k) + 1
         end if
      end do      ! nchan

   end do ! end loop pixel
 
   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (num_proc_domain)
   call da_proc_sum_ints (nrej_clw)
   call da_proc_sum_int (nrej_land)
   call da_proc_sum_int (nrej_glint)
   call da_proc_sum_ints (nrej_omb_abs)
   call da_proc_sum_ints (nrej_omb_std)
   call da_proc_sum_ints (nrej)
   call da_proc_sum_ints (ngood)

   if (rootproc) then
      if (num_fgat_time > 1) then
         write(filename,'(i2.2,a,i2.2)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)//'_',iv%time
      else
         write(filename,'(i2.2,a)') it,'_qcstat_'//trim(iv%instid(i)%rttovid_string)
      end if

      call da_get_unit(fgat_rad_unit)
      open(fgat_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         write(unit=message(1),fmt='(A,A)') 'error opening the output file ', filename
         call da_error("da_qc_amsr2.inc",171,message(1:1))
      end if

      write(fgat_rad_unit, fmt='(/a/)') ' Quality Control Statistics for '//iv%instid(i)%rttovid_string
      if(num_proc_domain > 0) write(fgat_rad_unit,'(a20,i7)') ' num_proc_domain  = ', num_proc_domain
      if(nrej_land > 0)write(fgat_rad_unit,'(a20,i7)') ' nrej_land        = ', nrej_land
      if(nrej_glint> 0)write(fgat_rad_unit,'(a20,i7)') ' nrej_glint       = ', nrej_glint
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_abs(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_abs(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_omb_std(:)  = '
      write(fgat_rad_unit,'(10i7)')     nrej_omb_std(:)
      write(fgat_rad_unit,'(a20)')    ' nrej_clw(:)      = '
      write(fgat_rad_unit,'(10i7)')     nrej_clw(:)
      write(fgat_rad_unit,'(a20)')    ' nrej(:)          = '
      write(fgat_rad_unit,'(10i7)')     nrej(:)
      write(fgat_rad_unit,'(a20)')    ' ngood(:)         = '
      write(fgat_rad_unit,'(10i7)')     ngood(:)

      close(fgat_rad_unit)
      call da_free_unit(fgat_rad_unit)
   end if
   if (trace_use) call da_trace_exit("da_qc_amsr2")

end subroutine da_qc_amsr2
subroutine da_write_iv_rad_ascii (it, ob, iv )

   !---------------------------------------------------------------------------
   ! Purpose: write out innovation vector structure for radiance data.
   !---------------------------------------------------------------------------

   implicit none

   integer      ,     intent(in)  :: it       ! outer loop count
   type (y_type),     intent(in)  :: ob       ! Observation structure.
   type (iv_type),    intent(in)  :: iv       ! O-B structure.

   integer                        :: n        ! Loop counter.
   integer                        :: i, k, l  ! Index dimension.
   integer                        :: nlevelss ! Number of obs levels.

   integer            :: ios, innov_rad_unit
   character(len=filename_len)  :: filename
   character(len=7)   :: surftype
   integer            :: ndomain
   logical            :: amsr2

   real, allocatable :: dtransmt(:,:), transmt_jac(:,:), transmt(:,:), lod(:,:), lod_jac(:,:)

   if (trace_use) call da_trace_entry("da_write_iv_rad_ascii")

   write(unit=message(1),fmt='(A)') 'Writing radiance OMB ascii file'
   call da_message(message(1:1))

   do i = 1, iv%num_inst
      if (iv%instid(i)%num_rad < 1) cycle

      ! count number of obs within the loc%proc_domain
      ! ---------------------------------------------
      ndomain = 0
      do n =1,iv%instid(i)%num_rad
         if (iv%instid(i)%info%proc_domain(1,n)) then
            ndomain = ndomain + 1
         end if
      end do
      if (ndomain < 1) cycle

      if (rtm_option==rtm_option_crtm .and. write_jacobian ) then
         allocate ( dtransmt(iv%instid(i)%nchan,iv%instid(i)%nlevels) )
         allocate ( transmt_jac(iv%instid(i)%nchan,iv%instid(i)%nlevels) )
         allocate ( transmt(iv%instid(i)%nchan,iv%instid(i)%nlevels) )
         allocate ( lod(iv%instid(i)%nchan,iv%instid(i)%nlevels) )
         allocate ( lod_jac(iv%instid(i)%nchan,iv%instid(i)%nlevels) )
      end if

      amsr2 = index(iv%instid(i)%rttovid_string,'amsr2') > 0

      write(unit=filename, fmt='(i2.2,a,i4.4)') it,'_inv_'//trim(iv%instid(i)%rttovid_string)//'.', myproc

      call da_get_unit(innov_rad_unit)
      open(unit=innov_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0 ) then
         call da_error("da_write_iv_rad_ascii.inc",58, &
            (/"Cannot open innovation radiance file"//filename/))
      Endif
      write(unit=innov_rad_unit,fmt='(a,a,i7,a,i5,a)') trim(iv%instid(i)%rttovid_string), &
                        ' number-of-pixels : ', ndomain, &
                        ' channel-number-of-each-pixel : ', iv%instid(i)%nchan, &
                        ' index-of-channels : '
      write(unit=innov_rad_unit,fmt='(10i5)') iv%instid(i)%ichan
      if ( amsr2 ) then
         write(unit=innov_rad_unit,fmt='(a)') ' pixel-info : i date scanpos landsea_mask  elv lat lon  satzen satazi clw'
      else
         write(unit=innov_rad_unit,fmt='(a)') ' pixel-info : i date scanpos landsea_mask  elv lat lon  satzen satazi'
      end if
      write(unit=innov_rad_unit,fmt='(a)') ' grid%xb-surf-info : i t2m mr2m(ppmv) u10 v10 ps ts smois tslb snowh isflg &
                    & soiltyp vegtyp vegfra elev clwp'
      ndomain = 0
      do n =1,iv%instid(i)%num_rad
         if (iv%instid(i)%info%proc_domain(1,n)) then
            ndomain=ndomain+1
            if ( amsr2 ) then ! write out clw
               write(unit=innov_rad_unit,fmt='(a,i7,2x,a,i6,i3,f6.0,4f8.2,f8.3)') 'INFO : ', ndomain, &
                                iv%instid(i)%info%date_char(n), &
                                iv%instid(i)%scanpos(n),        &
                                iv%instid(i)%landsea_mask(n),   &
                                iv%instid(i)%info%elv(n),       &
                                iv%instid(i)%info%lat(1,n),     &
                                iv%instid(i)%info%lon(1,n),     &
                                iv%instid(i)%satzen(n),         &
                                iv%instid(i)%satazi(n),         &
                                iv%instid(i)%clw(n)
            else ! no clw info
               write(unit=innov_rad_unit,fmt='(a,i7,2x,a,i6,i3,f6.0,4f8.2)') 'INFO : ', ndomain, &
                                iv%instid(i)%info%date_char(n), &
                                iv%instid(i)%scanpos(n),        &
                                iv%instid(i)%landsea_mask(n),   &
                                iv%instid(i)%info%elv(n),       &
                                iv%instid(i)%info%lat(1,n),     &
                                iv%instid(i)%info%lon(1,n),     &
                                iv%instid(i)%satzen(n),         &
                                iv%instid(i)%satazi(n)
            end if
            select case (iv%instid(i)%isflg(n))
            case (0) ;
               surftype = ' SEA : '
            case (1) ;
               surftype = ' ICE : '
            case (2) ;
               surftype = 'LAND : '
            case (3) ;
               surftype = 'SNOW : '
            case (4) ;
               surftype = 'MSEA : '
            case (5) ;
               surftype = 'MICE : '
            case (6) ;
               surftype = 'MLND : '
            case (7) ;
               surftype = 'MSNO : '
            end select
            write(unit=innov_rad_unit,fmt='(a,i7,9f10.2,3i3,f8.3,f10.2,f8.3)') surftype, n, &
                             iv%instid(i)%t2m(n), &
                             iv%instid(i)%mr2m(n),   &
                             iv%instid(i)%u10(n), &
                             iv%instid(i)%v10(n),  &
                             iv%instid(i)%ps(n),  &
                             iv%instid(i)%ts(n),  &
                             iv%instid(i)%smois(n),  &
                             iv%instid(i)%tslb(n),  &
                             iv%instid(i)%snowh(n), &
                             iv%instid(i)%isflg(n), &
                             nint(iv%instid(i)%soiltyp(n)), &
                             nint(iv%instid(i)%vegtyp(n)), &
                             iv%instid(i)%vegfra(n), &
                             iv%instid(i)%elevation(n), &
                             iv%instid(i)%clwp(n)

            write(unit=innov_rad_unit,fmt='(a)') 'OBS  : '
            write(unit=innov_rad_unit,fmt='(10f11.2)') ob%instid(i)%tb(:,n)
            write(unit=innov_rad_unit,fmt='(a)') 'BAK  : '
            write(unit=innov_rad_unit,fmt='(10f11.2)') iv%instid(i)%tb_xb(:,n)
            write(unit=innov_rad_unit,fmt='(a)') 'IVBC : '
            write(unit=innov_rad_unit,fmt='(10f11.2)')  iv%instid(i)%tb_inv(:,n)
            write(unit=innov_rad_unit,fmt='(a)') 'EMS  : '
            write(unit=innov_rad_unit,fmt='(10f11.2)')  iv%instid(i)%emiss(1:iv%instid(i)%nchan,n)
            if (rtm_option==rtm_option_crtm .and. write_jacobian) then
               write(unit=innov_rad_unit,fmt='(a)') 'EMS_JACOBIAN : '
               write(unit=innov_rad_unit,fmt='(10f10.3)') iv%instid(i)%emiss_jacobian(1:iv%instid(i)%nchan,n)
            end if
            write(unit=innov_rad_unit,fmt='(a)') 'ERR  : '
            write(unit=innov_rad_unit,fmt='(10f11.2)') iv%instid(i)%tb_error(:,n)
            write(unit=innov_rad_unit,fmt='(a)') 'QC   : '
            write(unit=innov_rad_unit,fmt='(10i11)') iv%instid(i)%tb_qc(:,n)

            if (write_profile) then
               nlevelss  = iv%instid(i)%nlevels
               if ( rtm_option == rtm_option_rttov ) then
               end if ! end if rtm_option_rttov

               if ( rtm_option == rtm_option_crtm ) then
                  write(unit=innov_rad_unit,fmt='(a)') &
                     'level fullp(mb) halfp(mb) t(k) q(g/kg) water(mm) ice(mm) rain(mm) snow(mm) graupel(mm) hail(mm)'
                  if (crtm_cloud) then
                     do k=1,iv%instid(i)%nlevels-1
                        write(unit=innov_rad_unit,fmt='(i3,2f10.2,f8.2,13f8.3)') &
                           k,  &
                           iv%instid(i)%pf(k,n), &
                           iv%instid(i)%pm(k,n), &
                           iv%instid(i)%tm(k,n), &
                           iv%instid(i)%qm(k,n), &
                           iv%instid(i)%qcw(k,n), &
                           iv%instid(i)%qci(k,n), &
                           iv%instid(i)%qrn(k,n), &
                           iv%instid(i)%qsn(k,n), &
                           iv%instid(i)%qgr(k,n), &
                           iv%instid(i)%qhl(k,n), &
                           iv%instid(i)%rcw(k,n), &
                           iv%instid(i)%rci(k,n), &
                           iv%instid(i)%rrn(k,n), &
                           iv%instid(i)%rsn(k,n), &
                           iv%instid(i)%rgr(k,n), &
                           iv%instid(i)%rhl(k,n)
                     end do ! end loop profile
                  else  ! no cloud
                     do k=1,iv%instid(i)%nlevels-1
                        write(unit=innov_rad_unit,fmt='(i3,2f10.2,f8.2,7f8.3)') &
                           k,  &
                           iv%instid(i)%pf(k,n), &
                           iv%instid(i)%pm(k,n), &
                           iv%instid(i)%tm(k,n), &
                           iv%instid(i)%qm(k,n), &
                           0.0, &
                           0.0, &
                           0.0, &
                           0.0, &
                           0.0, &
                           0.0
                     end do ! end loop profile
                  end if  ! end if crtm_cloud
               end if ! end if rtm_option_crtm

            end if  ! end if write_profile

            if ( rtm_option == rtm_option_crtm .and. write_jacobian) then

               if ( calc_weightfunc ) then
                  dtransmt(:,:)    = iv%instid(i)%der_trans(:,:,n)
                  transmt(:,:)     = iv%instid(i)%trans(:,:,n)
                  transmt_jac(:,:) = iv%instid(i)%trans_jacobian(:,:,n)
                  lod(:,:)         = iv%instid(i)%lod(:,:,n)
                  lod_jac(:,:)     = iv%instid(i)%lod_jacobian(:,:,n)
               else
                  dtransmt(:,:)    = 0.0
                  transmt(:,:)     = 0.0
                  transmt_jac(:,:) = 0.0
                  lod(:,:)         = 0.0
                  lod_jac(:,:)     = 0.0
               end if

               write(unit=innov_rad_unit,fmt='(a)') &
                  'channel level halfp(mb) t(k) q(g/kg) der_trans trans_jac trans lod_jac lod water(mm) ice(mm) rain(mm) snow(mm) graupel(mm) hail(mm)'
               if (crtm_cloud) then
                  do l=1,iv%instid(i)%nchan
                     do k=1,iv%instid(i)%nlevels-1
                        write(unit=innov_rad_unit,fmt='(i5,i3,f10.2,13f14.7,6f14.7)') &
                           l, k,  &
                           iv%instid(i)%pm(k,n), &
                           iv%instid(i)%t_jacobian(l,k,n), &
                           iv%instid(i)%q_jacobian(l,k,n), &
                           dtransmt(l,k),&
                           transmt_jac(l,k),&
                           transmt(l,k),&
                           lod_jac(l,k),&
                           lod(l,k),&
                           iv%instid(i)%water_jacobian(l,k,n), &
                           iv%instid(i)%ice_jacobian(l,k,n), & 
                           iv%instid(i)%rain_jacobian(l,k,n), &
                           iv%instid(i)%snow_jacobian(l,k,n), &
                           iv%instid(i)%graupel_jacobian(l,k,n), &
                           iv%instid(i)%hail_jacobian(l,k,n), &
                           iv%instid(i)%water_r_jacobian(l,k,n), &
                           iv%instid(i)%ice_r_jacobian(l,k,n), & 
                           iv%instid(i)%rain_r_jacobian(l,k,n), &
                           iv%instid(i)%snow_r_jacobian(l,k,n), &
                           iv%instid(i)%graupel_r_jacobian(l,k,n), &
                           iv%instid(i)%hail_r_jacobian(l,k,n)
                    end do ! end loop profile
                 end do ! end loop channels
               else  ! no cloud
                  do l=1,iv%instid(i)%nchan
                     do k=1,iv%instid(i)%nlevels-1
                        write(unit=innov_rad_unit,fmt='(i5,i3,f10.2,13f14.7,6f14.7)') &
                           l, k,  &
                           iv%instid(i)%pm(k,n), &
                           iv%instid(i)%t_jacobian(l,k,n), &
                           iv%instid(i)%q_jacobian(l,k,n), &
                           dtransmt(l,k),&
                           transmt_jac(l,k),&
                           transmt(l,k),&
                           lod_jac(l,k),&
                           lod(l,k),&
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0.
                     end do ! end loop profile
                  end do ! end loop channels
               end if  ! end if crtm_cloud
            end if !  end if write_jacobian

         end if ! end if proc_domain
      end do ! end do pixels
      if (rtm_option==rtm_option_crtm .and. write_jacobian ) then
         deallocate ( dtransmt )
         deallocate ( transmt_jac )
         deallocate ( transmt )
         deallocate ( lod )
         deallocate ( lod_jac )
      end if
      close(unit=innov_rad_unit)
      call da_free_unit(innov_rad_unit)
   end do ! end do instruments

   if (trace_use) call da_trace_exit("da_write_iv_rad_ascii")

end subroutine da_write_iv_rad_ascii

subroutine da_write_oa_rad_ascii (it, ob, iv, re )

   !---------------------------------------------------------------------------
   ! Purpose: write out OMB and OMA vector structure for radiance data.
   !---------------------------------------------------------------------------

   implicit none

   integer      ,     intent(in)  :: it       ! outer loop count
   type (y_type),     intent(in)  :: ob       ! Observation structure.
   type (iv_type),    intent(in)  :: iv       ! O-B structure.
   type (y_type),     intent(in)  :: re       ! O-A structure.

   integer                        :: n        ! Loop counter.
   integer                        :: i, k     ! Index dimension.
   integer                        :: nlevelss ! Number of obs levels.

   integer            :: ios, oma_rad_unit
   character(len=filename_len)  :: filename
   character(len=7)   :: surftype
   integer            :: ndomain
   logical            :: amsr2

   if (trace_use) call da_trace_entry("da_write_oa_rad_ascii")

   write(unit=message(1),fmt='(A)') 'Writing radiance OMA ascii file'
   call da_message(message(1:1))

   do i = 1, iv%num_inst
      if (iv%instid(i)%num_rad < 1) cycle

      ! count number of obs within the proc_domain
      !---------------------------------------------
      ndomain = 0
      do n =1,iv%instid(i)%num_rad
         if (iv%instid(i)%info%proc_domain(1,n)) then
            ndomain = ndomain + 1
         end if
      end do
      if (ndomain < 1) cycle

      amsr2 = index(iv%instid(i)%rttovid_string,'amsr2') > 0

      write(unit=filename, fmt='(i2.2,a,i4.4)') it,'_oma_'//trim(iv%instid(i)%rttovid_string)//'.', myproc

      call da_get_unit(oma_rad_unit)
      open(unit=oma_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_oa_rad_ascii.inc",49, &
            (/"Cannot open oma radiance file"//filename/))
      end if
      write(unit=oma_rad_unit,fmt='(a,a,i7,a,i5,a)') trim(iv%instid(i)%rttovid_string), &
                           ' number-of-pixels : ', ndomain, &
                           ' channel-number-of-each-pixel : ', iv%instid(i)%nchan, &
                           ' index-of-channels : '
      write(unit=oma_rad_unit,fmt='(10i5)') iv%instid(i)%ichan
      if ( amsr2 ) then
         write(unit=oma_rad_unit,fmt='(a)') ' pixel-info : i date scanpos landsea_mask  elv lat lon  satzen satazi clw'
      else
         write(unit=oma_rad_unit,fmt='(a)') ' pixel-info : i date scanpos landsea_mask  elv lat lon  satzen satazi '
      end if
      write(unit=oma_rad_unit,fmt='(a)') ' xb-surf-info : i t2m mr2m(ppmv) u10 v10 ps ts smois tslb snowh isflg &
                       & soiltyp vegtyp vegfra elev clwp'
      ndomain = 0
      do n=1,iv%instid(i)%num_rad
         if (iv%instid(i)%info%proc_domain(1,n)) then
            ndomain=ndomain+1
            if ( amsr2 ) then !write out clw
               write(unit=oma_rad_unit,fmt='(a,i7,2x,a,i6,i3,f6.0,4f8.2,f8.3)') 'INFO : ', ndomain, &
                                   iv%instid(i)%info%date_char(n), &
                                   iv%instid(i)%scanpos(n),        &
                                   iv%instid(i)%landsea_mask(n),   &
                                   iv%instid(i)%info%elv(n),       &
                                   iv%instid(i)%info%lat(1,n),     &
                                   iv%instid(i)%info%lon(1,n),     &
                                   iv%instid(i)%satzen(n),         &
                                   iv%instid(i)%satazi(n),         &
                                   iv%instid(i)%clw(n)
            else !no clw info
               write(unit=oma_rad_unit,fmt='(a,i7,2x,a,2i3,f6.0,4f8.2)') 'INFO : ', ndomain, &
                                   iv%instid(i)%info%date_char(n), &
                                   iv%instid(i)%scanpos(n),        &
                                   iv%instid(i)%landsea_mask(n),   &
                                   iv%instid(i)%info%elv(n),       &
                                   iv%instid(i)%info%lat(1,n),     &
                                   iv%instid(i)%info%lon(1,n),     &
                                   iv%instid(i)%satzen(n),         &
                                   iv%instid(i)%satazi(n)
            end if
            select case (iv%instid(i)%isflg(n))
            case (0) ;
               surftype = ' SEA : '
            case (1) ;
               surftype = ' ICE : '
            case (2) ;
               surftype = 'LAND : '
            case (3) ;
               surftype = 'SNOW : '
            case (4) ;
               surftype = 'MSEA : '
            case (5) ;
               surftype = 'MICE : '
            case (6) ;
               surftype = 'MLND : '
            case (7) ;
               surftype = 'MSNO : '
            end select
            write(unit=oma_rad_unit,fmt='(a,i7,9f10.2,3i3,f8.3,f10.2,f8.3)') surftype, n, &
                             iv%instid(i)%t2m(n), &
                             iv%instid(i)%mr2m(n),   &
                             iv%instid(i)%u10(n), &
                             iv%instid(i)%v10(n),  &
                             iv%instid(i)%ps(n),  &
                             iv%instid(i)%ts(n),  &
                             iv%instid(i)%smois(n),  &
                             iv%instid(i)%tslb(n),  &
                             iv%instid(i)%snowh(n), &
                             iv%instid(i)%isflg(n), &
                             nint(iv%instid(i)%soiltyp(n)), &
                             nint(iv%instid(i)%vegtyp(n)), &
                             iv%instid(i)%vegfra(n), &
                             iv%instid(i)%elevation(n), &
                             iv%instid(i)%clwp(n)

            write(unit=oma_rad_unit,fmt='(a)') 'OBS  : '
            write(unit=oma_rad_unit,fmt='(10f11.2)') ob%instid(i)%tb(:,n)
            write(unit=oma_rad_unit,fmt='(a)') 'BAK  : '
            write(unit=oma_rad_unit,fmt='(10f11.2)') iv%instid(i)%tb_xb(:,n)
            write(unit=oma_rad_unit,fmt='(a)') 'IVBC : '
            write(unit=oma_rad_unit,fmt='(10f11.2)')  iv%instid(i)%tb_inv(:,n)
            write(unit=oma_rad_unit,fmt='(a)') 'OMA  : '
            write(unit=oma_rad_unit,fmt='(10f11.2)')  re%instid(i)%tb(:,n)
            write(unit=oma_rad_unit,fmt='(a)') 'EMS  : '
            write(unit=oma_rad_unit,fmt='(10f11.2)')  iv%instid(i)%emiss(1:iv%instid(i)%nchan,n)
            write(unit=oma_rad_unit,fmt='(a)') 'ERR  : '
            write(unit=oma_rad_unit,fmt='(10f11.2)') iv%instid(i)%tb_error(:,n)
            write(unit=oma_rad_unit,fmt='(a)') 'QC   : '
            write(unit=oma_rad_unit,fmt='(10i11)') iv%instid(i)%tb_qc(:,n)

            if (write_profile) then
               nlevelss  = iv%instid(i)%nlevels
               if ( rtm_option == rtm_option_rttov ) then
               end if  ! end if rtm_option_rttov

               if ( rtm_option == rtm_option_crtm ) then
                  write(unit=oma_rad_unit,fmt='(a)') &
                     'level fullp(mb) halfp(mb) t(k) q(g/kg) water(mm) ice(mm) rain(mm) snow(mm) graupel(mm) hail(mm)'
                  if (crtm_cloud) then
                     do k=1,iv%instid(i)%nlevels-1
                        write(unit=oma_rad_unit,fmt='(i3,2f10.2,f8.2,13f8.3)') &
                           k,  &
                           iv%instid(i)%pf(k,n), &
                           iv%instid(i)%pm(k,n), &
                           iv%instid(i)%tm(k,n), &
                           iv%instid(i)%qm(k,n), &
                           iv%instid(i)%qcw(k,n), &
                           iv%instid(i)%qci(k,n), &
                           iv%instid(i)%qrn(k,n), &
                           iv%instid(i)%qsn(k,n), &
                           iv%instid(i)%qgr(k,n), &
                           iv%instid(i)%qhl(k,n), &
                           iv%instid(i)%rcw(k,n), &
                           iv%instid(i)%rci(k,n), &
                           iv%instid(i)%rrn(k,n), &
                           iv%instid(i)%rsn(k,n), &
                           iv%instid(i)%rgr(k,n), &
                           iv%instid(i)%rhl(k,n)
                     end do ! end loop profile
                  else  ! no cloud
                     do k=1,iv%instid(i)%nlevels-1
                        write(unit=oma_rad_unit,fmt='(i3,2f10.2,f8.2,7f8.3)') &
                           k,  &
                           iv%instid(i)%pf(k,n), &
                           iv%instid(i)%pm(k,n), &
                           iv%instid(i)%tm(k,n), &
                           iv%instid(i)%qm(k,n), &
                           0.0, &
                           0.0, &
                           0.0, &
                           0.0, &
                           0.0, &
                           0.0
                     end do ! end loop profile
                  end if ! end if crtm_cloud
               end if  ! end if crtm_option

            end if   ! end if write_profile
         end if    ! end if proc_domain
      end do     ! end do pixels
      close(unit=oma_rad_unit)
      call da_free_unit(oma_rad_unit)
   end do    !! end do instruments

   if (trace_use) call da_trace_exit("da_write_oa_rad_ascii")

end subroutine da_write_oa_rad_ascii


subroutine da_detsurtyp (snow,ice,landsea, vegtyp, soiltyp, is, ie, js, je, &
   i, j, dx, dy, dxm, dym, isflg,ob_vegtyp,ob_soiltyp, seap, icep, lndp, snop)

   !---------------------------------------------------------------------------
   ! Purpose: determine surface type at obs locations.
   !
   ! METHOD: using the background information 
   !      1. get the background landmask/snow/sea-ice
   !      2. calculate percentage of sea/ice/land/snow
   !      3. determine surface type at obs location
   !      4. find nearest grid vegtype and soil type 
   !
   !  HISTORY: 11/25/2008 - Use input surface type percentage      Tom Auligne
   !
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)  :: is, ie, js, je
   integer, intent(in)  :: i, j
   real   , intent(in)  :: dx, dxm, dy, dym
   real   , intent(in)  :: snow(is:ie,js:je)    ! Input variable
   real   , intent(in)  :: ice(is:ie,js:je)     ! Input variable
   real   , intent(in)  :: landsea(is:ie,js:je) ! Input variable
   integer, intent(in)  :: vegtyp(is:ie,js:je)  ! Input variable
   integer, intent(in)  :: soiltyp(is:ie,js:je) ! Input variable
   integer, intent(out) :: isflg                ! Output variable 
   real,    intent(out) :: ob_vegtyp            ! Output variable
   real,    intent(out) :: ob_soiltyp            ! Output variable
   real*8,  intent(out) :: seap, icep, lndp, snop ! percentage of surface type
   
   !     isflg    - surface flag
   !                0 sea
   !                1 sea ice
   !                2 land
   !                3 snow
   !                4 mixed predominately sea
   !                5 mixed predominately sea ice
   !                6 mixed predominately land
   !                7 mixed predominately snow

   !  local variables
!   integer   ::  n, xbflag(4)   ! surface type at xb location
!                               ! 0:sea  1:sea-ice 2:land  3:snow
   real      ::  w(4),minw      ! weight at 4 xb locations

   if (trace_use) call da_trace_entry("da_detsurtyp")    

   ! 1.0 determine surface type of xb at 4 location around obs
   !-------------------------------------------------------
!   if ( nint(landsea(i  ,j  )) == 0 ) xbflag(1) = 0 ! sea
!   if ( nint(landsea(i+1,j  )) == 0 ) xbflag(2) = 0 ! sea
!   if ( nint(landsea(i  ,j+1)) == 0 ) xbflag(3) = 0 ! sea
!   if ( nint(landsea(i+1,j+1)) == 0 ) xbflag(4) = 0 ! sea
!
!   if ( nint(landsea(i  ,j  )) == 1 ) xbflag(1) = 2 ! land/snow/sea-ice
!   if ( nint(landsea(i+1,j  )) == 1 ) xbflag(2) = 2 ! land/snow/sea-ice
!   if ( nint(landsea(i  ,j+1)) == 1 ) xbflag(3) = 2 ! land/snow/sea-ice
!   if ( nint(landsea(i+1,j+1)) == 1 ) xbflag(4) = 2 ! land/snow/sea-ice
!
!   if ( nint(snow(i  ,j  )) == 1 ) xbflag(1) = 3 ! snow
!   if ( nint(snow(i+1,j  )) == 1 ) xbflag(2) = 3 ! snow
!   if ( nint(snow(i  ,j+1)) == 1 ) xbflag(3) = 3 ! snow
!   if ( nint(snow(i+1,j+1)) == 1 ) xbflag(4) = 3 ! snow
!
!   if ( nint(ice(i  ,j  )) == 1 ) xbflag(1) = 1 ! sea-ice
!   if ( nint(ice(i+1,j  )) == 1 ) xbflag(2) = 1 ! sea-ice
!   if ( nint(ice(i  ,j+1)) == 1 ) xbflag(3) = 1 ! sea-ice
!   if ( nint(ice(i+1,j+1)) == 1 ) xbflag(4) = 1 ! sea-ice 

   ! 2.0 determine surface type percentage at obs location
   !------------------------------------------------------
   !  (i,j+1) 
   !    -----------------(i+1,j+1)
   !    |   w2          |
   !    |dym        w1  |
   !    |      obs      |
   !    |------- *      |
   !    |dy  w4  |  w3  |
   !    |   dx   | dxm  |
   !    |----------------
   !   (i,j)            (i+1,j)
   !
   !--------------------------------------------------------

!   seap = 0.0
!   icep = 0.0
!   snop = 0.0
!   lndp = 0.0
   w(1) = dym*dxm   ! weight for point (i,j)
   w(2) = dym*dx    ! weight for point (i+1,j)
   w(3) = dy *dxm   ! weight for point (i,j+1)
   w(4) = dy *dx    ! weight for point (i+1,j+1)

!   do n = 1, 4
!      if (xbflag(n) == 0) seap = seap+w(n)
!      if (xbflag(n) == 1) icep = icep+w(n)
!      if (xbflag(n) == 2) lndp = lndp+w(n)
!      if (xbflag(n) == 3) snop = snop+w(n)
!   end do
   
   lndp = 1.0/SUM(w) * ( w(1)*landsea(i,j)         + w(2)*landsea(i+1,j)          + &
                         w(3)*landsea(i,j+1)       + w(4)*landsea(i+1,j+1) )
   
   icep = 1.0/SUM(w) * ( w(1)*ice(i,j)              + w(2)*ice(i+1,j)              + &
                         w(3)*ice(i,j+1)            + w(4)*ice(i+1,j+1) )
   
   snop = 1.0/SUM(w) * ( w(1)*MIN(snow(i,j),   1.0) + w(2)*MIN(snow(i+1,j),   1.0) + &
                         w(3)*MIN(snow(i,j+1), 1.0) + w(4)*MIN(snow(i+1,j+1), 1.0) )
   
   if (icep > 0.0) then
      snop = MIN(snop, 1.0 - icep)
      seap = 1.0 - icep - snop
   else
      seap = 1.0 - lndp
   end if   
   lndp = MAX(1.0_8-seap-icep-snop, 0.0_8)
   
   ! fo2d   = dym*(dxm*fm2d(i,j  ) + dx*fm2d(i+1,j  )) &
   !         + dy *(dxm*fm2d(i,j+1) + dx*fm2d(i+1,j+1))

   ! 3.0 determine final surface flag at obs location
   !-----------------------------------------
   if (seap >= 0.99) isflg = 0
   if (icep >= 0.99) isflg = 1
   if (lndp >= 0.99) isflg = 2
   if (snop >= 0.99) isflg = 3
   if ( .not. (seap >= 0.99) .and. &
        .not. (icep >= 0.99) .and. &
        .not. (lndp >= 0.99) .and. &
        .not. (snop >= 0.99) ) then
      if (seap > lndp) then
         if (seap > icep) then
            if (seap > snop) then
               isflg = 4
            else
               isflg = 7
            end if
         else
            if (icep > snop) then
               isflg = 5
            else
               isflg = 7
            end if
         end if
      else
         if (lndp > icep) then
            if (lndp > snop) then
               isflg = 6
            else
               isflg = 7
            end if
         else
            if (icep > snop) then
               isflg = 5
            else
               isflg = 7
            end if
         end if
      end if 
   end if
   
   ! 4.0 find nearest grid vegtype and soil type
   !         at obs location
   !-----------------------------------------
   minw=min(w(1),w(2),w(3),w(4))
   if (minw == w(1)) then
      ob_vegtyp  = float(vegtyp (i+1,j+1))
      ob_soiltyp = float(soiltyp(i+1,j+1))
   else if (minw == w(2)) then
      ob_vegtyp  = float(vegtyp (i,j+1))
      ob_soiltyp = float(soiltyp(i,j+1))
   else if (minw == w(3)) then
      ob_vegtyp  = float(vegtyp (i+1,j))
      ob_soiltyp = float(soiltyp(i+1,j))
   else if (minw == w(4)) then
      ob_vegtyp  = float(vegtyp (i,j))
      ob_soiltyp = float(soiltyp(i,j))
   end if

   if (trace_use) call da_trace_exit("da_detsurtyp")    

end subroutine da_detsurtyp

subroutine da_cld_eff_radius(t,rho,qci,qrn,qsn,qgr,snow,xice,xland,method, &
                             reff_water,reff_ice,reff_rain,reff_snow,reff_grau)

   !---------------------------------------------------------------------------
   ! Purpose: compute effective radius of cloud liquid water and cloud ice
   !
   ! METHOD: 
   !   liquid water: adapted from WRFV2/phys/module_ra_cam.F, analytic formula following 
   !                 the formulation originally developed by J. T. Kiehl, and 
   !   ice (method 1): adapted from WRFV2/phys/module_ra_cam.F,  Kristjansson and Mitchell
   !   ice (method 2): WSM3, Hong et al., MWR 2004
   !   rain/snow/graupel: assume exponential particle size distribution and
   !                             spherical particles
   !                      use Gauss-Laguerre Quadrature for integration
   ! HISTORY: 12/15/2008 effective radius unit is micron.       Zhiquan Liu
   !---------------------------------------------------------------------------

   integer, intent(in)  :: method
   real*8,  intent(in)  :: t         ! temperature
   real,    intent(in)  :: rho       ! air density  kg/m3
   real*8,  intent(in)  :: qci       ! cloud ice mixing ratio kg/kg
   real*8,  intent(in)  :: qrn       ! cloud rain mixing ratio
   real*8,  intent(in)  :: qsn       ! cloud snow mixing ratio
   real*8,  intent(in)  :: qgr       ! cloud graupel mixing ratio
   real,    intent(in)  :: snow         ! snow water equivalent
   real*8,  intent(in)  :: xice         ! ice percentage
   real*8,  intent(in)  :: xland        ! landsea percentage
   real*8,  intent(out) :: reff_water   ! effective radius liquid water
   real*8,  intent(out) :: reff_ice     ! effective radius ice
   real*8,  intent(out) :: reff_rain    ! effective radius rain
   real*8,  intent(out) :: reff_snow    ! effective radius snow
   real*8,  intent(out) :: reff_grau    ! effective radius graupel

   !  local variables
   integer                      :: index, nk
   real                         :: snowh, corr
   real, parameter              :: rliqland  = 8.0     ! liquid drop size if over land
   real, parameter              :: rliqocean = 14.0    ! liquid drop size if over ocean
   real, parameter              :: rliqice   = 14.0    ! liquid drop size if over sea ice
   ! Tabulated values of re(T) in the temperature interval
   ! 180 K -- 274 K; hexagonal columns assumed:
   real, dimension(95), parameter  :: retab =                   &
         (/ 5.92779, 6.26422, 6.61973, 6.99539, 7.39234,        &
         7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,  &
         10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,  &
         15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,  &
         20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,  &
         27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,  &
         31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,  &
         34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,  &
         38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,  &
         42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,  &
         50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,  &
         65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,  &
         93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,  &
         124.954, 130.630, 136.457, 142.446, 148.608, 154.956,  &
         161.503, 168.262, 175.248, 182.473, 189.952, 197.699,  &
         205.728, 214.055, 222.694, 231.661, 240.971, 250.639 /)
   !
   ! constants from da_control.f:  pi,
   !
   real, parameter :: n0_rain  = 0.08      ! cm(-4)
   real, parameter :: n0_snow  = 0.04      ! cm(-4)
   real, parameter :: n0_grau  = 0.04      ! cm(-4)
   real, parameter :: rho_w    =  1000.0   ! kg m(-3)
   real, parameter :: rho_ice  =   900.0   ! kg m(-3)
   real, parameter :: rho_snow =   100.0   ! kg m(-3)
   real, parameter :: rho_grau =   400.0   ! kg m(-3)

   ! Abscissas of Gauss-Laguerre Integration
   real, dimension(32) :: xk = (/ 0.0444893658333, 0.23452610952,  &
      0.576884629302, 1.07244875382, 1.72240877644, 2.52833670643, &
      3.49221327285, 4.61645677223, 5.90395848335, 7.3581268086,   &
      8.98294126732, 10.783012089, 12.763745476, 14.9309117981,    &
      17.2932661372, 19.8536236493, 22.6357789624, 25.6201482024,  &
      28.8739336869, 32.3333294017, 36.1132042245, 40.1337377056,  &
      44.5224085362, 49.2086605665, 54.3501813324, 59.8791192845,  &
      65.9833617041, 72.6842683222, 80.1883747906, 88.735192639,   &
      98.8295523184, 111.751398227 /)

   ! total weights (weight*exp(xk)) of Gauss-Laguerre Integration
   real, dimension(32) :: totalw = (/ 0.114187105768, 0.266065216898, &
      0.418793137325, 0.572532846497, 0.727648788453, 0.884536718946, &
      1.04361887597, 1.20534920595, 1.37022171969, 1.53877595906,     &
      1.71164594592, 1.8895649683, 2.07318851235, 2.26590144444,      &
      2.46997418988, 2.64296709494, 2.76464437462, 3.22890542981,     &
      2.92019361963, 4.3928479809, 4.27908673189, 5.20480398519,      &
      5.11436212961, 4.15561492173, 6.19851060567, 5.34795780128,     &
      6.28339212457, 6.89198340969, 7.92091094244, 9.20440555803,     &
      11.1637432904, 15.3902417688 /)

   real, parameter :: limit = 1.0e-6
   real :: piover6   ! pi/6
   real :: sum1_rain, sum2_rain, sum1_snow, sum2_snow, sum1_grau, sum2_grau
   real :: lamda_rain, lamda_snow, lamda_grau
   real, dimension(32) :: psd_rain, psd_snow, psd_grau   ! partical size distribution

! initialize
   reff_water = 0.0
   reff_ice = 0.0
   reff_rain = 0.0
   reff_snow = 0.0
   reff_grau = 0.0
   lamda_rain = 0.0
   lamda_snow = 0.0
   lamda_grau = 0.0

! cloud liquid effective radius

   snowh = 0.001 * snow  ! here the snowh (in meter) is water-equivalent snow depth, 
                         ! which is different from the actually snow depth defined in the wrf output file
   
   ! Start with temperature-dependent value appropriate for continental air
   reff_water = rliqland + (rliqocean-rliqland) * min(1.0_8,max(0.0_8,(t_triple-t)*0.05_8))
   ! Modify for snow depth over land
   reff_water = reff_water + (rliqocean-reff_water) * min(1.0,max(0.0,snowh*10.))
   ! Ramp between polluted value over land to clean value over ocean.
   reff_water = reff_water + (rliqocean-reff_water) * min(1.0_8,max(0.0_8,1.0_8-xland))
   ! Ramp between the resultant value and a sea ice value in the presence of ice.
   reff_water = reff_water + (rliqice-reff_water) * min(1.0_8,max(0.0_8,xice))

! cloud ice effective radius

   if ( method == 1 ) then
      index = int(t-179.)
      index = min(max(index,1),94)
      corr = t - int(t)
      reff_ice = retab(index)*(1.-corr) + retab(index+1)*corr
   end if

! cloud ice effective radius
! rho*qci = 2.08*10**22 * Dice**8

   if ( method == 2 ) then
      ! 0.5 for diameter - radii conversion
      ! 1.0e6 for meter - micron conversion
      ! 0.125 = 1/8
      reff_ice = 1.0e6 * 0.5 * ( rho * qci / 2.08e22 ) ** 0.125
   end if
!
! cloud rain/snow/graupel effective radius
!
   piover6 = pi/6.
   if ( qrn > limit ) then
      lamda_rain = (piover6*rho_w*n0_rain*rho/qrn)**0.25
   end if
   if ( qsn > limit ) then
      lamda_snow = (piover6*rho_snow*n0_snow*rho/qsn)**0.25
   end if
   if ( qgr > limit ) then
      lamda_grau = (piover6*rho_grau*n0_grau*rho/qgr)**0.25
   end if
   sum1_rain = 0.0
   sum2_rain = 0.0
   sum1_snow = 0.0
   sum2_snow = 0.0
   sum1_grau = 0.0
   sum2_grau = 0.0
   if ( qrn > limit ) then
      do nk = 1, 32
         psd_rain(nk) = n0_rain*exp(-2.0*lamda_rain*xk(nk))
         sum1_rain = sum1_rain + totalw(nk) * (xk(nk)**3) * psd_rain(nk)
         sum2_rain = sum2_rain + totalw(nk) * (xk(nk)**2) * psd_rain(nk)
      end do
      reff_rain = 10000.0 * sum1_rain/sum2_rain    ! micron
   end if
   if ( qsn > limit ) then
      do nk = 1, 32
         psd_snow(nk) = n0_snow*exp(-2.0*lamda_snow*xk(nk))
         sum1_snow = sum1_snow + totalw(nk) * (xk(nk)**3) * psd_snow(nk)
         sum2_snow = sum2_snow + totalw(nk) * (xk(nk)**2) * psd_snow(nk)
      end do
      reff_snow = 10000.0 * sum1_snow/sum2_snow    ! micron
   end if
   if ( qgr > limit ) then
      do nk = 1, 32
         psd_grau(nk) = n0_grau*exp(-2.0*lamda_grau*xk(nk))
         sum1_grau = sum1_grau + totalw(nk) * (xk(nk)**3) * psd_grau(nk)
         sum2_grau = sum2_grau + totalw(nk) * (xk(nk)**2) * psd_grau(nk)
      end do
      reff_grau = 10000.0 * sum1_grau/sum2_grau    ! micron
   end if

end subroutine da_cld_eff_radius
subroutine da_ao_stats_rad ( stats_unit, iv, re )

   !---------------------------------------------------------------------------
   ! Purpose: Calculate statistics of obs minus analysis for radiance data.
   !
   ! METHOD:  average, rms, minimum, maximum of re
   !---------------------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! output unit for stats.
   type (iv_type), intent (inout) :: iv            ! innovation
   type (y_type),  intent (in)    :: re            ! o-a

   type (stats_rad_type), pointer  :: rad(:)
   integer                         :: n, k, i

   iv%nstats(radiance) = 0

   if ( iv%num_inst < 1 ) return

   if (trace_use) call da_trace_entry("da_ao_stats_rad")

   allocate ( rad(1:iv%num_inst) )

   do i = 1, iv%num_inst                          !! loop for sensors

      allocate (rad(i)%ichan(1:iv%instid(i)%nchan))
      rad(i)%ichan(1:iv%instid(i)%nchan)%num  = 0
      rad(i)%ichan(1:iv%instid(i)%nchan)%ave  = 0.0
      rad(i)%ichan(1:iv%instid(i)%nchan)%rms  = 0.0
      rad(i)%ichan(1:iv%instid(i)%nchan)%minimum%value  = -missing_r
      rad(i)%ichan(1:iv%instid(i)%nchan)%maximum%value  =  missing_r
      rad(i)%ichan(1:iv%instid(i)%nchan)%minimum%n      = 1
      rad(i)%ichan(1:iv%instid(i)%nchan)%maximum%n      = 1
      do k=1,iv%instid(i)%nchan
         rad(i)%ichan(k)%minimum%l      = k
         rad(i)%ichan(k)%maximum%l      = k
      end do

      if (iv%instid(i)%num_rad < 1) cycle

      do k=1, iv%instid(i)%nchan               !! loop for channels
         do n=1, iv%instid(i)%num_rad              !! loop for pixels
            if (iv%instid(i)%info%proc_domain(1,n)) then
               call da_stats_calculate( n,k,iv%instid(i)%tb_qc(k,n), &
                                 re%instid(i)%tb(k,n), rad(i)%ichan(k)%num, &
                                 rad(i)%ichan(k)%minimum, rad(i)%ichan(k)%maximum, &
                                 rad(i)%ichan(k)%ave, rad(i)%ichan(k)%rms)

            end if          ! end if( oi%sound(n)%loc%proc_domain )
         end do                                 !! end loop for pixels
      end do                        !  end loop for channels
   end do                         !  end loop for sensor

   do i = 1, iv%num_inst                          !! loop for sensors
      do k=1, iv%instid(i)%nchan               !! loop for channels
         ! FIX? generate 1D array to allow a da_proc_sum_ints call here
         ! Do inter-processor communication to gather statistics.
         call da_proc_sum_int ( rad(i)%ichan(k)%num )
         call da_proc_stats_combine(rad(i)%ichan(k)%ave, rad(i)%ichan(k)%rms, &
                           rad(i)%ichan(k)%minimum%value, rad(i)%ichan(k)%maximum%value, &
                           rad(i)%ichan(k)%minimum%n, rad(i)%ichan(k)%maximum%n, &
                           rad(i)%ichan(k)%minimum%l, rad(i)%ichan(k)%maximum%l )
 
         iv%nstats(radiance) = iv%nstats(radiance) + rad(i)%ichan(k)%num
      end do                        !  end loop for channels

      if (rootproc) then
         if (any(rad(i)%ichan(:)%num /= 0)) then
            write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for radiance         '//iv%instid(i)%rttovid_string
            call da_print_stats_rad( stats_unit, iv%instid(i)%nchan, rad(i) )
         end if
      end if
   end do                         !  end loop for sensor

   do i = 1, iv%num_inst           ! loop for sensors
      deallocate (rad(i)%ichan)
   end do

   deallocate (rad)

   if (trace_use) call da_trace_exit("da_ao_stats_rad")

end subroutine da_ao_stats_rad


subroutine da_oi_stats_rad ( stats_unit, iv )

   !---------------------------------------------------------------------------
   ! Purpose: Calculate statistics of obs minus background for radiance data.
   !
   ! METHOD:  average, rms, minimum, maximum of iv
   !---------------------------------------------------------------------------

   implicit none

   integer,        intent (in)      :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in)      :: iv            ! Innovation

   type (stats_rad_type),    pointer  :: rad(:)
   integer                          :: n, k, i

   if (trace_use) call da_trace_entry("da_oi_stats_rad")

   allocate ( rad(1:iv%num_inst) )

   do i = 1, iv%num_inst                          !! loop for sensors
      allocate ( rad(i)%ichan(1:iv%instid(i)%nchan) )
      rad(i)%ichan(1:iv%instid(i)%nchan)%num  = 0
      rad(i)%ichan(1:iv%instid(i)%nchan)%ave  = 0.0
      rad(i)%ichan(1:iv%instid(i)%nchan)%rms  = 0.0
      rad(i)%ichan(1:iv%instid(i)%nchan)%minimum%value  = -missing_r
      rad(i)%ichan(1:iv%instid(i)%nchan)%maximum%value  =  missing_r
      rad(i)%ichan(1:iv%instid(i)%nchan)%minimum%n      = 1
      rad(i)%ichan(1:iv%instid(i)%nchan)%maximum%n      = 1
      do k=1,iv%instid(i)%nchan
         rad(i)%ichan(k)%minimum%l      = k
         rad(i)%ichan(k)%maximum%l      = k
      end do

      if (iv%instid(i)%num_rad < 1) cycle
      do k=1, iv%instid(i)%nchan               !! loop for channels
         do n=1, iv%instid(i)%num_rad              !! loop for pixels
            if (iv%instid(i)%info%proc_domain(1,n)) then
               call da_stats_calculate( n,k,iv%instid(i)%tb_qc(k,n), &
                  iv%instid(i)%tb_inv(k,n), rad(i)%ichan(k)%num, &
                  rad(i)%ichan(k)%minimum, rad(i)%ichan(k)%maximum, &
                  rad(i)%ichan(k)%ave, rad(i)%ichan(k)%rms)

            end if          ! end if( iv%sound(n)%loc%proc_domain )
         end do                                 !! end loop for pixels
      end do                        !  end loop for channels
   end do                         !  end loop for sensor

   do i = 1, iv%num_inst                          !! loop for sensors
      do k=1, iv%instid(i)%nchan               !! loop for channels
         ! Do inter-processor communication to gather statistics.
         call da_proc_sum_int ( rad(i)%ichan(k)%num )
         call da_proc_stats_combine(rad(i)%ichan(k)%ave, rad(i)%ichan(k)%rms, &
            rad(i)%ichan(k)%minimum%value, rad(i)%ichan(k)%maximum%value, &
            rad(i)%ichan(k)%minimum%n, rad(i)%ichan(k)%maximum%n, &
            rad(i)%ichan(k)%minimum%l, rad(i)%ichan(k)%maximum%l )

      end do                        !  end loop for channels

      if (rootproc) then
         if ( any( rad(i)%ichan(:)%num /= 0 ) ) then
            write(unit=stats_unit, fmt='(/a/)') &
               ' Diagnostics of OI for radiance         '//iv%instid(i)%rttovid_string 
            call da_print_stats_rad( stats_unit, iv%instid(i)%nchan, rad(i) )
         end if
      end if
   end do                         !  end loop for sensor

   do i = 1, iv%num_inst           ! loop for sensors
      deallocate (rad(i)%ichan)
   end do

   deallocate (rad)

   if (trace_use) call da_trace_exit("da_oi_stats_rad")

end subroutine da_oi_stats_rad


subroutine da_print_stats_rad( stats_unit, nchan, rad )

   !---------------------------------------------------------------------------
   !  Purpose: print out statistics of omb, oma for radiance data.
   !
   !  METHOD:  print out average, rms, minimum, maximum of iv, re
   !---------------------------------------------------------------------------

   implicit none

   integer,           intent(in)    :: stats_unit, nchan
   type (stats_rad_type), intent(in)    :: rad
   
   integer    :: k, n, nmin, nmax
   integer    :: used_nchan 

   if (trace_use) call da_trace_entry("da_print_stats_rad")

   used_nchan = 0
   do k=1, nchan                  !! loop for channels
      if(rad%ichan(k)%num > 0) used_nchan = used_nchan + 1
   end do
   write(unit=stats_unit, fmt='((a,i5))')  ' used_nchan: ', used_nchan 

   write(unit=stats_unit, fmt='(6a)') &
        ' Channel ', &
        ' num  ', &
        ' ave  ', &
        ' rms  ', &
        ' min  ', &
        ' max  '

   do k=1, nchan                  !! loop for channels
      if (rad%ichan(k)%num > 0) then
         n    = rad%ichan(k)%num
         nmin = rad%ichan(k)%minimum%n
         nmax = rad%ichan(k)%maximum%n

         write(unit=stats_unit, fmt='((i3,i7,4f8.2))') &
            k, rad%ichan(k)%num, rad%ichan(k)%ave/real(n), &
            sqrt(rad%ichan(k)%rms/real(n)), &
            rad%ichan(k)%minimum%value, rad%ichan(k)%maximum%value
      end if
   end do

   if (trace_use) call da_trace_exit("da_print_stats_rad")

end subroutine da_print_stats_rad


subroutine da_qc_rad (it, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: perform quality control for radiance data.
   !
   ! METHOD:  separated QC for each sensor
   !---------------------------------------------------------------------------

   implicit none

   integer      ,  intent(in)      :: it         ! outer loop count
   type (y_type),  intent(in)      :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.

   integer :: i, nchan,p,j
   logical   :: amsua, amsub, hirs, msu,airs, hsb, ssmis, mhs, iasi, seviri
   logical   :: mwts, mwhs, atms, amsr2

   integer, allocatable :: index(:)
   integer :: num_tovs_avg
   integer, allocatable :: excess_count(:)
   integer, allocatable :: spare_count(:)
   integer :: transfer
   logical :: copy_found
   integer :: temp(num_procs)

   if (trace_use) call da_trace_entry("da_qc_rad")

   if ( .not. allocated(num_tovs_before) )  allocate (num_tovs_before(iv%num_inst,num_procs))
   if ( .not. allocated(num_tovs_after) )   allocate (num_tovs_after(iv%num_inst,num_procs))

   ! Cannot be more total send,receives than combination of processors
   if ( .not. allocated(tovs_copy_count) )  allocate (tovs_copy_count(iv%num_inst))
   if ( .not. allocated(tovs_send_pe) )     allocate (tovs_send_pe(iv%num_inst,num_procs*num_procs))
   if ( .not. allocated(tovs_recv_pe) )     allocate (tovs_recv_pe(iv%num_inst,num_procs*num_procs))
   if ( .not. allocated(tovs_send_start) )  allocate (tovs_send_start(iv%num_inst,num_procs*num_procs))
   if ( .not. allocated(tovs_send_count) )  allocate (tovs_send_count(iv%num_inst,num_procs*num_procs))
   if ( .not. allocated(tovs_recv_start) )  allocate (tovs_recv_start(iv%num_inst,num_procs*num_procs))

   call da_trace("da_qc_rad", message="allocated tovs redistibution arrays")

   if ( .not. allocated(index) )            allocate (index(num_procs))
   if ( .not. allocated(excess_count) )     allocate (excess_count(num_procs))
   if ( .not. allocated(spare_count) )      allocate (spare_count(num_procs))

   do i = 1, iv%num_inst

      !if (iv%instid(i)%info%n2 < iv%instid(i)%info%n1) cycle

      nchan    = iv%instid(i)%nchan

      amsua = trim(rttov_inst_name(rtminit_sensor(i))) == 'amsua'
      amsub = trim(rttov_inst_name(rtminit_sensor(i))) == 'amsub'
      hirs  = trim(rttov_inst_name(rtminit_sensor(i))) == 'hirs'
      msu   = trim(rttov_inst_name(rtminit_sensor(i))) == 'msu'
      airs  = trim(rttov_inst_name(rtminit_sensor(i))) == 'airs'
      hsb   = trim(rttov_inst_name(rtminit_sensor(i))) == 'hsb'
      ssmis = trim(rttov_inst_name(rtminit_sensor(i))) == 'ssmis'
      mhs   = trim(rttov_inst_name(rtminit_sensor(i))) == 'mhs'
      iasi  = trim(rttov_inst_name(rtminit_sensor(i))) == 'iasi'	  
      mwts  = trim(rttov_inst_name(rtminit_sensor(i))) == 'mwts'
      mwhs  = trim(rttov_inst_name(rtminit_sensor(i))) == 'mwhs'
      atms  = trim(rttov_inst_name(rtminit_sensor(i))) == 'atms'
      seviri = trim(rttov_inst_name(rtminit_sensor(i))) == 'seviri'
      amsr2 = trim(rttov_inst_name(rtminit_sensor(i))) == 'amsr2'
      if (hirs) then
         ! 1.0 QC for HIRS
         call da_qc_hirs(it, i,nchan,ob,iv)
      else if (airs) then
         call da_qc_airs(it, i,nchan,ob,iv)
      else if ( hsb ) then
         ! call da_qc_hsb(it, i,nchan,ob,iv)
         call da_warning("da_qc_rad.inc",73,(/'QC Not implemented for HSB'/))
      else if (amsua) then
         call da_qc_amsua(it,i,nchan,ob,iv)
      else if ( amsub ) then
         call da_qc_amsub(it,i,nchan,ob,iv)
      else if (msu) then
         ! call da_qc_msu(it, i,nchan, ob,iv)
         call da_warning("da_qc_rad.inc",80,(/'QC Not implemented for MSU'/))
      else if (ssmis) then
         call da_qc_ssmis(it, i,nchan,ob,iv)
      else if (mhs) then
         call da_qc_mhs(it,i,nchan,ob,iv)
      else if (iasi) then
         call da_qc_iasi(it,i,nchan,ob,iv)		 
      else if (mwhs) then
         call da_qc_mwhs(it,i,nchan,ob,iv)
      else if (mwts) then
         call da_qc_mwts(it,i,nchan,ob,iv)
      else if (atms) then
         call da_qc_atms(it,i,nchan,ob,iv)
      else if (seviri) then
         call da_qc_seviri(it,i,nchan,ob,iv)
      else if (amsr2) then
         call da_qc_amsr2(it,i,nchan,ob,iv)
      else
         write(unit=message(1),fmt='(A,A)') &
            "Unrecognized instrument",trim(rttov_inst_name(rtminit_sensor(i)))
         call da_error("da_qc_rad.inc",100,message(1:1))
      end if

      ! Report number of observations to other processors via rootproc

      num_tovs_before(i,:) = 0
      num_tovs_before(i,myproc+1)=iv%instid(i)%num_rad
      temp(:)= num_tovs_before(i,:)
      call da_proc_sum_ints(temp(:))

      call wrf_dm_bcast_integer(temp(:),num_procs)
      num_tovs_before(i,:) = temp(:)

      num_tovs_after(i,:) = num_tovs_before(i,:)

      if (rootproc .and. print_detail_rad) then
         write(unit=message(1),fmt='(A,I1,A)') "Instrument ",i, &
            " initial tovs distribution"
         write(unit=message(2),fmt=*) num_tovs_before(i,:)
         call da_message(message(1:2))
      end if

      ! Decide how to reallocate observations

      num_tovs_avg=sum(num_tovs_before(i,:))/num_procs

      call da_trace_int_sort(num_tovs_before(i,:),num_procs,index)

      do p=1,num_procs
         excess_count(p)=num_tovs_before(i,index(p))-num_tovs_avg
         spare_count(p)=num_tovs_avg-num_tovs_before(i,index(p))
      end do

      tovs_copy_count(i) = 0
      tovs_send_start(i,:) = 0
      tovs_send_count(i,:) = 0

      do
         copy_found = .false.
         do p=1,num_procs
            if (spare_count(p) > tovs_min_transfer) then
               do j=num_procs,1,-1
                  if (excess_count(j) > tovs_min_transfer) then
                     copy_found = .true.
                     tovs_copy_count(i)=tovs_copy_count(i)+1
                     tovs_send_pe(i,tovs_copy_count(i)) = index(j)-1
                     tovs_recv_pe(i,tovs_copy_count(i)) = index(p)-1
                     transfer=min(spare_count(p),excess_count(j))
                     tovs_send_count(i,tovs_copy_count(i)) = transfer
                     tovs_recv_start(i,tovs_copy_count(i)) = num_tovs_after(i,index(p))+1
                     num_tovs_after(i,index(p))=num_tovs_after(i,index(p))+transfer
                     num_tovs_after(i,index(j))=num_tovs_after(i,index(j))-transfer
                     tovs_send_start(i,tovs_copy_count(i)) = num_tovs_after(i,index(j))+1
                     spare_count(p)=spare_count(p)-transfer
                     excess_count(j)=excess_count(j)-transfer
                     exit
                  end if   
               end do
            end if
         end do
         if (.not. copy_found) exit
      end do   

      if (print_detail_rad) then
         write(unit=message(1),fmt='(A,I1,A)') "Instrument ",i," final tovs distribution"
         write(unit=message(2),fmt=*) num_tovs_after(i,:)
         call da_message(message(1:2))
      end if

      iv % instid(i) % num_rad_glo = sum(num_tovs_after(i,:))
   end do

   deallocate (index)
   deallocate (excess_count)
   deallocate (spare_count)

   if (trace_use) call da_trace_exit("da_qc_rad")

end subroutine da_qc_rad


subroutine da_setup_satcv(iv, be)

   !-----------------------------------------------------------------------
   ! Purpose: Set up satellite control variable
   !-----------------------------------------------------------------------

   implicit none
   
   type (iv_type), intent(inout) :: iv          ! Obs. increment structure.
   type (be_type), intent(inout) :: be          ! Background error structure.
   
   integer                      :: i, n, k, size_js, js_start, nclouds, ncv, kts_100hPa(1)
   integer                      :: satcv_size
   
   if (trace_use) call da_trace_entry("da_setup_satcv")
   
      size_js  = 0                        
      js_start = be % cv % size_jb + be % cv % size_je + be % cv % size_jp
      
      do i= 1, iv % num_inst
         allocate( iv%instid(i)%cv_index(iv%instid(i)%num_rad) )
         do n = 1, iv%instid(i)%num_rad
	    satcv_size = 0
	   
         ! Tskin
	 !------
	   if (use_satcv(1)) then
	      iv%instid(i)%cv_index(n)%ts = js_start + size_js + 1
              satcv_size = satcv_size + 1
	   end if
	      
	 ! Cloud Cover(s)  
	 !---------------
	   if (use_satcv(2)) then
              kts_100hPa = MAXLOC(iv%instid(i)%pm(kts:kte,n), &
                           MASK = iv%instid(i)%pm(kts:kte,n) < 100.0)
              nclouds    = kte - kts_100hPa(1) + 1
	      ncv        = nclouds !4
              allocate(iv%instid(i)%cv_index(n)%vtox(nclouds,nclouds))
              allocate(iv%instid(i)%cv_index(n)%cc(ncv))
	      iv%instid(i)%cv_index(n)%cc(:)   = (/ (js_start+size_js+satcv_size+k, k=1,ncv) /)
              iv%instid(i)%cv_index(n)%nclouds = nclouds
              iv%instid(i)%cv_index(n)%ncv     = ncv
              satcv_size = satcv_size + ncv
	   end if
	   
	   size_js = size_js + satcv_size
        end do
      end do
      
      be % cv % size_js = size_js
      cv_size_domain_js = size_js   
   
   if (trace_use) call da_trace_exit ("da_setup_satcv")

end subroutine da_setup_satcv
subroutine da_mspps_emis(tb, nchan, em)

! http://www.star.nesdis.noaa.gov/corp/scsb/mspps/algorithms.html
! land algorithms for emissivities at three AMSU channels (23.8, 31.4, 50.3 GHz)

   implicit none

   integer,                intent(in)  :: nchan
   real, dimension(nchan), intent(in)  :: tb
   real, dimension(nchan), intent(out) :: em

   real, dimension(3), parameter :: b0 = (/ -2.5404E-1,-2.2606E-1, 8.9494E-2 /)
   real, dimension(3), parameter :: b1 = (/  1.1326E-2, 3.4481E-3,-3.6615E-3 /)
   real, dimension(3), parameter :: b2 = (/ -1.9479E-5,-9.7185E-6,-4.2390E-7 /)
   real, dimension(3), parameter :: b3 = (/ -4.5763E-3, 4.3299E-3, 1.0636E-2 /)
   real, dimension(3), parameter :: b4 = (/  1.7833E-5, 5.3281E-6,-6.4559E-6 /)
   real, dimension(3), parameter :: b5 = (/  3.2324E-3, 1.8668E-3,-4.2449E-4 /)
   real, dimension(3), parameter :: b6 = (/ -1.9056E-5,-1.5369E-5,-6.6878E-6 /)
   real, parameter :: f1  = 23.8
   real, parameter :: f2  = 31.4
   real, parameter :: f3  = 50.3
   real, parameter :: f4  = 52.8
   real, parameter :: f5  = 53.596
   real, parameter :: f15 = 89.0
   real, parameter :: rmiss =   0.0
   real, parameter :: tbmin =  50.0
   real, parameter :: tbmax = 550.0

   integer :: k

   em = rmiss  ! initialize

   if ( tb(1) > tbmin .and. tb(1) < tbmax .and.  &
        tb(2) > tbmin .and. tb(2) < tbmax .and.  &
        tb(3) > tbmin .and. tb(3) < tbmax ) then
      do k = 1, 2
      ! do k = 1, 3
         em(k) = b0(k)+b1(k)*tb(1)+b2(k)*tb(1)**2+b3(k)*tb(2)+b4(k)*tb(2)**2 &
                +b5(k)*tb(3)+b6(k)*tb(3)**2
      end do
      em(3)=em(2)+(f3-f2)*(em(2)-em(1))/(f2-f1)   !linear interpolated
      em(4)=em(3)+(f4-f3)*(em(3)-em(2))/(f3-f2)   !linear interpolated
      em(5)=em(4)+(f5-f4)*(em(4)-em(3))/(f4-f3)
      em(15)=em(4)+(f15-f4)*(em(4)-em(3))/(f4-f3)
   end if

end subroutine da_mspps_emis

subroutine da_mspps_ts(tb, nchan, satzen, ts)
 
! http://www.star.nesdis.noaa.gov/corp/scsb/mspps/algorithms.html
! land algorithm for surface temperature
 
   use gsi_constants, only: deg2rad

   implicit none

   integer,                intent(in)  :: nchan
   real, dimension(nchan), intent(in)  :: tb
   real,                   intent(in)  :: satzen
   real,                   intent(out) :: ts

   real, parameter :: rmiss = -999.0
   real, parameter :: tbmin = 50.0
   real, parameter :: tbmax = 550.0

   real :: cza

   ts = rmiss  ! initialize

   if ( tb(1) > tbmin .and. tb(1) < tbmax .and.   &
        tb(2) > tbmin .and. tb(2) < tbmax .and.   &
        tb(3) > tbmin .and. tb(3) < tbmax .and.   &
        satzen >= 0.0 .and. satzen <= 90.0 ) then
      cza = COS(satzen*deg2rad)
      ts = 2.9079E2-(8.5059E-1-1.9821E-3*tb(1))*tb(1)+         &
           (6.1433E-1-2.3579E-3*tb(2))*tb(2)-                  &
           (1.1493-5.4709E-3*tb(3))*tb(3)-1.50E1*(cza-5.40E-1)
   end if

end subroutine da_mspps_ts



end module da_radiance1

