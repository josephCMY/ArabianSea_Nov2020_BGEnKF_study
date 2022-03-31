!=========================================================================
! Subroutine to construct all of the simulatated observations
! Written by Man-Yau Chan.
!=========================================================================
! This subroutine calls all of the other xb subroutines.
subroutine compute_yf(wrf_file, g_comm, s_comm, gid, sid, iid, jid, proj)

  ! Load modules
  use constants
  use namelist_define
  use wrf_field_indices_define
  use obs_define
  use mpi_module
  use netcdf
  use radar 
  use common_variables
  use boundary_variables

  implicit none

  ! Input variable definitions
  integer, intent(in) :: g_comm, s_comm, gid,sid, iid, jid
  character (len=10), intent(in) :: wrf_file
  type(proj_info), intent(in) :: proj

  ! Local subroutine variable definitions
  integer   :: ngx, ngz, iob_radmin, iob_radmax, iob_MWmin, iob_MWmax
  integer   :: i, j, k, m, n, iob, iiob, nob, ie, iunit=80011, ounit, ii, jj, kk, is, it, ig, iv, i1,j1, itot
  integer   :: ist,ied,jst,jed,kst,ked, istart,iend,jstart,jend, uist,uied,ujst,ujed
  integer   :: sist,sied,sjst,sjed, sistart,siend,sjstart,sjend,sis,sie,sjs,sje
  real      :: gaussdev, error, xb
  real, dimension(ni,nj,nk,nv)                   :: x_t
  real, dimension(3,3,nk,nv,nm) ::  xob

  real, dimension(nk,nv,nm)         :: xob_column
  real, dimension(3,3,nk,nv) :: tmp
  real, dimension(3):: center_xb

  character (len=10) :: filename, date, time, zone, varname
  character (len=10)  :: obstype
  double precision :: timer,t0, time1, time2, time_nonIR, time_IR
  real, dimension(obs%num,nm) :: yasend_tb
  real, dimension(obs%num) :: yasend_mw
  character (len=20) :: format1
  real, dimension( obs%num ) :: obs_lvl_send, obs_lvl_recv


  ! Output variable communicated via common_variables module: ya

  ! ---------------------------------------------------------------------------------------------


  ! Determine slab domain bounding indices (in the domain's coordinate system)
  ! --------------------------------------------------------------------------
  istart=iid*ni+1
  iend=(iid+1)*ni
  jstart=jid*nj+1
  jend=(jid+1)*nj
  if(iid==(nicpu-1)) iend=ix+1
  if(jid==(njcpu-1)) jend=jx+1
  if(iend<istart .or. jend<jstart) then
    if(m==1.and.gid==0) write(*,'(a,i6,a)') '*******domain is too small to be decomposed! set nicpu,njcpu smaller.********'
    stop
  endif


  ! Initialization of ya.
  ya=0.;  yasend=0.

  !! Get truth from fort.80010 if needed
  !! -----------------------------------
  !! this part not fully tested.
  !filename="fort.80010"
  !if ( use_simulated .or. use_ideal_obs ) then
  !   call rd_truth_wrf(filename,ix,jx,kx,ni,nj,nk,nv,nm,gid,sid,iid,jid,x_t)
  !   if ( expername .eq. 'hurricane ' ) &
  !      call hurricane_center_shift(filename,ix,jx,kx,xlong,xlat,znu,znw,proj,times)
  !endif

  time1 = MPI_Wtime()
  if (my_proc_id==0) then
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30)')'Running h(x) for gts obs:'
    close(17)
  endif

  !-------------------------------------------------------------------------------
  ! PHASE 1: Generating non-BT simulated observations
  !-------------------------------------------------------------------------------
  ! Parallelization scheme:
  !------------------------
  ! The relevant xb_to_* subroutines use 3x3 columns of the model to compute the
  ! simulated observations. The 3x3 columns for obs iob is centered on (in
  ! domain coordinates):
  !   ctr_x = int(obs%position(iob,1))+1
  !   ctr_y = int(obs%position(iob,2))+1
  !
  ! Easiest parallelization is for each slab to handle obs with (ctr_x, ctr_y)
  ! within the slab. When all slabs are done, we can use MPI_ALLREDUCE to gather
  ! the ensemble of simulated obs together and disseminate to all processes
  ! 
  ! Challenge is what to do with points that fall on the edge of the slab. The
  ! 3x3 columns will then overlap multiple slabs. 
  ! To overcome this issue, we can simply extend each slab's north-south and
  ! east-west boundaries by 1. 
  ! As a result, processes handling the edge obs do not need to talk to other
  ! processes to obtain their stuff. This eliminates synchronization.
  !
  ! Procedure:
  ! 1) Extend the slabs' boundaries by 1 on each side. Do this for both mean
  !    and member fields.
  ! 2) In each slab: perform level assignment and/or osse obs generation based
  !                  on the extended mean field.
  ! 3) In each slab: run the various xb_to_* routines using the extended member
  !                  fields
  ! 4) Communicate derived ya and obs%position(:,3) across all processes.
  !-------------------------------------------------------------------------------

  ! Note: step 1) has been written directly into read_ensemble_bcast.
  ! ABEI also inflated the extended slab boundary.

  time1=MPI_Wtime()

  ! 2) Assign obs model level, and generate OSSE obs (if needed)
  ! ------------------------------------------------------------
  ! Find obs location in z direction, return as obs%position(:,3)
  ! For simulated obs, calculate obs value.
  ! Note that no communication between threads occur here.

  ! This model level assignment is only executed once per EnKF execution. 
  ! I.e., during the batchwise data assimilation process, the level assignment
  ! should not be run. 
  level_assignment_switch: if ( flag_level_assignment ) then

    ! Initialize model lvls 
    do iob = 1, obs%num
      obstype = obs%type(iob)
      if (obstype .ne. 'Radiance  ') obs%position(iob,3) = 0.
    enddo

    iob_lvl_loop: do iob=1,obs%num

       obstype = obs%type(iob)

       ! Examine location of obs
       i1=int(obs%position(iob,1))
       j1=int(obs%position(iob,2))
       tmp=0.

       ! If center of 3x3 columns is not in slab, skip.
       if ( i1+1 < istart .or. i1+1 > iend    &
           .or. j1+1 < jstart .or. j1+1 > jend ) &
           cycle iob_lvl_loop


       ! If proc gets here, slab can handle obs site.

       ! Grab information from slab and form the 3x3 cols
       sn_lvl_loop: do j=j1,j1+2
         we_lvl_loop: do i=i1,i1+2
          
           ! If i and j are within the slab proper
           if (i>=istart.and.i<=iend.and.j>=jstart.and.j<=jend) &
             tmp(i-i1+1,j-j1+1,:,:)=xm(i-istart+1,j-jstart+1,:,:)

           ! Western boundary case
           if ( i == istart -1 .and. j>=jstart.and.j<=jend ) &
             tmp(i-i1+1,j-j1+1,:,:)=xm_we_bdy(1,j-jstart+1,:,:)

           ! Eastern boundary case
           if ( i == iend +1 .and. j>=jstart.and.j<=jend ) &
             tmp(i-i1+1,j-j1+1,:,:)=xm_we_bdy(2,j-jstart+1,:,:)

           ! Southern boundary case
           if ( j == jstart -1 .and. i>=istart.and.i<=iend ) &
             tmp(i-i1+1,j-j1+1,:,:)=xm_sn_bdy(1,i-istart+1,:,:)

           ! Northern boundary case
           if ( j == jend +1 .and. i>=istart.and.i<=iend ) &
             tmp(i-i1+1,j-j1+1,:,:)=xm_sn_bdy(2,i-istart+1,:,:)

           ! West-south corner
           if ( i == istart - 1 .and. j == jstart - 1 ) &
             tmp(i-i1+1,j-j1+1,:,:)=xm_crnrs(1,1,:,:)

           ! East-south corner
           if ( i == iend   + 1 .and. j == jstart - 1 ) &
             tmp(i-i1+1,j-j1+1,:,:)=xm_crnrs(2,1,:,:)

           ! West-north corner
           if ( i == istart - 1 .and. j == jend   + 1 ) &
             tmp(i-i1+1,j-j1+1,:,:)=xm_crnrs(1,2,:,:)

           ! East-north corner
           if ( i == iend   + 1 .and. j == jend   + 1 ) &
             tmp(i-i1+1,j-j1+1,:,:)=xm_crnrs(2,2,:,:)
           
         enddo we_lvl_loop
       enddo sn_lvl_loop

       ! Proceeding to assign model height        
       write(filename,'(a5,i5.5)') wrf_file(1:5), iunit+numbers_en+1-1
       if ( obstype(1:5) == 'Radar' ) then
          call xb_to_rv(filename,proj,tmp,ix,jx,kx,nv,iob,xlong,znw,xb,0) 
       else if ( obstype(1:1) == 'P' .or. obstype(1:1) == 'H'  ) then
          call xb_to_sounding(filename,proj,tmp,ix,jx,kx,nv,iob,xlong,znu,znw,p_top,xb,1,0)
       else if (obstype .ne. 'Radiance  ') then
          obs%position(iob,3) = 1.
       endif
    enddo iob_lvl_loop
 
    ! Inserting safety measure.
    call MPI_BARRIER( comm, ierr )

    if (my_proc_id==0) then
      write(*,'(a,x,f5.1,x,a)') 'Found gts zlvl in', MPI_Wtime()- time1, 'seconds.'
      open( 17, file='enkf_time.log', status='old',position='append', action='write')
      write(17,'(a30, x,f8.1,x,a)')'    Found gts zlvl in', MPI_Wtime()- time1, 'seconds.'
      close(17)
    endif

  endif level_assignment_switch


  time1=MPI_Wtime()
  ! 2) Compute simulated observations for ensemble
  ! ----------------------------------------------
  ! Each slab runs obs op for sites within the slab.

  if(my_proc_id==0) write(*,*) 'Calculating Hx...'

  timer = MPI_Wtime()
  
  iob_run_loop: do iob = 1, obs%num


    ! Special treatment for hurricane position and intensity obs (HPI)
    ! Hurricane position and intensity obs 
    obstype = obs%type(iob)
    HPI_obs_switch: if ( obstype == 'longitude' .or. obstype == 'latitude  ' .or. obstype == 'min_slp   ' ) then
      if (nprocs >= numbers_en+1 .and. my_proc_id+1 <= numbers_en ) then
        ie = my_proc_id+1
        write(filename,'(a5,i5.5)') wrf_file(1:5), 80010 + ie
        call hurricane_center_assimilation(filename,ix,jx,kx,int(obs%position(iob,1)),int(obs%position(iob,2)), &
                                           center_xb,znu,znw,xlong,xlat,proj)
        if ( obstype == 'longitude ' ) yasend(iob,ie) = xlong(center_xb(1),center_xb(2))
        if ( obstype == 'latitude  ' ) yasend(iob,ie) = xlat(center_xb(1),center_xb(2))
        if ( obstype == 'min_slp   ' ) yasend(iob,ie) = center_xb(3) - pr/100.
      endif
    endif HPI_obs_switch


    ! Prepare x profiles near obs position (xob) for xb calculation
    xob = 0.
    
    obstype = obs%type(iob)

    ! Examine location of obs
    i1=int(obs%position(iob,1))
    j1=int(obs%position(iob,2))

    



    ! If center of 3x3 columns is not in slab, skip.
    if ( i1+1 < istart .or. i1+1 > iend    &
        .or. j1+1 < jstart .or. j1+1 > jend ) &
        cycle iob_run_loop

    ! If proc gets here, slab can handle obs site.

    ! Grab information from slab and form the 3x3 cols
    sn_run_loop: do j=j1,j1+2
      we_run_loop: do i=i1,i1+2
       
        ! If i and j are within the slab proper
        if (i>=istart.and.i<=iend.and.j>=jstart.and.j<=jend) &
          xob(i-i1+1,j-j1+1,:,:,:)=x(i-istart+1,j-jstart+1,:,:,:)

        ! Western boundary case
        if ( i == istart -1 .and. j>=jstart.and.j<=jend ) &
          xob(i-i1+1,j-j1+1,:,:,:)=x_we_bdy(1,j-jstart+1,:,:,:)

        ! Eastern boundary case
        if ( i == iend +1 .and. j>=jstart.and.j<=jend ) &
          xob(i-i1+1,j-j1+1,:,:,:)=x_we_bdy(2,j-jstart+1,:,:,:)

        ! Southern boundary case
        if ( j == jstart -1 .and. i>=istart.and.i<=iend ) &
          xob(i-i1+1,j-j1+1,:,:,:)=x_sn_bdy(1,i-istart+1,:,:,:)

        ! Northern boundary case
        if ( j == jend +1 .and. i>=istart.and.i<=iend ) &
          xob(i-i1+1,j-j1+1,:,:,:)=x_sn_bdy(2,i-istart+1,:,:,:)

        ! West-south corner
        if ( i == istart - 1 .and. j == jstart - 1 ) &
          xob(i-i1+1,j-j1+1,:,:,:)=x_crnrs(1,1,:,:,:)

        ! East-south corner
        if ( i == iend   + 1 .and. j == jstart - 1 ) &
          xob(i-i1+1,j-j1+1,:,:,:)=x_crnrs(2,1,:,:,:)

        ! West-north corner
        if ( i == istart - 1 .and. j == jend   + 1 ) &
          xob(i-i1+1,j-j1+1,:,:,:)=x_crnrs(1,2,:,:,:)

        ! East-north corner
        if ( i == iend   + 1 .and. j == jend   + 1 ) &
          xob(i-i1+1,j-j1+1,:,:,:)=x_crnrs(2,2,:,:,:)
        
      enddo we_run_loop
    enddo sn_run_loop

    if(print_detail>4 .and. gid==0) &
       write(*,'(a,i6,3a,f10.2,a,3f10.2)') 'No.', iob,' ',obstype,'=', obs%dat(iob), ' at (x,y,z)=', obs%position(iob,1:3)

    ! Run the observation operator.

    ! Iterate over ensemble members held by proc
    ens_run_loop: do n = 1, nm
      ie=(n-1)*nmcpu+gid+1
      if (ie<=numbers_en+1) then

        ! Write file name for member
        write( filename, '(a5,i5.5)') wrf_file(1:5), iunit+ie-1

        ! Radar obs
        if ( obstype(1:5) == 'Radar' ) then
          call xb_to_rv(filename,proj,xob(:,:,:,:,n),ix,jx,kx,nv,iob,xlong,znw,yasend(iob,ie),obs%position(iob,3),1)


        ! Upper air in-situ obs
        else if ( obstype(1:1) == 'P' .or. obstype(1:1) == 'H'  ) then
          call xb_to_sounding (filename,proj,xob(:,:,:,:,n),ix,jx,kx,nv,iob,xlong,znu,znw,p_top,yasend(iob,ie),1,1)

        ! Surface obs
        else if ( obstype(1:1) == 'S' ) then
          call xb_to_surface(filename,proj,xob(:,:,:,:,n),ix,jx,kx,nv,iob,times,yasend(iob,ie))
        ! Sea-level obs
        else if ( obstype(1:3) == 'slp' ) then
          call xb_to_slp(filename,xob(:,:,:,:,n),ix,jx,kx,nv,iob,znu,znw,yasend(iob,ie))

        ! Precipitable water obs
        else if ( obstype(1:2) == 'pw' ) then
          call xb_to_pw(filename,xob(:,:,:,:,n),ix,jx,kx,nv,iob,znu,znw,yasend(iob,ie))

        ! Idealized obs
        else if ( obstype(1:5) == 'ideal' ) then
          call xb_to_idealsound(filename,xob(:,:,:,:,n),ix,jx,kx,nv,iob,yasend(iob,ie))
        endif
      endif
    enddo ens_run_loop

  enddo iob_run_loop

  call MPI_BARRIER(COMM, IERR)
  if (my_proc_id==0) then
    write(*,'(a,x,f5.1,x,a)') 'Ran GTS h(x) in', MPI_Wtime()- time1, 'seconds.'
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30, x,f8.1,x,a)')'    Ran GTS h(x) in', MPI_Wtime()- time1, 'seconds.'
    close(17)
  endif
  time1=MPI_Wtime()

  ! 4) Communicate ya and obs%position(:,3) across all processes
  ! -----------------------------------------------------------
  ! ya is communicated across all procs
  call MPI_ALLREDUCE( yasend, ya, obs%num*(numbers_en+1), MPI_REAL, &
                      MPI_SUM, comm, ierr )


  ! Obs position is communicated across procs with the same gid.
  level_assignment_communication_switch: if ( flag_level_assignment ) then
      obs_lvl_send = 0.
      do iob = 1, obs%num
        obstype = obs%type(iob)
        if (obstype .ne. 'Radiance  ') obs_lvl_send(iob) = obs%position(iob,3)
      enddo
      obs_lvl_recv = 0.

      call MPI_ALLREDUCE( obs_lvl_send, obs_lvl_recv, obs%num, MPI_REAL, MPI_SUM, &
                          s_comm, ierr)
      ! Fill in all non-Radiance obs lvls
      do iob = 1, obs%num
        obstype = obs%type(iob)
        if (obstype .ne. 'Radiance  ') obs%position(iob,3) = obs_lvl_recv(iob)
      enddo
  endif level_assignment_communication_switch

  if (my_proc_id==0) then
    write(*,'(a,x,f5.1,x,a)') 'Communicated GTS h(x) in', MPI_Wtime()- time1, 'seconds.'
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30, x,f8.1,x,a,/)')'    Communicated GTS in', MPI_Wtime()- time1, 'seconds.'
    close(17)
  endif
  time1=MPI_Wtime()


  ! Printing HPI
  if (my_proc_id == 0) then
    do iob = 1, obs%num
      obstype = obs%type(iob)
      if  ( obstype == 'longitude' .or. obstype == 'latitude  ' .or. obstype == 'min_slp   ' ) &
        write(*,*) obstype, ya(iob,:) 
    enddo
  endif



  !-------------------------------------------------------------------------------
  ! PHASE 2: Generating ensemble's IR-BT
  !-------------------------------------------------------------------------------

  ! Computing infrared brightness temperatures
  if(raw%radiance%num.ne.0) then
    ! Each thread runs CRTM for obs locations within slab
    ! Note: each thread can contain multiple members (potentially)
    yasend_tb=0.
    call xb_to_radiance(x,proj,ix,jx,kx,ni,nj,nk,nv,nm,xlong,xlat,xland,iob_radmin,iob_radmax,yasend_tb,iid,jid,nicpu,njcpu,.true.)
    if (my_proc_id == 0) then
      write(*,*) 'Info about CRTM'
      write(*,'(a,i5,i5)') 'iob_radmin, iob_radmax', iob_radmin, iob_radmax
    endif
    do n = 1, nm
      ie = (n-1)*nmcpu+gid+1 
      yasend(iob_radmin:iob_radmax,ie) = yasend_tb(iob_radmin:iob_radmax,n)
    enddo
    

    call MPI_Barrier( comm, ierr)
    if (  my_proc_id == 0 ) then
      write(*,*) "Survived running xb_to_radiance"
      write(*,'(a,x,f5.1,x,a)') 'Finished IR h(x) in', MPI_Wtime()- time1, 'seconds.'
      open( 17, file='enkf_time.log', status='old',position='append', action='write')
      write(17,'(a30, x,f8.1,x,a,/)')'Finished IR h(x) in', MPI_Wtime()- time1, 'seconds.'
      close(17)
    endif
  endif

  ! Computing microwave brightness temperatures
  if(raw%microwave%num.ne.0) then
    yasend_mw=0.
    do ie = 1, numbers_en+1
      yasend_mw = 0.0
      write( filename, '(a5,i5.5)') wrf_file(1:5), iunit+ie-1
      call xb_to_microwave(filename,proj,ix,jx,kx,xlong,xlat,xland,times,iob_MWmin,iob_MWmax,yasend_mw,.true.,.true.)
      yasend(iob_MWmin:iob_MWmax,ie) = yasend_mw(iob_MWmin:iob_MWmax)
    enddo
    call MPI_Barrier( comm, ierr)
    if (  my_proc_id == 0 ) then
      write(*,*) "Survived running xb_to_microwave"
      write(*,'(a,x,f5.1,x,a)') 'Finished MW h(x) in', MPI_Wtime()- time1, 'seconds.'
      open( 17, file='enkf_time.log', status='old',position='append', action='write')
      write(17,'(a30, x,f8.1,x,a,/)')'Finished MW h(x) in', MPI_Wtime()- time1, 'seconds.'
      close(17)
    endif
  endif
  
  call MPI_Allreduce(yasend,ya,obs%num*(numbers_en+1),MPI_REAL,MPI_SUM,comm,ierr)
  time1=MPI_Wtime()
  !exit
  !if (use_elf) call
  !MPI_Allreduce(yasend_nc,ya_nc,obs%num,MPI_REAL,MPI_SUM,comm,ierr)
  !ya_ca = (abs(obs%dat-ya_nc)+abs(ya(:,numbers_en+1)-ya_nc))/2.
  
  call MPI_Barrier(comm, ierr)
  
  !write(*,'(a,i3,a,f,i)') 'Proc ', my_proc_id, ' radiance ', ya(iob_radmin+7000,
  !2), numbers_en
  !
  !
  call MPI_Barrier( comm, ierr)
  
  if ( my_proc_id == 0 ) then
    write( *,*) 'Finished computing Hx. Now computing <Hx>.'
  endif
  ! Compute linear mean of simulated observations
  ! Man-Yau Chan, Dec 29th 2019
  ! Blanking out the means to be safe
  ya( :, numbers_en+1) = 0.0
  ! Compute linear mean
  do ie=1, numbers_en
    ya( :, numbers_en+1) = ya( :,numbers_en+1) + ya( : , ie) / float(  numbers_en )
  enddo
  
  if ( my_proc_id == 0 ) then
    write( *,*) 'Finished computing <Hx>.'
    write(*,'(a,x,f5.1,x,a)') 'Finished computing <Hx>', MPI_Wtime()- time1, 'seconds.'
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30, x,f8.1,x,a,/)')'Finished computing <Hx>', MPI_Wtime()- time1, 'seconds.'
    close(17)
  endif


  !calcurate mean of ya (radmean) by Minamide 2015.9.25 > 2016.12.6
  !yfm_radiance = 0
  !do ie = 1, numbers_en
  !  yfm_radiance = yfm_radiance + ya(:,ie)/float(numbers_en)
  !enddo
  !yam_radiance = yfm_radiance
  
  write(format1,'(a,i4,a)') '(i5,a,',numbers_en+1,'f10.2)'
  do iob=1,obs%num
    if(print_detail>4 .and. my_proc_id==0) write(*,format1) iob, ' '//obs%type(iob)//' ya=',ya(iob,:)!numbers_en+1)
  enddo
  if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' Calculation of y=Hx tooks ', MPI_Wtime()-timer, ' seconds.'


end subroutine compute_yf







!=======================================================================================
subroutine xb_to_surface(inputfile,proj,xa,ix,jx,kx,nv,iob,times,xb)
use constants
use namelist_define
use obs_define
use netcdf
use wrf_tools
use map_utils
use mpi_module
use wrf_field_indices_define
implicit none
character(len=10), intent(in)           :: inputfile
character(len=19), intent(in)           :: times
character(len=10) :: obstype
type(proj_info), intent(in)             :: proj                   ! 1st guest map info
real, dimension(3,3,kx+1,nv), intent(in) :: xa                     ! 1st guest
integer, intent(in)                     :: ix, jx, kx, nv, iob
real, intent(out)                       :: xb
real, dimension(3, 3      )           :: t2m, th2m, q2m, u10m, v10m, psfc
integer                                 :: i1, j1, k1, i, j, k, m, ii, jj, kk, obs_ii,obs_jj
integer                                 :: i_t2, i_th2, i_q2, i_u10, i_v10, i_psfc
real                                    :: dx, dxm, dy, dym, dz, dzm
real                                    :: u10, v10, t2, q2, th2, ps

! Getting obs info and interpolation weights
obstype=obs%type(iob)
obs_ii=obs%position(iob,1)
obs_jj=obs%position(iob,2)
i1 = int( obs_ii )
j1 = int( obs_jj )
dx  = obs_ii-real(i1)
dxm = real(i1+1)-obs_ii
dy  = obs_jj-real(j1)
dym = real(j1+1)-obs_jj

! Recentering coordiantes
i1 = 1
j1 = 1

i_t2 = ind_t2
i_th2 = ind_th2
i_q2 = ind_q2
i_u10 = ind_u10
i_v10 = ind_v10
i_psfc = ind_psfc
!do m = 1, nv
!   if ( enkfvar(m) == 'T2        ' ) i_t2=m
!   if ( enkfvar(m) == 'TH2       ' ) i_th2=m
!   if ( enkfvar(m) == 'Q2        ' ) i_q2=m
!   if ( enkfvar(m) == 'U10       ' ) i_u10=m
!   if ( enkfvar(m) == 'V10       ' ) i_v10=m
!   if ( enkfvar(m) == 'PSFC      ' ) i_psfc=m
!enddo
if(i_t2>0)  t2m (i1:i1+1, j1:j1+1)  = xa( 1:2, 1:2, 1, i_t2 )
if(i_th2>0) th2m (i1:i1+1, j1:j1+1) = xa( 1:2, 1:2, 1, i_th2 )
if(i_q2>0)  q2m (i1:i1+1, j1:j1+1)  = xa( 1:2, 1:2, 1, i_q2 )
if(i_u10>0) u10m (i1:i1+1, j1:j1+1) = xa( 1:2, 1:2, 1, i_u10 )
if(i_v10>0) v10m (i1:i1+1, j1:j1+1) = xa( 1:2, 1:2, 1, i_v10 )
if(i_psfc>0) psfc (i1:i1+1, j1:j1+1)= xa( 1:2, 1:2, 1, i_psfc )
!if ( i_t2== 0  ) call get_variable2d(inputfile, 'T2        ', ix, jx, 1,    t2m )
!if ( i_th2== 0 ) call get_variable2d(inputfile, 'TH2       ', ix, jx, 1,    th2m)
!if ( i_q2== 0  ) call get_variable2d(inputfile, 'Q2        ', ix, jx, 1,    q2m )
!if ( i_u10== 0 ) call get_variable2d(inputfile, 'U10       ', ix, jx, 1,    u10m)
!if ( i_v10== 0 ) call get_variable2d(inputfile, 'V10       ', ix, jx, 1,    v10m)
!if ( i_psfc== 0) call get_variable2d(inputfile, 'PSFC      ', ix, jx, 1,    psfc)

t2  = dym*(dx*t2m (i1+1,j1) + dxm*t2m (i1,j1)) + dy*(dx*t2m (i1+1,j1+1) + dxm*t2m (i1, j1+1))
th2 = dym*(dx*th2m(i1+1,j1) + dxm*th2m(i1,j1)) + dy*(dx*th2m(i1+1,j1+1) + dxm*th2m(i1, j1+1))
q2  = dym*(dx*q2m (i1+1,j1) + dxm*q2m (i1,j1)) + dy*(dx*q2m (i1+1,j1+1) + dxm*q2m (i1, j1+1))
u10 = dym*(dx*v10m(i1+1,j1) + dxm*u10m(i1,j1)) + dy*(dx*u10m(i1+1,j1+1) + dxm*u10m(i1, j1+1))
v10 = dym*(dx*u10m(i1+1,j1) + dxm*v10m(i1,j1)) + dy*(dx*v10m(i1+1,j1+1) + dxm*v10m(i1, j1+1))
ps  = dym*(dx*psfc(i1+1,j1) + dxm*psfc(i1,j1)) + dy*(dx*psfc(i1+1,j1+1) + dxm*psfc(i1, j1+1))


!-------------------------------------------------
if ( obstype(10:10) == 'P' ) then
     xb = ps
else if ( obstype(10:10) == 'U' ) then
     xb = u10
else if ( obstype(10:10) == 'V' ) then
     xb = v10
else if ( obstype(10:10) == 'T' ) then
     xb = th2
else if ( obstype(10:10) == 'K' ) then
     xb = t2
else if ( obstype(10:10) == 'Q' ) then
     xb = q2*1000.
else if ( obstype(10:10) == 'D' ) then
     xb = mixrat_to_tdew(q2, ps)
else if ( obstype(10:10) == 'R' ) then
     xb = rel_humidity(q2, t2, ps)
else if ( obstype(10:10) == 'S' ) then   !wind speed
     xb = sqrt(u10**2.+v10**2.)
endif


end subroutine xb_to_surface

!=======================================================================================
subroutine xb_to_sounding (inputfile,proj,xa,ix,jx,kx,nv,iob,xlong,znu,znw,p_top,xb,itruth,isimulated )
   use constants
   use namelist_define
   use obs_define
   use netcdf
   use wrf_tools
   use map_utils
   use mpi_module
   use wrf_field_indices_define

   implicit none

   character(len=10), intent(in)           :: inputfile
   type(proj_info), intent(in)             :: proj                   ! 1st guest map info
   real, dimension(3,3,kx+1,nv), intent(in)  :: xa                     ! 1st guest
   integer, intent(in)                     :: ix, jx, kx, nv, iob  
   real, dimension(ix, jx ), intent(in)    :: xlong
   real, dimension(kx+1), intent(in)       :: znw
   real, dimension(kx), intent(in)         :: znu
   real,  intent(in)                       :: p_top
   integer, intent(in)                     :: itruth, isimulated

   real :: obs_ii, obs_jj, obs_kk, obs_pp !
   real, intent(out)   :: xb
   character(len=10)   :: obstype

   !real, dimension(ix, jx, kx+1)           :: ph, phb
   !real, dimension(ix, jx, kx  )           :: pt, pb, qv, qc, qr
   !real, dimension(ix, jx      )           :: mu, mub
   !real, dimension(ix+1, jx, kx)           :: u
   !real, dimension(ix, jx+1, kx)           :: v
   !real, dimension(kx+1)                   :: znw0
   !real, dimension(kx)                     :: znu0

   real, dimension(2, 2, kx+1)           :: ph, phb
   real, dimension(2, 2, kx  )           :: pt, pb, qv, qc, qr
   real, dimension(2, 2      )           :: mu, mub
   real, dimension(3, 2, kx)           :: u
   real, dimension(2, 3, kx)           :: v
   real, dimension(kx+1)                   :: znw0
   real, dimension(kx)                     :: znu0

   real, dimension(2, 2, kx+1)             :: p
   real, dimension(kx+1)                   :: work, worku, workv

   real, dimension(kx)                     :: pres, ptt, qvt
   real                                    :: mu1, mub1, long, grid_u, grid_v, true_u, true_v
   real                                    :: dx, dxm, dy, dym, dz, dzm
   integer                                 :: i1, j1, k1, i, j, k, m, ii, jj, kk

   obstype=obs%type(iob)

   ! Computing obs position and the interpolation weights
   obs_ii=obs%position(iob,1)
   obs_jj=obs%position(iob,2)
   obs_kk=obs%position(iob,3)
   obs_pp=obs%position(iob,4)
   i1 = int( obs_ii )
   j1 = int( obs_jj ) 
   dx  = obs_ii-real(i1)
   dxm = real(i1+1)-obs_ii
   dy  = obs_jj-real(j1)
   dym = real(j1+1)-obs_jj

   ! Switching over to indices relative to first guess (xa)
   i1 = 1
   j1 = 1

!.. calculate pressure from geopotential height
    if ( itruth == 0 )then
         call get_variable1d(inputfile, 'ZNW       ', kx+1, 1, znw0)
         call get_variable1d(inputfile, 'ZNU       ', kx  , 1, znu0)
    else
         znw0 = znw 
         znu0 = znu
    endif 

!.. Calculate obs%position(iob,3)
    if ( isimulated == 0 .or. obstype(10:10) == 'T' .or. obstype(10:10) == 'D' .or.  &
                              obstype(10:10) == 'R' .or. obstype(10:10) == 'Q' ) then
   
!..    get data from background
       if(ind_ph>0) &
         ph(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,ind_ph)
       if(ind_phb>0) &
         phb(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,ind_phb)
       if(ind_t>0) &
         pt(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,ind_t)
       if(ind_qvapor>0) &
         qv(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,ind_qvapor)
       if(ind_qcloud>0) &
         qc(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,ind_qcloud)
       if(ind_qrain>0) & 
         qr(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,ind_qrain)
       if(ind_mu>0)  &
         mu(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,ind_mu)
       if(ind_mub>0) &
         mub(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,ind_mub)

       ! Interesting, we consider qv as the total mixing ratio 
       qv (i1:i1+1, j1:j1+1, 1:kx) = qv (i1:i1+1, j1:j1+1, 1:kx) + qc (i1:i1+1, j1:j1+1, 1:kx) + qr (i1:i1+1, j1:j1+1, 1:kx)

!...... get mut mu, mub at obs' position from horizontal interpolation
       mu1 = dym*(dx*mu(i1+1,j1  ) + dxm*mu(i1,j1  )) + dy*(dx*mu(i1+1,j1+1) + dxm*mu(i1,j1+1))
       mub1 = dym*(dx*mub(i1+1,j1  ) + dxm*mub(i1,j1  )) + dy*(dx*mub(i1+1,j1+1) + dxm*mub(i1,j1+1))

!........ get qvt (qvapor) at obs' position from horizontal interpolation
!........ get ptt(theta)  at obs' position from horizontal interpolation
       qvt(1:kx) = dym*(dx*qv(i1+1,j1,1:kx) + dxm*qv(i1,j1,1:kx)) + dy*(dx*qv(i1+1,j1+1,1:kx) + dxm*qv(i1,j1+1,1:kx))
       ptt(1:kx) = dym*(dx*pt(i1+1,j1,1:kx) + dxm*pt(i1,j1,1:kx)) + dy*(dx*pt(i1+1,j1+1,1:kx) + dxm*pt(i1,j1+1,1:kx))

       p(1:2,1:2,1:kx+1) = ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1)
       p(1,1,1:kx+1) = dym*(dx*p(2,1,1:kx+1) + dxm*p(1,1,1:kx+1)) + dy*(dx*p(2,2,1:kx+1) + dxm*p(1,2,1:kx+1))

       call eta_to_pres(znw0(1:kx+1), mu1+mub1, qvt(1:kx), p(1,1,1:kx+1), ptt(1:kx)+to, kx, pres(1:kx))
!   if( print_detail==3 )write(*,'(3x,a,4f)')'obs location =', obs%position(iob,1:4)
!   if( print_detail==3 )write(*,'(3x,a,f)')'          mu =', mu(1,1)
!   if( print_detail==3 )write(*,'(3x,a,10f)')'          p  =', p(1,1,1:10)
!       call cal_press_from_q( kx, znu0, znw0, qvt, mu1, mub1, p_top, pres )
!       if( print_detail > 100 ) write(*,'(a,100f10.0)')'Pressure = ', pres
!       if( print_detail > 100 ) write(*,'(a,100f10.4)')'qvt = ', qvt

       if ( isimulated == 0 ) then
          if ( obstype(1:1) == 'H' ) then
!........... get geopotential height profile around obs%position, then horizontal interpolated to obs's position
             call to_zk(obs%position(iob,4), p(1,1,1:kx)/g, obs%position(iob,3), kx) 
          else if ( obstype(1:1) == 'P' ) then
             call to_zk(obs%position(iob,4), pres(1:kx), obs%position(iob,3), kx)
          endif
!          if ( obs%position(iob,3) .ge. real(kx-1) ) obs%position(iob,3) = real(kx-1)
          if ( obs%position(iob,3) .lt. 1. ) obs%position(iob,3) = 1.
          obs_kk = obs%position(iob,3)
       endif
    endif



!.. Get xb from background
   if ( obs_kk .gt. real(kx-1) .or. obs_kk .lt. 1. ) then
      xb = -999999.
      return
   endif

   k1  = int( obs_kk )
   dz  = obs_kk-real(k1)
   dzm = real(k1+1)-obs_kk

!.. U, V
    if ( obstype(10:10) == 'U' .or. obstype(10:10) == 'V' .or. obstype(10:10) == 'S' ) then
       u(i1:i1+2,j1:j1+1,k1:k1+1)=xa(1:3,1:2,k1:k1+1,ind_u)
       v(i1:i1+1,j1:j1+2,k1:k1+1)=xa(1:2,1:3,k1:k1+1,ind_v)

       worku = -88888.
       workv = -88888.
       if ( dx == 0.500 ) then
          worku(k1:k1+1) = dym*u(i1+1,j1,k1:k1+1) + dy*u(i1+1,j1+1,k1:k1+1)
       else if ( dx > 0.500 ) then
          worku(k1:k1+1) = dym*( u(i1+1,j1,k1:k1+1)*(dxm+0.5)+u(i1+2,j1,k1:k1+1)*(dx-0.5) ) + &
                           dy *( u(i1+1,j1+1,k1:k1+1)*(dxm+0.5)+u(i1+2,j1+1,k1:k1+1)*(dx-0.5) )
       else if ( dx < 0.500 ) then
          worku(k1:k1+1) = dym*( u(i1,j1,k1:k1+1)*(dxm-0.5)+u(i1+1,j1,k1:k1+1)*(dx+0.5) ) + &
                           dy *( u(i1,j1+1,k1:k1+1)*(dxm-0.5)+u(i1+1,j1+1,k1:k1+1)*(dx+0.5) )
       endif
       if ( dy == 0.500 ) then
          workv(k1:k1+1) = v(i1,j1+1,k1:k1+1)*dxm + v(i1+1,j1+1,k1:k1+1)*dx
       else if ( dy > 0.500 ) then
          workv(k1:k1+1) = dxm*( v(i1,j1+1,k1:k1+1)*(dym+0.5) + v(i1,j1+2,k1:k1+1)*(dy-0.5) ) + &
                           dx *( v(i1+1,j1+1,k1:k1+1)*(dym+0.5) + v(i1+1,j1+2,k1:k1+1)*(dy-0.5) )
       else if ( dy < 0.500 ) then
          workv(k1:k1+1) = dxm*( v(i1,j1,k1:k1+1)*(dym-0.5) + v(i1,j1+1,k1:k1+1)*(dy+0.5) ) + &
                           dx *( v(i1+1,j1,k1:k1+1)*(dym-0.5) + v(i1+1,j1+1,k1:k1+1)*(dy+0.5) )
       endif

!       if( print_detail > 100 ) write(*,*)'i1,j1,k1, v=',i1,j1,k1,v(i1,j1,k1),v(i1,j1+1,k1),v(i1,j1+2,k1),v(i1+1,j1,k1)
!       if( print_detail > 100 ) write(*,'(a,100f10.2)')'U profile =', worku(k1:k1+1)
!       if( print_detail > 100 ) write(*,'(a,100f10.2)')'V profile =', workv(k1:k1+1)

       if ( obs_kk .le. 1. ) then
          grid_u = worku(k1)
          grid_v = workv(k1)
       else
          grid_u = dzm*worku(k1)+dz*worku(k1+1)
          grid_v = dzm*workv(k1)+dz*workv(k1+1)
       endif
!       if( print_detail > 100 ) write(*,'(a,f10.2)')'grid U =', grid_u, 'grid V =', grid_v

!------011108 changed obs. wind to gridwind
       long = ( xlong(i1  ,j1)*dym + xlong(i1  ,j1+1)*dy ) * dxm +          & 
              ( xlong(i1+1,j1)*dym + xlong(i1+1,j1+1)*dy ) * dx
       call gridwind_to_truewind( long, proj, grid_u, grid_v, true_u, true_v )
    
       if ( obstype(10:10) == 'U' ) then
           xb = grid_u
       else if ( obstype(10:10) == 'V' ) then
           xb = grid_v
       else if ( obstype(10:10) == 'S' ) then   !wind speed
           xb = sqrt(grid_u**2.+grid_v**2.)
       endif

!       if( print_detail > 100 ) then
!           write(*,'(12f10.3)')obs_kk, grid_u, grid_v, dz, dzm, long, worku(k1:k1+1), workv(k1:k1+1), true_u, true_v
!       endif

!.. T
    else if ( obstype(10:10) == 'T' ) then
       do k = k1, k1+1
          work(k) = theta_to_temp(ptt(k)+to, pres(k))
       enddo
       if ( obs_kk .le. 1. ) then
          xb = work(1)
       else
          xb = dzm*work(k1)+dz*work(k1+1)
       endif

!.. TD(if used TD, not RH)
    else if ( obstype(10:10) == 'D' ) then
       do k = k1, k1+1
          work(k) = mixrat_to_tdew(qvt(k), pres(k))
       enddo
       if ( obs_kk .le. 1. ) then
          xb = work(1)
       else
          xb = dzm*work(k1) + dz*work(k1+1)
       endif

!.. RH
    else if ( obstype(10:10) == 'R' )then
       do k = k1, k1+1
          work(k) = theta_to_temp(ptt(k)+to, pres(k))
       enddo
       do k = k1, k1+1
          work(k) = rel_humidity(qvt(k), work(k), pres(k))
       enddo
       if ( obs_kk .le. 1. ) then
          xb = work(1)
       else
          xb = dzm*work(k1) + dz*work(k1+1)
       endif
       if ( xb > 100. )xb = 100.

!.. Q
    else if ( obstype(10:10) == 'Q' )then
       if ( obs_kk .le. 1. ) then
          xb = qvt(k1)*1000.
       else
          xb = (dzm*qvt(k1)+dz*qvt(k1+1))*1000.
       endif
       if ( xb .le. 0.0 ) xb = -99999.

!.. HG
    else if ( obstype(10:10) == 'H' )then
       call destag_zstag(znu0, znw0, kx, p(1,1,1:kx+1), work(1:kx))
       if ( obs_kk .le. 1. ) then
          xb = work(1)/g
       else
          xb = (dzm*work(k1) + dz*work(k1+1))/g
       endif 

    endif

!   if( print_detail > 100 )write(*,'(a,3f8.2, f10.1)')'xb_to_sounding '//obstype//' obs position :', obs_ii, obs_jj, obs%position(iob,3:4)
    
end subroutine xb_to_sounding
!
!=======================================================================================
   subroutine xb_to_idealsound(inputfile,xa,ix,jx,kx,nv,iob,xb)
   use constants
   use namelist_define
   use obs_define 
   use netcdf 
   use wrf_tools
   use map_utils 

   implicit none

   character(len=10), intent(in)           :: inputfile
   real, dimension(3,3,kx+1,nv), intent(in) :: xa 
   integer, intent(in) :: ix, jx, kx, nv, iob
   real :: obs_ii, obs_jj, obs_kk 
   character(len=10) :: obstype
   real, intent(out) :: xb
   integer :: i1, j1, k1, i, j, k, m, ii, jj, kk

   obstype=obs%type(iob)
   obs_ii=obs%position(iob,1)
   obs_jj=obs%position(iob,2)
   obs_kk=obs%position(iob,3)
   i1 = int( obs_ii )
   j1 = int( obs_jj )
   k1 = int( obs_kk )
    
   if (obstype .eq. 'idealU    ' ) then
      do m = 1, nv
         if ( enkfvar(m) == 'U         ' ) then
            xb = xa(i1,j1,k1,m)
         endif
      enddo
   else if (obstype .eq. 'idealV    ' ) then
      do m = 1, nv
         if (enkfvar(m) == 'V         ' ) then
            xb = xa(i1,j1,k1,m)
         endif
      enddo 
   else if (obstype .eq. 'idealPT   ' ) then
      do m = 1, nv
         if (enkfvar(m) == 'T         ' ) then
            xb = xa(i1,j1,k1,m) 
         endif
      enddo
   else if (obstype .eq. 'idealQV   ' ) then
      do m = 1, nv
         if (enkfvar(m) == 'QVAPOR    ' ) then
            xb = xa(i1,j1,k1,m)*1000.
         endif
      enddo
   else if (obstype .eq. 'idealPH   ' ) then
      xb = 0.
      do m = 1, nv
         if ( enkfvar(m) == 'PH        ' .or. enkfvar(m) == 'PHB       ' ) then
            xb = xb + xa(i1,j1,k1,m)
         endif
      enddo
   endif

   end subroutine xb_to_idealsound

!=======================================================================================
  subroutine xb_to_rv(inputfile,proj,xa,ix,jx,kx,nv,iob,xlong,znw,xb,kkflag)
  use constants 
  use namelist_define
  use mapinfo_define
  use obs_define
  use map_utils 
  use netcdf
  use wrf_tools
  use radar
  implicit none
  character (len=10), intent(in) :: inputfile
  type(proj_info), intent(in) :: proj
  integer, intent(in)         :: ix, jx, kx, nv, iob
  real, dimension(3,3,kx+1,nv), intent(in)  :: xa
  real, dimension(ix, jx ), intent(in) :: xlong
  real, dimension(kx+1), intent(in)    :: znw
  real, intent(out)                    :: xb
  integer, intent(in)                  :: kkflag
  real  :: obs_ii, obs_jj, obs_kk, obs_hh, radar_ii, radar_jj, radar_elv

  real, dimension(ix+1, jx  , kx)        :: u
  real, dimension(ix  , jx+1, kx)        :: v
  real, dimension(ix  , jx  , kx+1)        :: w, ph, phb
  real, dimension(ix  , jx  , kx  )        :: t, qv, p, pb
  real, dimension(ix  , jx      )        :: mu, mub, terrain
    
  integer :: m, ii, jj, kk, i, j, k, i1, j1
  integer :: i_u, i_v, i_w, i_t, i_qv, i_mu, i_ph, i_phb, i_mub, i_pb

  obs_ii=obs%position(iob,1)
  obs_jj=obs%position(iob,2)
  obs_kk=obs%position(iob,3)
  obs_hh=obs%position(iob,4)

  radar_ii=obs%sta(iob,1)
  radar_jj=obs%sta(iob,2)
  radar_elv=obs%sta(iob,4)

  i1 = int( obs_ii )
  j1 = int( obs_jj )
    
  i_u  = 0 ; i_v   = 0 ; i_w  = 0 ; i_t = 0 ; i_qv = 0 ; i_mu = 0; i_mub = 0;
  i_ph = 0 ; i_phb = 0 ; i_pb = 0 ;

  do m=1,nv
     if ( enkfvar(m) == 'U         ' ) i_u=m
     if ( enkfvar(m) == 'V         ' ) i_v=m
     if ( enkfvar(m) == 'W         ' ) i_w=m
  enddo
  if(i_u>0) u(i1:i1+2,j1:j1+1,1:kx)=xa(1:3,1:2,1:kx,i_u)
  if(i_v>0) v(i1:i1+1,j1:j1+2,1:kx)=xa(1:2,1:3,1:kx,i_v)
  if(i_w>0) w(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_w)
  if ( i_u == 0 ) call get_variable3d(inputfile, 'U         ', ix+1, jx, kx, 1, u )
  if ( i_v == 0 ) call get_variable3d(inputfile, 'V         ', ix, jx+1, kx, 1, v )
  if ( i_w == 0 ) call get_variable3d(inputfile, 'W         ', ix, jx, kx+1, 1, w )
  
  if ( kkflag == 0 ) then
     do m = 1, nv
        if ( enkfvar(m) == 'T         ' ) i_t=m
        if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
        if ( enkfvar(m) == 'PH        ' ) i_ph=m
        if ( enkfvar(m) == 'PHB       ' ) i_phb=m
        if ( enkfvar(m) == 'PB        ' ) i_pb=m
        if ( enkfvar(m) == 'MU        ' ) i_mu=m
        if ( enkfvar(m) == 'MUB       ' ) i_mub=m
     end do
     if(i_t>0) t(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_t)
     if(i_qv>0) qv(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qv)
     if(i_pb>0) pb(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_pb)
     if(i_ph>0) ph(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_ph)
     if(i_phb>0) phb(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_phb)
     if(i_mu>0) mu(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mu)
     if(i_mub>0) mub(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mub)
     if( i_t <1  ) call get_variable3d(inputfile, 'T         ', ix, jx, kx, 1, t )
     if( i_qv <1 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx, 1, qv )
     if( i_pb <1 ) call get_variable3d(inputfile, 'PB        ', ix, jx, kx, 1, pb )
     if( i_ph <1 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
     if( i_phb <1) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb )
     if( i_mu <1 ) call get_variable2d(inputfile, 'MU        ', ix, jx, 1, mu )
     if( i_mub <1) call get_variable2d(inputfile, 'MUB       ', ix, jx, 1, mub )
     do j = j1, j1+1
     do i = i1, i1+1
        call cal_ph( kx, znw, t(i,j,1:kx), qv(i,j,1:kx), pb(i,j,1:kx), mu(i,j), mub(i,j), phb(i,j,1:kx+1), ph(i,j,1:kx+1) )
     enddo
     enddo

!     if (print_detail > 1000) write(*,'(a,70f8.0)') 'phb+ph = ',(phb(i1,j1,1:kx+1)+ph(i1,j1,1:kx+1))/g
!     if (print_detail > 1000) write(*,'(a,f8.0)') 'obs_hh =', obs_hh
     call calc_radar_data_point_kk ( ix, jx, kx, ph, phb, obs_ii, obs_jj, obs_hh, obs%position(iob,3) )
     obs_kk=obs%position(iob,3)
!     if (print_detail > 1000) write(*,'(a, 4f8.2)')'obs_ijhk =', obs_ii, obs_jj, obs_hh, obs_kk
  endif

  call calculate_rv ( ix, jx, kx, u, v, w, xlong, proj, obs_ii, obs_jj, obs_kk, obs_hh, radar_ii, radar_jj, radar_elv, xb )
!  write(*,*)'i,j,k =', obs_ii, obs_jj, obs_kk
!  write(*,*)'u,v,w =', u(int(obs_ii),int(obs_jj),int(obs_kk)), v(int(obs_ii),int(obs_jj),int(obs_kk)), w(int(obs_ii),int(obs_jj),int(obs_kk))

  end subroutine xb_to_rv

!=======================================================================================
subroutine xb_to_pw(inputfile,xa,ix,jx,kx,nv,iob,znu,znw,xb)
  use constants
  use namelist_define
  use obs_define
  use netcdf
  use wrf_tools
  use wrf_field_indices_define
  implicit none
  character(len=10), intent(in)           :: inputfile
  character(len=10) :: obstype
  integer, intent(in)                     :: ix, jx, kx, nv, iob
  real, dimension(3,3,kx+1,nv), intent(in) :: xa                     ! 1st guest
  real, dimension(kx), intent(in)         :: znu
  real, dimension(kx+1), intent(in)       :: znw
  real, intent(out)                       :: xb
  real, dimension(2, 2)                   :: pw,temp,vtemp,zdiff
  real, dimension(2+1, 2  , kx  )         :: u
  real, dimension(2  , 2+1, kx  )         :: v
  real, dimension(2  , 2  , kx  )           :: t,p,pb,z,qv
  real, dimension(2  , 2  , kx+1)           :: ph,phb
  integer                                 :: itot, m, i, j, k, i1, j1, i2, j2, search_radiu
  integer                                 :: i_t, i_u, i_v, i_ph, i_phb, i_qv, i_p, i_pb
  real                                    :: obs_ii, obs_jj, dx,dxm,dy,dym

! Loading obs info and interpolation factor
obstype=obs%type(iob)
obs_ii=obs%position(iob,1)
obs_jj=obs%position(iob,2)
i1 = int( obs_ii )
j1 = int( obs_jj )
dx  = obs_ii-real(i1)
dxm = real(i1+1)-obs_ii
dy  = obs_jj-real(j1)
dym = real(j1+1)-obs_jj


! Recenter coordinates
i1 = 1
j1 = 1
i_qv = ind_qvapor 
i_ph = ind_ph
i_phb = ind_phb
i_t = ind_t
i_p = ind_p
i_pb = ind_pb
!do m = 1, nv
!   if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
!   if ( enkfvar(m) == 'P         ' ) i_p=m
!   if ( enkfvar(m) == 'PB        ' ) i_pb=m
!   if ( enkfvar(m) == 'PH        ' ) i_ph=m
!   if ( enkfvar(m) == 'PHB       ' ) i_phb=m
!   if ( enkfvar(m) == 'T         ' ) i_t=m
!enddo
if(i_qv>0) qv (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qv )
if(i_ph>0) ph (i1:i1+1, j1:j1+1, 1:kx+1) = xa( 1:2,1:2,1:kx+1,i_ph )
if(i_phb>0) phb (i1:i1+1, j1:j1+1, 1:kx+1) = xa( 1:2,1:2,1:kx+1,i_phb )
if(i_p>0) p (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_p )
if(i_pb>0) pb (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_pb )
if(i_t>0) t (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_t )
     
!if ( i_t  == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, t)
!if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
!if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
!if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
!if ( i_p  == 0 ) call get_variable3d(inputfile, 'P         ', ix, jx, kx,   1, p  )
!if ( i_pb == 0 ) call get_variable3d(inputfile, 'PB        ', ix, jx, kx,   1, pb )

ph(i1:i1+1, j1:j1+1, 1:kx+1) = (ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1))/g
p(i1:i1+1, j1:j1+1, 1:kx)  = p(i1:i1+1, j1:j1+1, 1:kx) + pb(i1:i1+1, j1:j1+1, 1:kx)
do j=j1,j1+1
do i=i1,i1+1
  z(i, j, 1:kx) = (ph(i, j, 1:kx)*(znw(2:kx+1)-znu(1:kx))+ph(i, j, 2:kx+1)*(znu(1:kx)-znw(1:kx)))/(znw(2:kx+1)-znw(1:kx))
enddo
enddo
do j=j1,j1+1
do i=i1,i1+1
do k=1,kx
  if( p(i,j,k)/100. >= 1200. ) write(*,'(a,3i4,2f8.0)')'P error at:',i,j,k,p(i,j,k)/100.,pb(i,j,k)/100.
enddo
enddo
enddo
pw = 0.
do k=1,kx-1
  temp=(t(i1:i1+1,j1:j1+1,k)+300.)*(p(i1:i1+1,j1:j1+1,k)/Pr)**(rd/cp)
  vtemp=(1+0.61*qv(i1:i1+1,j1:j1+1,k))*temp
  zdiff=z(i1:i1+1,j1:j1+1,k+1)-z(i1:i1+1,j1:j1+1,k)
  pw=pw+p(i1:i1+1,j1:j1+1,k)/(rd*vtemp)*zdiff*qv(i1:i1+1,j1:j1+1,k)
enddo
xb = dym*(dx*pw(2,1) + dxm*pw(1,1)) + dy*(dx*pw(2,2) + dxm*pw(1,2))
xb = xb/10

end subroutine xb_to_pw

!=======================================================================================
subroutine xb_to_slp(inputfile,xa,ix,jx,kx,nv,iob,znu,znw,xb )

!---------------------
! slp_center subroutine calculates sea level pressure and find the hurricane center
!---------------------
  use constants
  use namelist_define
  use obs_define
  use netcdf
  use wrf_tools
  use wrf_field_indices_define

  implicit none

  character(len=10), intent(in)           :: inputfile
  character(len=10) :: obstype
  integer, intent(in)                     :: ix, jx, kx, nv, iob
  real, dimension(3,3,kx+1,nv), intent(in) :: xa                     ! 1st guest
  real, dimension(kx), intent(in)         :: znu
  real, dimension(kx+1), intent(in)       :: znw
  real, intent(out)                       :: xb

  real, dimension(2, 2)                   :: slp
  !real, dimension(ix+1, jx, kx  )         :: u
  !real, dimension(ix, jx+1, kx  )         :: v
  !real, dimension(ix, jx, kx  )           :: t
  !real, dimension(ix, jx, kx+1)           :: ph
  !real, dimension(ix, jx, kx+1)           :: phb
  !real, dimension(ix, jx, kx  )           :: p
  !real, dimension(ix, jx, kx  )           :: z
  !real, dimension(ix, jx, kx  )           :: pb
  !real, dimension(ix, jx, kx  )           :: qv, qc, qr
  real, dimension(2+1, 2  , kx  )           :: u
  real, dimension(2  , 2+1, kx  )           :: v
  real, dimension(2  , 2  , kx  )           :: t
  real, dimension(2  , 2  , kx+1)           :: ph
  real, dimension(2  , 2  , kx+1)           :: phb
  real, dimension(2  , 2  , kx  )           :: p
  real, dimension(2  , 2  , kx  )           :: z
  real, dimension(2  , 2  , kx  )           :: pb
  real, dimension(2  , 2  , kx  )           :: qv, qc, qr

  integer                                 :: itot, m, i, j, k, i1, j1, i2, j2, search_radiu
  real                                    :: obs_ii, obs_jj, dx,dxm,dy,dym

! Load basic info and interpolation weights
obstype=obs%type(iob)
obs_ii=obs%position(iob,1)
obs_jj=obs%position(iob,2)
i1 = int( obs_ii )
j1 = int( obs_jj )
dx  = obs_ii-real(i1)
dxm = real(i1+1)-obs_ii
dy  = obs_jj-real(j1)
dym = real(j1+1)-obs_jj

! Switching over to indices relative to first guess
i1=1
j1=1


if(ind_qvapor>0) qv (i1:i1+1, j1:j1+1, 1:kx)      = xa( 1:2, 1:2, 1:kx, ind_qvapor )
if(ind_qcloud>0) qc (i1:i1+1, j1:j1+1, 1:kx)      = xa( 1:2, 1:2, 1:kx, ind_qcloud )
if(ind_qrain>0)  qr (i1:i1+1, j1:j1+1, 1:kx)      = xa( 1:2, 1:2, 1:kx, ind_qrain )
if(ind_ph>0)     ph (i1:i1+1, j1:j1+1, 1:kx+1)    = xa( 1:2,1:2,1:kx+1, ind_ph )
if(ind_phb>0)    phb (i1:i1+1, j1:j1+1, 1:kx+1)   = xa( 1:2,1:2,1:kx+1, ind_phb )
if(ind_p>0)      p (i1:i1+1, j1:j1+1, 1:kx)       = xa( 1:2,1:2,1:kx,   ind_p )
if(ind_pb>0)     pb (i1:i1+1, j1:j1+1, 1:kx)      = xa( 1:2,1:2,1:kx,   ind_pb )
if(ind_t>0)      t (i1:i1+1, j1:j1+1, 1:kx)       = xa( 1:2,1:2,1:kx,   ind_t )
     
!if ( i_t  == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, t)
!if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
!if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
!if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
!if ( i_qc == 0 ) call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx,   1, qc )
!if ( i_qr == 0 ) call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx,   1, qr )
!if ( i_p  == 0 ) call get_variable3d(inputfile, 'P         ', ix, jx, kx,   1, p  )
!if ( i_pb == 0 ) call get_variable3d(inputfile, 'PB        ', ix, jx, kx,   1, pb )

qv (i1:i1+1, j1:j1+1, 1:kx) = qv (i1:i1+1, j1:j1+1, 1:kx) + qc (i1:i1+1, j1:j1+1, 1:kx) + qr (i1:i1+1, j1:j1+1, 1:kx)
ph(i1:i1+1, j1:j1+1, 1:kx+1) = (ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1))/g
p(i1:i1+1, j1:j1+1, 1:kx)  = p(i1:i1+1, j1:j1+1, 1:kx) + pb(i1:i1+1, j1:j1+1, 1:kx)

do j=j1,j1+1
do i=i1,i1+1
  z(i, j, 1:kx) = (ph(i, j, 1:kx)*(znw(2:kx+1)-znu(1:kx))+ph(i, j, 2:kx+1)*(znu(1:kx)-znw(1:kx)))/(znw(2:kx+1)-znw(1:kx))
enddo
enddo
do j=j1,j1+1
do i=i1,i1+1
do k=1,kx
  if( p(i,j,k)/100. >= 1200. ) write(*,'(a,3i4,2f8.0)')'P error at:', i,j,k,p(i,j,k)/100.,pb(i,j,k)/100.
enddo
enddo
enddo
   call compute_seaprs(2, 2, kx, z(i1:i1+1,j1:j1+1,:), t(i1:i1+1,j1:j1+1,:), p(i1:i1+1,j1:j1+1,:), qv(i1:i1+1,j1:j1+1,:), slp, 0)

   xb = dym*(dx*slp(2,1) + dxm*slp(1,1)) + dy*(dx*slp(2,2) + dxm*slp(1,2))
   xb = xb*100.

end subroutine xb_to_slp






!=======================================================================================
subroutine xb_to_radiance(x,proj,ix,jx,kx,ni,nj,nk,nv,nm,xlong,xlat,landmask,iob_radmin,iob_radmax,xb_tb,iid,jid,nslabs_i,nslabs_j,use_cloud)

!---------------------
! radiance subroutine calculates brightness temperature for satellite channels
!---------------------

  USE constants
  USE netcdf
  USE mpi_module
  USE CRTM_Module
  use namelist_define
  use obs_define
  use wrf_tools
  use wrf_field_indices_define


  implicit none

  integer, intent(in)                      :: ix, jx, kx
  integer, intent(in)                      :: nslabs_i, nslabs_j ! Number of slabs in i and j directions
  integer, intent(in)                      :: ni, nj, nk, nv, nm
  integer, intent(out)                     :: iob_radmin,iob_radmax
  real, dimension(ni,nj,nk,nv,nm), intent(in) :: x
  integer, intent(in)                      :: iid, jid
  type(proj_info), intent(in)              :: proj                   ! 1st guestmap info
  logical, intent(in)                      :: use_cloud
  real, dimension(obs%num,nm), intent(out) :: xb_tb
  real, dimension(ix, jx ), intent(in)     :: xlong
  real, dimension(ix, jx ), intent(in)     :: xlat
  real, dimension(ix, jx ), intent(in)     :: landmask
  integer                                  :: iob,irad, comm_size
  real                                     :: dx,dxm,dy,dym
  integer                                  :: istart, iend, jstart, jend
  integer                                  :: slab_ii, slab_jj, obs_ii, obs_jj

  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'ctrm'
!  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/'
!  CHARACTER(*), PARAMETER :: CRTM_code ='/work/03154/tg824524/tools/EnKF_crtm/code/CRTM/crtm_wrf/'
  REAL, PARAMETER :: P1000MB=100000.D0
  REAL, PARAMETER :: R_D=287.D0
  REAL, PARAMETER :: Cpd=7.D0*R_D/2.D0
  REAL, PARAMETER :: Re=6378000.0
  !====================
  !setup for geostationary satellites
   REAL, PARAMETER :: sat_h=35786000.0
   REAL :: sat_lon
  !====================
!  INTEGER, intent(in) :: ix = ix  !total number of the x-grid
!  INTEGER, parameter, intent(in) :: jx = jx  !total number of the y-grid
!  INTEGER, parameter, intent(in) :: kx = kx        !level range
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 1 
!  INTEGER, PARAMETER :: N_LAYERS    = kx
  INTEGER, PARAMETER :: N_ABSORBERS = 2 
!  INTEGER, PARAMETER :: N_CLOUDS    = kx*5
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  INTEGER, PARAMETER :: N_SENSORS = 1
  REAL(fp) :: ZENITH_ANGLE, SCAN_ANGLE, sat_dis

  ! Variables
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  CHARACTER(256) :: FILE_NAME
  CHARACTER(256) :: obstype
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: l, m, irec, yend, ystart, nyi
  integer :: ncid,ncrcode
  character(LEN=16) :: var_name
  character(LEN=3)  :: file_ens
  character(LEN=13) :: tb_fname
  integer :: tt,v,z,n,reci,ens,n_ec,num_radgrid
  INTEGER :: ncl,icl,k1,k2
  integer :: lat_radiance(obs%num)  ! latitude
  integer :: lon_radiance(obs%num) ! longitude
  real :: lat(ix,jx)   ! in radian
  real :: lon(ix,jx)   ! in radian
  real :: pres(ni,nj,kx,nm)
  real :: tk(ni,nj,kx,nm)
  real :: delz(kx)
  real :: tv(ni,nj,kx,nm)
  real :: air_density(ni,nj,kx,nm)
  real, allocatable, dimension(:,:,:,:) :: Tbsend, Tb

  ! Variables regarding sensors detected from radiance_so file
  character(256), dimension(10) :: sensor_list 
  integer                       :: num_sensors, flag_sensor, isen
  integer                       :: n_Channels, iob2


  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_Options_type)                 :: Options(N_PROFILES)
  ! ============================================================================



  xb_tb = 0.
  ! ============================================================================
  ! 1.1 **** make a loop to get the number of satellite-radiance-iob ****
  !
  num_radgrid = 0
  num_sensors = 0
  check_cycle:do iob=1,obs%num
    obstype = obs%type(iob)
    if ( obstype(1:8) == 'Radiance' ) then

      ! Special handling for first radiance observation
      if(num_radgrid == 0) then
        ! Dealing with observation locations
        num_radgrid = num_radgrid + 1
        lon_radiance(num_radgrid) = int(obs%position(iob,1))
        lat_radiance(num_radgrid) = int(obs%position(iob,2))
        iob_radmin = iob
        iob_radmax = iob
        ! Dealing with sensor name
        num_sensors = num_sensors + 1
        sensor_list(num_sensors) = trim(adjustl(obs%sat(iob)))

      ! Non-first radiance obs 
      else
        iob_radmax = iob

        !! Iterative loop to check for double-dipping
        !! Skip to next obs if obs location is wrong.
        !do irad = 1,num_radgrid
        !  if((lon_radiance(irad)==int(obs%position(iob,1))).and.(lat_radiance(irad)==int(obs%position(iob,2))))cycle check_cycle
        !enddo

        ! Loading observation location and incrementing obs number
        num_radgrid = num_radgrid + 1
        lon_radiance(num_radgrid) = int(obs%position(iob,1))
        lat_radiance(num_radgrid) = int(obs%position(iob,2))

        ! Checking if sensor is new
        flag_sensor=1
        do irad = 1, num_sensors
          ! If sensor found, break loop
          if ( sensor_list(irad) == trim(adjustl(obs%sat(iob))) ) then
            flag_sensor=0
            exit
          endif
        enddo
        ! If sensor is new
        if (flag_sensor == 1) then
          num_sensors = num_sensors + 1
          sensor_list(num_sensors) = trim(adjustl(obs%sat(iob)))
        endif

      endif ! End of if statement for num_radgrid ?= 0
    endif ! End of if radiance statement.
  enddo check_cycle


  ! ================================================================
  ! 2. **** Generate useful variables ***
  ! ================================================================
  lat = xlat/180.0*3.14159
  lon = xlong/180.0*3.14159
  pres=0
  tk=0

  pres(:,:,:,:) = x(:,:,1:kx,ind_p,:) + x(:,:,1:kx,ind_pb,:)
  tk = (x(:,:,1:kx,ind_t,:) + 300.0) * ( (pres / P1000MB) ** (R_D/Cpd) )
  tv = tk * (1 + 0.61*x(:,:,:,ind_qvapor,:)/(1+x(:,:,:,ind_qvapor,:)))
  air_density = pres / (tv * 287.2)

  ! Subdomain's start and end indices 
  ! (in terms of full domain indices)
  istart=iid*ni+1
  iend=(iid+1)*ni
  jstart=jid*nj+1
  jend=(jid+1)*nj
  if(iid==(nslabs_i-1)) iend=ix+1
  if(jid==(nslabs_j-1)) jend=jx+1


  ! =================================================================
  ! 3. Performing CRTM for the various available satellites
  ! =================================================================
  sensor_loop: do isen= 1, num_sensors 
    
    ! ===============================================================
    ! 3.1 Sensor ID and related properties
    ! ===============================================================

    Sensor_ID = trim(sensor_list(isen))

    ! Satellite longitude
    if (Sensor_Id == 'ahi_h8' ) then
       sat_lon = 140.7/180.0*3.14159
    else if (Sensor_Id == 'abi_gr') then
       sat_lon = -89.5/180.0*3.14159
    else if (Sensor_Id == 'abi_g16') then
       sat_lon = -75.2/180.0*3.14159
    else if (Sensor_Id == 'imgr_g13') then
       sat_lon = -60.0/180.0*3.14159
    else if (Sensor_Id == 'mviri_m07') then
       sat_lon = 57.0/180.0*3.14159
    else if ( Sensor_Id == 'seviri_m08' ) then
       sat_lon = 41.5/180.0*3.14159
    endif

    ! Note that the channel is specified after the initialization of
    ! CRTM

    ! ===============================================================
    ! 3.2 Initialize CRTM!
    !
    ! Initialization is done using the sensor predefined in Sensor_ID
    !
    ! IMPORTANT: 
    ! The directory that enkf.mpi is being used in 
    ! (e.g., run/201110160000/enkf/d01) must contain a symbolic link
    ! to the directory containing CRTM coefficients. This symbolic
    ! link must be called "coefficients".
    ! ===============================================================

    ! 3.2.1 -- Initialize
    CALL CRTM_Version( Version )
    Error_Status = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hencethe (/../)
                              ChannelInfo  , &  ! Output
                              IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                              IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                              File_Path='coefficients/', &
                              Quiet=.true.)
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error initializing CRTM'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF

    ! 3.2.2 -- Dealing with channel specifics
    ! Himawari-8 AHI or GOES ABI
    if (Sensor_Id == 'abi_gr' .or. Sensor_Id == 'ahi_h8' &
        .or. Sensor_Id == 'abi_g16' ) then
      Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), Channel_Subset =(/8,9,10/) )
      IF ( Error_Status /= SUCCESS ) THEN
        Message = 'Error initializing CRTM'
        CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
        STOP
      END IF
    endif
    ! Meteosat-8 SEVIRI
    if (Sensor_Id == 'seviri_m08' ) then
      Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), Channel_Subset =(/5/) )
      IF ( Error_Status /= SUCCESS ) THEN
        Message = 'Error initializing CRTM'
        CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
        STOP
      END IF
    endif
    ! Meteosat-7 MVIRI
    if (Sensor_Id == 'mviri_m07' ) then
      Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), Channel_Subset =(/2/) )!,3/) )
      IF ( Error_Status /= SUCCESS ) THEN
        Message = 'Error initializing CRTM for Met7'
        CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
        STOP
      END IF
    endif

    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

    ! Initialize container
    allocate(Tb(ni,nj,n_Channels,nm))
    Tb = 0.

    ! 3.2.3 -- Allocating stucture array
    ! -----------------------
    ! Note that only those structure arrays with a channel
    ! dimension are allocated here because we've parameterized
    ! the number of profiles in the N_PROFILES parameter.
    !
    ! Users can make the 
    ! then the INPUT arrays (Atm, Sfc) will also have to be allocated.
    ALLOCATE( RTSolution( n_Channels, N_PROFILES ), STAT=Allocate_Status )
    IF ( Allocate_Status /= 0 ) THEN
      Message = 'Error allocating structure arrays'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF
  
    ! 3.2.4 -- Allocate the STRUCTURES
    ! ---------------------------
    ! The input FORWARD structure
    CALL CRTM_Atmosphere_Create( Atm, kx, N_ABSORBERS, kx*5, N_AEROSOLS)
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
      Message = 'Error allocating CRTM Atmosphere structures'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF


    ! ===============================================================
    ! 3.3 Run CRTM pointwise
    ! ===============================================================
    ! PARALLELIZATION METHOD
    ! ----------------------
    ! Each thread contains a unique combination of slab and ensemble
    ! members. Ie, each thread has to process a unique combination of
    ! obs sites and members
    !
    ! Thus, the only thing any thread has to do is to run CRTM on
    ! the obs sites within said thread's slab and members.
    ! 
    ! ===============================================================



    ! 3.3.2 --- Run pointwise CRTM in parallel
    obs_loop: do iob = 1 , num_radgrid 

      obs_ii=int(lon_radiance(iob))
      obs_jj=int(lat_radiance(iob))

      ! Skip obs if obs not in slab
      if (  ( ( istart > obs_ii ) .or. ( iend < obs_ii ) ) &
          .or. ( ( jstart > obs_jj ) .or. ( jend < obs_jj ) ) ) then
        cycle obs_loop
      endif

      ! Compute subdomain-relative coordinates
      slab_ii = obs_ii - istart + 1
      slab_jj = obs_jj - jstart + 1
     
      iob2 = iob_radmin + iob -1
      ! Skip over obs from different sensor
      if( trim(adjustl(obs%sat(iob2))) /= trim(Sensor_ID ) ) then
        cycle obs_loop
      endif

      ! 3.3.3 -- Converting WRF data for CRTM structure
    
      ! satellite information    
      sat_dis=sqrt(Re**2.0+(Re+sat_h)**2.0-2.0*Re*(Re+sat_h)*cos(lon(obs_ii,obs_jj)-sat_lon)*cos(lat(obs_ii,obs_jj)))
      SCAN_ANGLE=180.0/3.14159*asin(Re/sat_dis*sqrt(1-(cos(lon(obs_ii,obs_jj)-sat_lon)*cos(lat(obs_ii,obs_jj)))**2))
      ZENITH_ANGLE=SCAN_ANGLE+180.0/3.14159*acos(cos(lon(obs_ii,obs_jj)-sat_lon)*cos(lat(obs_ii,obs_jj)))

     
      ! Iterating over the individual members within thread
      members_loop: do n = 1, nm
       
        ! load WRF data into CRTM structures
        !--- calcurating height interval btwn levels (delz)
        do z=1,kx
          if(z.eq.1) then
            delz(z) = ( x(slab_ii, slab_jj, z+1, ind_ph, n) &
                        + x(slab_ii, slab_jj, z+1, ind_phb, n) &
                      )/9.806 &
                      - x( slab_ii, slab_jj, 1, ind_hgt, n)
          else
            delz(z) = (  ( x(slab_ii, slab_jj, z+1, ind_ph, n) &
                           + x(slab_ii, slab_jj, z+1, ind_phb, n) ) &
                       - ( x(slab_ii, slab_jj, z, ind_ph, n) &
                           + x(slab_ii, slab_jj, z, ind_phb, n) ) &
                      )/9.806
          endif
        enddo
        if (delz(1) <= 0.) delz(1) = delz(2)



        !---Atmospheric Profile
        atm(1)%Climatology         = TROPICAL
        atm(1)%Absorber_Id(1:2)    = (/ H2O_ID, O3_ID /)
        atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS,VOLUME_MIXING_RATIO_UNITS /)

        ! Fill in pressure levels in terms of hPa
        atm(1)%Level_Pressure(0) &
        = ( pres(slab_ii,slab_jj,kx,n)*3.0/2.0 &
            - pres(slab_ii,slab_jj,kx-1,n)/2.0 )/100.0

        do z=kx,1,-1
          if(z.eq.1) then
            atm(1)%Level_Pressure(kx-z+1) &
            = max( x(slab_ii, slab_jj, 1, ind_psfc, n),   &
                   pres( slab_ii, slab_jj, 1, n )*3./2.       &
                  -pres( slab_ii, slab_jj, 2, n )/2. )/100.
          !x(slab_ii, slab_jj, 1, ind_psfc, n) 
          else
            atm(1)%Level_Pressure(kx-z+1) &
            = (  ( pres(slab_ii, slab_jj,z-1,n) &
                   + pres(slab_ii,slab_jj,z,n) )/2.0 &
              )/100.0  ! convert from Pa to hPA
          endif
          atm(1)%Pressure(kx-z+1)    = pres(slab_ii, slab_jj, z,n)/100.
          atm(1)%Temperature(kx-z+1) =   tk(slab_ii, slab_jj, z,n)
          ! Dealing with negative water vapor
          if (x(slab_ii, slab_jj, z, ind_qvapor, n) > 0) then
            atm(1)%Absorber(kx-z+1,1)  &
             = x(slab_ii, slab_jj, z, ind_qvapor, n) *1000.
          else
            atm(1)%Absorber(kx-z+1,1)=1e-8*1e3
          endif
         !qvapor(x,y,z)*1000.0
        enddo
        atm(1)%Absorber(:,2) = 5.0E-02 
        !if (my_proc_id == 0) write(*,*) atm(1)%Pressure


        !---Cloud Profile
        do z=1,kx*5
         atm(1)%Cloud(z)%Type = 0
         atm(1)%Cloud(z)%Effective_Radius = 0.0
         atm(1)%Cloud(z)%Water_Content = 0.0
        enddo

        ncl = 0
        icl = 0
        !--calculating # of clouds (cloud and rain)
        do z=kx,1,-1
          if( x(slab_ii,slab_jj,z,ind_qcloud,n).gt.0.0 ) & 
            ncl = ncl + 1
          if( x(slab_ii,slab_jj,z,ind_qrain ,n).gt.0.0 ) & 
            ncl = ncl + 1
          if( x(slab_ii,slab_jj,z,ind_qice  ,n).gt.0.0 ) & 
            ncl = ncl + 1
          if( x(slab_ii,slab_jj,z,ind_qsnow ,n).gt.0.0 ) & 
            ncl = ncl + 1
          if( x(slab_ii,slab_jj,z,ind_qgraup,n).gt.0.0 ) & 
            ncl = ncl + 1
        enddo


        !--Data for cloud
        atm(1)%n_Clouds         = ncl
        IF ( atm(1)%n_Clouds > 0 ) THEN
        do z=kx,1,-1
          if(x(slab_ii,slab_jj,z,ind_qcloud,n).gt.0.0) then
            icl = icl + 1
            k1 = kx-z+1
            k2 = kx-z+1
            atm(1)%Cloud(icl)%Type = WATER_CLOUD
            atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 16.8_fp
            atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
              x( slab_ii, slab_jj, z, ind_qcloud, n )   &
              * air_density( slab_ii, slab_jj, z, n ) * delz(z)
          endif
        enddo
        do z=kx,1,-1
          if(x(slab_ii,slab_jj,z,ind_qrain ,n).gt.0.0) then
            icl = icl + 1
            k1 = kx-z+1
            k2 = kx-z+1
            atm(1)%Cloud(icl)%Type = RAIN_CLOUD
            atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1000.0_fp
            atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
              x( slab_ii, slab_jj, z, ind_qrain, n )   &
              * air_density( slab_ii, slab_jj, z, n ) * delz(z)
          endif
        enddo
        do z=kx,1,-1
          if(x(slab_ii,slab_jj,z,ind_qice,n).gt.0.0) then
            icl = icl + 1
            k1 = kx-z+1
            k2 = kx-z+1
            atm(1)%Cloud(icl)%Type = ICE_CLOUD
            atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 25.0_fp
            atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
              x( slab_ii, slab_jj, z, ind_qice, n )   &
              * air_density( slab_ii, slab_jj, z, n ) * delz(z)
          endif
        enddo
        do z=kx,1,-1
          if(x(slab_ii,slab_jj,z,ind_qsnow,n).gt.0.0) then
            icl = icl + 1
            k1 = kx-z+1
            k2 = kx-z+1
            atm(1)%Cloud(icl)%Type = SNOW_CLOUD
            atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 750.0_fp
            atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
              x( slab_ii, slab_jj, z, ind_qsnow, n )   &
              * air_density( slab_ii, slab_jj, z, n ) * delz(z)
          endif
        enddo
        do z=kx,1,-1
          if(x(slab_ii,slab_jj,z,ind_qgraup,n).gt.0.0) then
            icl = icl + 1
            k1 = kx-z+1
            k2 = kx-z+1
            atm(1)%Cloud(icl)%Type = GRAUPEL_CLOUD
            atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1500.0_fp
            atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
              x( slab_ii, slab_jj, z, ind_qgraup, n )   &
              * air_density( slab_ii, slab_jj, z, n ) * delz(z)
          endif
        enddo
        ENDIF

   
        !---Surface data
        if(landmask(obs_ii, obs_jj).eq.1.0) then
         sfc(1)%Water_Coverage = 0.0_fp
         sfc(1)%Land_Coverage = 1.0_fp
         sfc(1)%Land_Temperature = x(slab_ii, slab_jj, 1, ind_tsk, n) 
         sfc(1)%Soil_Temperature = x(slab_ii, slab_jj, 1, ind_tsk, n) 

        else
         sfc(1)%Water_Coverage = 1.0_fp
         sfc(1)%Land_Coverage = 0.0_fp
         sfc(1)%Water_Type = 1  ! Sea water
         sfc(1)%Water_Temperature= x(slab_ii, slab_jj, 1, ind_tsk, n)

        endif


    
        ! 3.3.4 GeometryInfo input
        ! ----------------------
        ! All profiles are given the same value
        !  The Sensor_Scan_Angle is optional.
        CALL CRTM_Geometry_SetValue( Geometry, &
                                     Sensor_Zenith_Angle = ZENITH_ANGLE, &
                                     Sensor_Scan_Angle   = SCAN_ANGLE )
    
    
        ! 3.3.5 Run CRTM  
        ! --------------------------------------------
        Options%RT_Algorithm_ID = RT_SOI
 
        Error_Status = CRTM_Forward( Atm        , &
                                     Sfc        , &
                                     Geometry   , &
                                     ChannelInfo, &
                                     RTSolution , &
                                     Options = Options )
        IF ( Error_Status /= SUCCESS ) THEN
          Message = 'Error in CRTM Forward Model'
          CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
          STOP
        END IF
   
        ! 3.3.6 Collecting outputs
        !
        ! User should read the user guide or the source code of the routine
        ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
        ! select the needed variables for outputs.  These variables are contained
        ! in the structure RTSolution.
    
        !---for file output, edited 2014.9.26
        do l = 1, n_Channels
            Tb(slab_ii,slab_jj,l,n) &
              = real(RTSolution(l,1)%Brightness_Temperature)

            if( Tb(slab_ii,slab_jj,l,n) /= Tb(slab_ii,slab_jj,l,n) .or. &
                Tb(slab_ii,slab_jj,l,n)>HUGE(Tb(slab_ii,slab_jj,l,n)) .or. &
                Tb(slab_ii,slab_jj,l,n) < 100 .or. Tb(slab_ii,slab_jj,l,n) > 400 ) then
              Tb(slab_ii,slab_jj,l,n)=-888888.
            endif
        enddo
      enddo members_loop !-- end of loop over ensemble members
    enddo obs_loop !--- end of iob(x,y)-loop

    ! 3.3.7 --- writing the output
    Tb_out_loop: do iob =1, num_radgrid 

      obs_ii=int(lon_radiance(iob))
      obs_jj=int(lat_radiance(iob))

      ! Skip obs if obs not in slab
      if (  ( ( istart > obs_ii ) .or. ( iend < obs_ii ) ) &
          .or. ( ( jstart > obs_jj ) .or. ( jend < obs_jj ) ) ) then
        cycle Tb_out_loop 
      endif

      ! Compute subdomain-relative coordinates
      slab_ii = obs_ii - istart + 1
      slab_jj = obs_jj - jstart + 1

      iob2 = iob + iob_radmin - 1

      do n = 1, nm
        ! Dealing with the various sensors.
        ! Dealing with GOES-16 ABI or Him8 AHI
        if (Sensor_Id == 'abi_gr' .or. Sensor_Id == 'ahi_h8' &
            .or. Sensor_Id == 'abi_g16' ) then
            if (obs%ch(iob2) .eq. 8)  &
              xb_tb(iob2,n) = Tb(slab_ii,slab_jj,1,n) !6.19um
            if (obs%ch(iob2) .eq. 9)  &
              xb_tb(iob2,n) = Tb(slab_ii,slab_jj,2,n) !6.95um
            if (obs%ch(iob2) .eq. 10)  &
              xb_tb(iob2,n) = Tb(slab_ii,slab_jj,3,n) !7.34um
        elseif ( Sensor_Id == 'seviri_m08' ) then
            if (obs%ch(iob2) .eq. 5) &
              xb_tb(iob2,n) = Tb(slab_ii,slab_jj,1,n) !6.19um
        elseif (Sensor_Id == 'mviri_m07' ) then
            if (obs%ch(iob2) .eq. 2) &
              xb_tb(iob2,n) = Tb(slab_ii,slab_jj,1,n)   !Meteosat7 ch-3 WV absorb band
        else
            do l = 1, n_Channels
               if (obs%ch(iob2) .eq. RTSolution(l,1)%Sensor_Channel) &
                 xb_tb(iob2,n) = Tb(slab_ii, slab_jj, l,n)
            enddo
        endif
      enddo
    enddo Tb_out_loop

    call MPI_BARRIER( comm, ierr)

    ! ============================================================================
    ! 7. **** DESTROY THE CRTM ****
    !
    deallocate(Tb)
    Error_Status = CRTM_Destroy( ChannelInfo )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error destroying CRTM'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF
    
    deallocate( RTSolution )

    ! ============================================================================
  enddo sensor_loop

end subroutine xb_to_radiance








!=======================================================================================

subroutine xb_to_microwave(inputfile,proj,ix,jx,kx,xlong,xlat,landmask,times,iob_radmin,iob_radmax,xb_tb,write_file_override,use_cloud)

!---------------------
! radiance subroutine calculates brightness temperature for satellite channels
!---------------------

  USE constants
  USE netcdf
  USE mpi_module
  USE CRTM_Module
  use namelist_define
  use obs_define
  use wrf_tools
  use wrf_field_indices_define

  implicit none

  integer, intent(in)                      :: ix, jx, kx
  integer, intent(out)                     :: iob_radmin,iob_radmax
  character(len=10), intent(in)            :: inputfile
  type(proj_info), intent(in)              :: proj                   ! 1st guess map info
  logical, intent(in)                      :: use_cloud
  real, dimension(obs%num), intent(out)    :: xb_tb
  real, dimension(ix, jx ), intent(in)     :: xlong
  real, dimension(ix, jx ), intent(in)     :: xlat
  real, dimension(ix, jx ), intent(in)     :: landmask
  character(len=80), intent(in)            :: times
  logical, intent(in)                      :: write_file_override
  integer                                  :: iob,isensor,n_sensors_used,ich,ich2,ich3
  real                                     :: obs_ii, obs_jj, dx,dxm,dy,dym

  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'ctrm'
  INTEGER, PARAMETER :: N_SUPPORTED_SENSOR_CLASSES = 9
  CHARACTER(*), PARAMETER, DIMENSION(N_SUPPORTED_SENSOR_CLASSES) :: SENSOR_CLASS_NAMES = &
      (/'ssmis','ssmi','gmi_gpm_hf','gmi_gpm_lf','saphir_meghat','amsr2_gcom-w1','atms_npp','atms_n20','mhs'/)
  INTEGER, DIMENSION(N_SUPPORTED_SENSOR_CLASSES) :: N_CHANNELS = &
      (/24, 7, 13, 13, 6, 14, 22, 22, 5/)

  REAL, PARAMETER :: P1000MB=100000.D0
  REAL, PARAMETER :: R_D=287.D0
  REAL, PARAMETER :: Cpd=7.D0*R_D/2.D0
  REAL, PARAMETER :: Re=6378000.0
  REAL, PARAMETER :: MY_PI = 4*ATAN(1.)
  REAL, PARAMETER :: RADS_PER_DEGREE = MY_PI / 180.
  !====================
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 1
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  INTEGER, PARAMETER :: N_SENSORS = 1
  INTEGER, PARAMETER :: N_STREAMS = 16
  INTEGER, PARAMETER :: x_coarse = 1, y_coarse = 1  ! how many wrf gridpoints should be averaged together in a CRTM profile

  ! Variables
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_ID
  CHARACTER(256) :: CRTM_OUT_DIR
  CHARACTER(256) :: FILE_NAME
  CHARACTER(4)   :: YEAR_STR
  CHARACTER(2)   :: MONTH_STR, DAY_STR
  CHARACTER(3)   :: DOMAIN
  CHARACTER(256) :: obstype
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: sum_n_channels        ! total number of channels across all included sensors
  INTEGER :: n_Channels_all        ! total number of channels with observations across all sensors
  INTEGER :: n_Channels_crtm       ! total number of channels with observations across all sensors that need CRTM calculations
  INTEGER :: l, m, irec, grand_count, yend, ystart, xend, xstart
  integer :: ncid,ncrcode,fid
  character(LEN=16) :: var_name
  character(LEN=3)  :: file_ens
  integer :: x,y,i,j,k,tt,v,z,n,reci,ens,n_ec,num_radgrid
  INTEGER :: ncl,icl,k1,k2
  real :: sigma, search_radius
  real :: cputime
  integer, allocatable, dimension(:) :: lat_radiance_min, lat_radiance_max ! latitude
  integer, allocatable, dimension(:) :: lon_radiance_min, lon_radiance_max ! longitude
  integer :: search_lat_min, search_lat_max
  integer :: search_lon_min, search_lon_max
  real :: beam_conv_simple, beam_conv_gaussian_simple
  real(fp) :: my_scan_angle, my_zenith_angle, my_azimuth_angle, my_azimuth_sine, my_azimuth_cosine, my_azimuth_angle_est
  integer :: my_angle_sum_count, scan_angle_search_radius
  INTEGER :: YEAR, MONTH, DAY

  logical, DIMENSION(:), allocatable :: is_channel_used, is_channel_crtm   ! logical array, length sum_n_channels
  integer, DIMENSION(:), allocatable :: channel_numbers_crtm, sensor_index_crtm, channel_search   ! length n_Channels_crtm
  integer, DIMENSION(:,:), allocatable :: channel_numbers_all_eachSensor   ! size n_sensors_used X n_Channels_all (using a total of only n_Channels_crtm elements)
  integer, allocatable, dimension(:) :: n_channels_all_eachSensor  ! length number of sensors
  INTEGER :: n_Channels_all_mySensor   ! total number of channels with observations for the sensor of interest
  integer, allocatable, dimension(:) :: channel_numbers_all_mySensor  ! length n_channels_crtm_mySensor
  INTEGER :: mySensor_Tb_offset  ! offset of channel numbers in the Tb array for sensor of interest
  integer, DIMENSION(:,:), allocatable :: channel_numbers_crtm_eachSensor   ! size n_sensors_used X n_Channels_crtm (using a total of only n_Channels_crtm elements)
  integer, allocatable, dimension(:) :: n_channels_crtm_eachSensor  ! length number of sensors
  integer :: n_Channels_crtm_mySensor  ! total number of channels with observations for the sensor of interest that need CRTM calculations
  integer, allocatable, dimension(:) :: channel_numbers_crtm_mySensor  ! length n_channels_crtm_mySensor
  integer, allocatable, dimension(:) :: sensors

  logical :: crtm_out_files
  integer, allocatable, dimension(:) :: numx_crtm, numy_crtm  ! length number of sensors
  real :: lat(ix,jx)   ! in radian
  real :: lon(ix,jx)   ! in radian
  real :: p(ix,jx,kx)
  real :: pb(ix,jx,kx)
  real :: pres(ix,jx,kx)
  real :: ph(ix,jx,kx+1)
  real :: phb(ix,jx,kx+1)
  real :: delz(kx)
  real :: t(ix,jx,kx)
  real :: tk(ix,jx,kx)
  real :: qvapor(ix,jx,kx)
  real :: qcloud(ix,jx,kx)
  real :: qrain(ix,jx,kx)
  real :: qice(ix,jx,kx)
  real :: qsnow(ix,jx,kx)
  real :: qgraup(ix,jx,kx)
  real :: nrain(ix,jx,kx)
  real :: psfc(ix,jx)
  real :: hgt(ix,jx)
  real :: tsk(ix,jx)
  real :: xland(ix,jx)
  real :: u10(ix,jx)
  real :: v10(ix,jx)
  real :: windspeed(ix,jx)
  real :: westwind(ix,jx)
  real :: winddir(ix,jx)
  real, allocatable, dimension(:,:,:) :: Tbsend, Tb_crtm, Tb  ! size # sensors X spatial grid
  real, allocatable, dimension(:,:,:) :: scan_angle, zenith_angle, azimuth_angle, azimuth_sine, azimuth_cosine  ! size # sensors X spatial grid
  real :: azimuth_angle_theta
  integer, allocatable, dimension(:,:,:) :: angle_sum_count  ! size # sensors X spatial grid
  logical :: is_sensor_used(N_SUPPORTED_SENSOR_CLASSES)
  integer, allocatable, dimension(:) :: sensors_used, sensor_ch_offset  ! length number of sensors
  CHARACTER(256), allocatable, dimension(:) :: Sensor_IDs_used
  CHARACTER(256) :: CRTM_Init_Sensor_ID
  integer, allocatable, dimension(:) :: sensor_indices   ! length number of obs
  integer :: my_sensor_index
  integer :: max_n_channels_all_sensor
  integer :: scan_angle_search_number
  integer :: scan_angle_search_lonmin, scan_angle_search_lonmax
  integer :: scan_angle_search_latmin, scan_angle_search_latmax

  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_Options_type)                 :: Options(N_PROFILES)

  integer :: mp_scheme = 0  ! wsm6, gfdlfv3, thompson08
  logical :: l_use_default_re = .false.
  integer, parameter :: N_CLOUDS=5

  ! ============================================================================

  ! ============================================================================
  ! 1.5. **** make a loop to get the number of satellite-microwave-iob ****
  !      **** and lat-lon bounds for CRTM calculations                 ****

  if(my_proc_id==0) then
    call cpu_time(cputime)
    write(*,*) 'entered xb.f at ', cputime
  endif


  call open_file(inputfile, nf_nowrite, fid)
  ncrcode = nf_get_att_real(fid, nf_global,'DX', dx)
  ncrcode = nf_get_att_real(fid, nf_global,'DY', dy)
  ! convert to km
  dx = dx / 1000.
  dy = dy / 1000.
  call close_file(fid)


  ! determine the number of unique sensors of observation data in the file
  ! determine the total number of microwave observations
  is_sensor_used = .FALSE.
  num_radgrid = 0
  check_cycle1:do iob=1,obs%num
    obstype = obs%type(iob)
    if ( obstype(1:9) == 'Microwave' ) then
      if(num_radgrid == 0) then
        iob_radmin = iob
      endif
      iob_radmax = iob
      num_radgrid = num_radgrid + 1
      ! Get sensor id from observation object
      Sensor_ID = trim(adjustl(obs%sat(iob)))
      do isensor = 1,N_SUPPORTED_SENSOR_CLASSES
        if(INDEX(trim(Sensor_ID), trim(SENSOR_CLASS_NAMES(isensor)))>0) exit
      enddo
      is_sensor_used(isensor) = .TRUE.
    endif
  enddo check_cycle1
  n_sensors_used = count(is_sensor_used)
  if (my_proc_id==0) write(*,*) 'n_sensors_used:', n_sensors_used

  ! determine number of channels across all sensors with observations
  allocate(sensors_used(n_sensors_used), STAT = Allocate_Status)
  allocate(sensor_ch_offset(n_sensors_used), STAT = Allocate_Status)
  sensor_ch_offset = 0
  sum_n_channels = 0
  max_n_channels_all_sensor = 0
  i = 1
  do isensor = 1,N_SUPPORTED_SENSOR_CLASSES
    if (is_sensor_used(isensor)) then
      sum_n_channels = sum_n_channels + N_CHANNELS(isensor)
      max_n_channels_all_sensor = max(max_n_channels_all_sensor, N_CHANNELS(isensor))
      sensors_used(i) = isensor   ! will result in sensor indices of increasing value from 1 to number of sensors
      if (i < n_sensors_used) then
        sensor_ch_offset(i+1) = sum(sensor_ch_offset(1:i)) + N_CHANNELS(isensor)
        i = i + 1
      else
        exit
      endif
    endif
  enddo
  if(my_proc_id==0) write(*,*) 'sensors_used:',sensors_used

  ! allocate arrays. 
  ! the channels are arranged by the order of sesnors in SENSOR_CLASS_NAMES

  allocate(lon_radiance_min(n_sensors_used), STAT = Allocate_Status)
  allocate(lon_radiance_max(n_sensors_used), STAT = Allocate_Status)
  allocate(lat_radiance_min(n_sensors_used), STAT = Allocate_Status)
  allocate(lat_radiance_max(n_sensors_used), STAT = Allocate_Status)
  lon_radiance_min = ix
  lon_radiance_max = 0
  lat_radiance_min = jx
  lat_radiance_max = 0

  allocate(scan_angle(n_sensors_used,ix,jx), STAT = Allocate_Status)
  allocate(zenith_angle(n_sensors_used,ix,jx), STAT = Allocate_Status)
  allocate(azimuth_angle(n_sensors_used,ix,jx), STAT = Allocate_Status)
  allocate(azimuth_sine(n_sensors_used,ix,jx), STAT = Allocate_Status)
  allocate(azimuth_cosine(n_sensors_used,ix,jx), STAT = Allocate_Status)
  allocate(angle_sum_count(n_sensors_used,ix,jx), STAT = Allocate_Status)
  scan_angle = 0.d0
  zenith_angle = 0.d0
  azimuth_angle = 0.d0
  azimuth_sine = 0.d0
  azimuth_cosine = 0.d0
  angle_sum_count = 0

  allocate(is_channel_used(sum_n_channels), STAT = Allocate_Status)
  allocate(is_channel_crtm(sum_n_channels), STAT = Allocate_Status)
  is_channel_used = .FALSE.
  is_channel_crtm = .FALSE.

  allocate(sensor_indices(num_radgrid), STAT = Allocate_Status)

  allocate(Sensor_IDs_used(n_sensors_used), STAT = Allocate_Status)

  ! determine locations of observations for each sensor
  ! get Sensor_ID for each sensor used
  check_cycle2:do i=1,num_radgrid
    iob = i + iob_radmin - 1

    ! Get sensor id index for each observation, save them for easy reference later
    Sensor_ID = trim(adjustl(obs%sat(iob)))
    do isensor = 1,N_SUPPORTED_SENSOR_CLASSES
      if(INDEX(trim(Sensor_ID), trim(SENSOR_CLASS_NAMES(isensor)))>0) exit
    enddo
    do my_sensor_index = 1,n_sensors_used
      if(sensors_used(my_sensor_index) .EQ. isensor) exit
    enddo
    sensor_indices(i) = my_sensor_index
    Sensor_IDs_used(my_sensor_index) = trim(Sensor_ID)

    ! identify that this channel will be used
    is_channel_used(obs%ch(iob) + sensor_ch_offset(my_sensor_index)) = .TRUE.

    ! determine satellite beam convolution area for the observation
    !if (my_proc_id == 0) write(*,*) 'about to calculate sigma:\n', iob, obs%efov_aScan(iob), obs%efov_cScan(iob)
    sigma = (0.5 * ( ( ( ( obs%efov_aScan(iob) + obs%efov_cScan(iob) ) /2 ) / 1.18) ) )
    search_radius = 1 + (2 * sigma)

    search_lon_min = max(1, nint(obs%position(iob,1)) - ceiling(1+(search_radius / dx)))
    search_lon_max = min(ix,nint(obs%position(iob,1)) + ceiling(1+(search_radius / dx)))
    search_lat_min = max(1, nint(obs%position(iob,2)) - ceiling(1+(search_radius / dy)))
    search_lat_max = min(jx,nint(obs%position(iob,2)) + ceiling(1+(search_radius / dy)))

    ! include specifications for this observation in determining the optimal
    ! scan and zenith angles for CRTM simulations
    scan_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    scan_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + abs(obs%scan_angle(iob))

    zenith_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    zenith_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + abs(obs%zenith_angle(iob))

    azimuth_angle_theta = 90 - obs%azimuth_angle(iob)
    if (azimuth_angle_theta < 0) azimuth_angle_theta = azimuth_angle_theta + 360
    azimuth_angle_theta = azimuth_angle_theta * RADS_PER_DEGREE

    azimuth_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    azimuth_angle(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + azimuth_angle_theta

    azimuth_sine(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    azimuth_sine(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + sin(azimuth_angle_theta)

    azimuth_cosine(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    azimuth_cosine(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + cos(azimuth_angle_theta)

    angle_sum_count(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) = &
    angle_sum_count(my_sensor_index,search_lon_min:search_lon_max,search_lat_min:search_lat_max) + 1

    if (.NOT. write_file_override) then
      lon_radiance_min(my_sensor_index) = min(lon_radiance_min(my_sensor_index), nint(obs%position(iob,1)) &
                                        - ceiling(1+(search_radius / dx)))
      lon_radiance_max(my_sensor_index) = max(lon_radiance_max(my_sensor_index), nint(obs%position(iob,1)) &
                                        + ceiling(1+(search_radius / dx)))
      lat_radiance_min(my_sensor_index) = min(lat_radiance_min(my_sensor_index), nint(obs%position(iob,2)) &
                                        - ceiling(1+(search_radius / dy)))
      lat_radiance_max(my_sensor_index) = max(lat_radiance_max(my_sensor_index), nint(obs%position(iob,2)) &
                                        + ceiling(1+(search_radius / dy)))
    endif

  enddo check_cycle2
  if (my_proc_id==0 .AND. inputfile == 'fort.80011') write(*,*) "Sensor_IDs_used: ", Sensor_IDs_used

  if (write_file_override) then
    lon_radiance_min = 1
    lon_radiance_max = ix
    lat_radiance_min = 1
    lat_radiance_max = jx
  else
    do isensor = 1,n_sensors_used
      lon_radiance_min(isensor) = max(1,lon_radiance_min(isensor))
      lon_radiance_max(isensor) = min(ix,lon_radiance_max(isensor))
      lat_radiance_min(isensor) = max(1,lat_radiance_min(isensor))
      lat_radiance_max(isensor) = min(jx,lat_radiance_max(isensor))
    enddo
  endif

  allocate(numx_crtm(n_sensors_used), STAT = Allocate_Status)
  allocate(numy_crtm(n_sensors_used), STAT = Allocate_Status)
  numx_crtm = lon_radiance_max - lon_radiance_min + 1
  if(my_proc_id==0)  write(*,*) "numx_crtm: ", numx_crtm
  numy_crtm = lat_radiance_max - lat_radiance_min + 1
  if(my_proc_id==0)  write(*,*) "numy_crtm: ", numy_crtm

  ! ============================================================================
  ! 2. **** INITIALIZE THE CRTM ****
  !
  ! 2a. Determine which if any channels that 
  !     brightness temperatures have already
  !     been calculated, and load those
  ! ------------------------------------------

  n_Channels_all = count(is_channel_used)
  ALLOCATE(n_channels_all_eachSensor(n_sensors_used), STAT = Allocate_status)
  ALLOCATE(channel_numbers_all_eachSensor(n_sensors_used,max_n_channels_all_sensor), STAT = Allocate_status)
  n_channels_all_eachSensor = 0

  ! BTs of all cahnnels from all sensors is sotred in a single array.
  ! allocate array for entire WRF domain, but use only locations relevant
  ! to the observations from each sensor.
  allocate(Tb(ix,jx,n_Channels_all))
  Tb = -8888.

  ! get array of channels used and the sensor they are with, and
  ! load brightness temperature file if it exists. if no file exists,
  ! then identify the channel for CRTM simulation
  CALL GET_ENVIRONMENT_VARIABLE('CRTM_OUT_DIR',CRTM_OUT_DIR)
  ich2 = 0
  crtm_out_files = .FALSE.

  do isensor = 1,n_sensors_used
    ich3 = 0
    do ich = 1,N_CHANNELS(sensors_used(isensor))
      if (is_channel_used(sensor_ch_offset(isensor) + ich)) then
        ich2 = ich2 + 1
        ich3 = ich3 + 1
        channel_numbers_all_eachSensor(isensor,ich3) = ich
        write(FILE_NAME,'(a,a,a,a,i2.2,a,a)') trim(CRTM_OUT_DIR), '/Tb_', trim(SENSOR_CLASS_NAMES(sensors_used(isensor))), '_ch', ich,'_', trim(inputfile)
        INQUIRE(FILE=trim(FILE_NAME),EXIST=crtm_out_files)
        if (crtm_out_files) then ! read file of brightness temperatures
          if(my_proc_id==0)  WRITE(*,*) '  Loading BT data from this file: ', trim(FILE_NAME)
          if (numx_crtm(isensor) .NE. ix .OR. numy_crtm(isensor) .NE. jx) then
            lon_radiance_min(isensor) = 1
            lon_radiance_max(isensor) = ix
            lat_radiance_min(isensor) = 1
            lat_radiance_max(isensor) = jx
            numx_crtm(isensor) = ix
            numy_crtm(isensor) = jx
            if(my_proc_id==0)  write(*,*) "  new numx_crtm(isensor) and numy_crtm(isensor): ", numx_crtm(isensor), numy_crtm(isensor)
          endif

          open(ens,file=trim(FILE_NAME),access='direct',recl=4)
          irec = 0
          do j = 1, numy_crtm(isensor)
          do i = 1, numx_crtm(isensor)
            irec= irec +1
            read( ens, rec=irec) Tb(i,j,ich3)
          enddo
          enddo
          close (ens)
        else
          is_channel_crtm(sensor_ch_offset(isensor) + ich) = .TRUE.
        endif
      endif
    enddo
    n_channels_all_eachSensor(isensor) = ich3
    if (ich2 .eq. n_Channels_all) exit  ! stop when all the channels have been found
  enddo
  deallocate(is_channel_used)

  ! 2b. Determine the total number of channels
  !     for which the CRTM was initialized
  ! ------------------------------------------

  n_Channels_crtm = count(is_channel_crtm)

  ! determine sensor indices and channel numbers for crtm
  ALLOCATE(channel_numbers_crtm(n_Channels_crtm), STAT = Allocate_Status)
  ALLOCATE(sensor_index_crtm(n_Channels_crtm), STAT = Allocate_Status)
  ALLOCATE(channel_numbers_crtm_eachSensor(n_sensors_used,n_Channels_crtm), STAT = Allocate_status)
  ALLOCATE(n_channels_crtm_eachSensor(n_sensors_used), STAT = Allocate_status)
  channel_numbers_crtm_eachSensor = 0
  ich2 = 0
  do isensor = 1,n_sensors_used
    ich3 = 0
    do ich = 1,N_CHANNELS(sensors_used(isensor))
      !if(my_proc_id==0)  write(*,*) "is channel ", ich, " of sensor ",
      !sensors_used(isensor), " for crtm?"
      if (is_channel_crtm(sensor_ch_offset(isensor) + ich)) then
        !if(my_proc_id==0)  write(*,*) "yes"
        ich2 = ich2 + 1
        channel_numbers_crtm(ich2) = ich
        sensor_index_crtm(ich2) = isensor

        ich3 = ich3 + 1
        channel_numbers_crtm_eachSensor(isensor,ich3) = ich
      endif
    enddo
    n_channels_crtm_eachSensor(isensor) = ich3
    if (ich2 .eq. n_Channels_crtm) exit  ! stop when all the channels have been found
  enddo
  deallocate(is_channel_crtm)

  ! ============================================================================
  ! 3. **** GET INPUT DATA ****
  !
  ! For fill the Atm structure array.
  !
  ! 4a1. Loading Atmosphere and Surface input
  ! --------------------------------
  call get_variable3d(inputfile,'P',ix,jx,kx,1,p)
  call get_variable3d(inputfile,'PB',ix,jx,kx,1,pb)
  call get_variable3d(inputfile,'PH',ix,jx,kx+1,1,ph)
  call get_variable3d(inputfile,'PHB',ix,jx,kx+1,1,phb)
  call get_variable3d(inputfile,'T',ix,jx,kx,1,t)
  call get_variable3d(inputfile,'QVAPOR',ix,jx,kx,1,qvapor)
  if(my_proc_id==0)  WRITE(*,*) 'use cloud: ', use_cloud
  if (use_cloud) then
     call get_variable3d(inputfile,'QCLOUD',ix,jx,kx,1,qcloud)
     call get_variable3d(inputfile,'QRAIN',ix,jx,kx,1,qrain)
     call get_variable3d(inputfile,'QICE',ix,jx,kx,1,qice)
     call get_variable3d(inputfile,'QSNOW',ix,jx,kx,1,qsnow)
     call get_variable3d(inputfile,'QGRAUP',ix,jx,kx,1,qgraup)
  else
     qcloud = 0.
     qrain = 0.
     nrain = 0.
     qice = 0.
     qsnow = 0.
     qgraup = 0.
  endif
  call get_variable2d(inputfile,'PSFC',ix,jx,1,psfc)
  call get_variable2d(inputfile,'TSK',ix,jx,1,tsk)
  call get_variable2d(inputfile,'XLAND',ix,jx,1,xland)
  call get_variable2d(inputfile,'HGT',ix,jx,1,hgt)
  call get_variable2d(inputfile,'U10',ix,jx,1,u10)
  call get_variable2d(inputfile,'V10',ix,jx,1,v10)

  lat = xlat/180.0*MY_PI
  lon = xlong/180.0*MY_PI
  pres = P + PB
  tk = (T + 300.0) * ( (pres / P1000MB) ** (R_D/Cpd) )
  where(qvapor.lt.0.0) qvapor=1.0e-8
  where(qcloud.lt.0.0) qcloud=0.0
  where(qice.lt.0.0) qice=0.0
  where(qrain.lt.0.0) qrain=0.0
  where(nrain.lt.0.0) nrain=0.0
  where(qsnow.lt.0.0) qsnow=0.0
  where(qgraup.lt.0.0) qgraup=0.0
  windspeed = sqrt(U10**2 + V10**2)
  where (ISNAN(windspeed)) windspeed = 0.0d0
  westwind = (U10 .lt. 0)
  ! corrected for CRTM version 2.3.0
  ! winddir = -180*(westwind) + ( 90 - ( atan(V10/U10)/RADS_PER_DEGREE) )
  winddir = 180*(1 - westwind) + (90 - ( atan(V10/U10)/RADS_PER_DEGREE) )
  where (ISNAN(winddir)) winddir = 0.0d0

  ! ============================================================================
  ! 4. **** INITIALIZE CRTM AND ALLOCATE STRUCTURE ARRAYS ****
  !
  ! 4a. This initializes the CRTM for the sensors
  !     predefined in the example SENSOR_ID parameter.
  !     NOTE: The coefficient data file path is hard-
  !           wired for this example.
  ! --------------------------------------------------

  CALL CRTM_Version( Version )

  do isensor = 1,n_sensors_used

    ! determine if CRTM gets run for any channels of this sensor
    n_channels_crtm_mySensor = n_channels_crtm_eachSensor(isensor)
    if (n_channels_crtm_mySensor > 0) then

    if(my_proc_id==0) then
      call cpu_time(cputime)
      WRITE(*,*) 'starting CRTM for sensor ', isensor,' at ', cputime
      WRITE(*,*) 'sensor ', isensor,' is ', trim(Sensor_IDs_used(isensor))
    endif

    ! If sensor class is either gpm_gmi_lf or gmi_gpm_hf (low- or high-frequency),
    ! then initialize the CRTM with 'gpm_gmi'. Otherwise, use the actual Sensor_ID
    ! specified by the observation.
    ! 'gpm_gmi_lf' is channels 1-9, while 'gpm_gmi_hg' is channels 10-13. 
    ! These two sets of channels have different scan and zenith angles. If observations
    ! from both of these two sets of channels are assimilated, then treating these sets of
    ! channels as different classes of sensors is necessary to have the calculated average 
    ! scan and zenith angles of observations in the vicinity of a given model grid point 
    ! to be correct for each observation. However, the CRTM still needs to be told 'gmi_gpm'
    ! for either.
    if(INDEX(trim(Sensor_IDs_used(isensor)), 'gmi')>0) then
      CRTM_Init_Sensor_ID = 'gmi_gpm'
    else
      CRTM_Init_Sensor_ID = trim(Sensor_IDs_used(isensor))
    endif

    if(my_proc_id==0)  write(*,*)
    if(my_proc_id==0)  write(*,*) "CRTM ver.",TRIM(Version)
    if (my_proc_id==0 .AND. (inputfile == 'fort.80011' .OR. inputfile == 'fort.80071')) then
      Error_Status = CRTM_Init( (/CRTM_Init_Sensor_ID/), &  ! Input... must be an array, hencethe (/../)
                                ChannelInfo  , &  ! Output
                                IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                                IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                                File_Path='coefficients/', &
                                CloudCoeff_File_rain  = 'WSM6_RainLUT_-109z-1.bin',&
                                CloudCoeff_File_snow  = 'WSM6_SnowLUT_-109z-1.bin',&
                                CloudCoeff_File_graup = 'WSM6_GraupelLUT_-109z-1.bin',&
                                Quiet=.false.)
    else
      Error_Status = CRTM_Init( (/CRTM_Init_Sensor_ID/), &  ! Input... must be an array, hencethe (/../)
                                ChannelInfo  , &  ! Output
                                IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                                IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                                File_Path='coefficients/', &
                                CloudCoeff_File_rain  = 'WSM6_RainLUT_-109z-1.bin',&
                                CloudCoeff_File_snow  = 'WSM6_SnowLUT_-109z-1.bin',&
                                CloudCoeff_File_graup = 'WSM6_GraupelLUT_-109z-1.bin',&
                                Quiet=.true.)
    end if

    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error initializing CRTM'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF

    n_channels_crtm_mySensor = n_channels_crtm_eachSensor(isensor)
    ALLOCATE(channel_numbers_crtm_mySensor(n_channels_crtm_mySensor), STAT=Allocate_Status )
    channel_numbers_crtm_mySensor = channel_numbers_crtm_eachSensor(isensor,1:n_channels_crtm_mySensor)

    Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), Channel_Subset = channel_numbers_crtm_mySensor )

    ! 4b. Allocate the ARRAYS
    ! -----------------------
    ! Note that only those structure arrays with a channel
    ! dimension are allocated here because we've parameterized
    ! the number of profiles in the N_PROFILES parameter.
    !
    ! Users can make the 
    ! then the INPUT arrays (Atm, Sfc) will also have to be allocated.
    ALLOCATE( RTSolution( n_channels_crtm_mySensor, N_PROFILES ), STAT=Allocate_Status )
    IF ( Allocate_Status /= 0 ) THEN
      Message = 'Error allocating structure arrays'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF

    ! 4c. Allocate the STRUCTURES
    ! ---------------------------
    ! The input FORWARD structure
    CALL CRTM_Atmosphere_Create( Atm, kx, N_ABSORBERS, N_CLOUDS, N_AEROSOLS)
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
      Message = 'Error allocating CRTM Atmosphere structures'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF

    allocate(Tbsend(numx_crtm(isensor),numy_crtm(isensor),n_channels_crtm_mySensor))
    allocate(Tb_crtm(numx_crtm(isensor),numy_crtm(isensor),n_channels_crtm_mySensor))
    Tbsend = 0.
    Tb_crtm = 0.

    ! 4a2. Parallerization with grids
    ! --------------------------------
    grand_count = 0
    do j=1, numy_crtm(isensor)
       y = (j-1) + lat_radiance_min(isensor)
    do i=1, numx_crtm(isensor)
       x = (i-1) + lon_radiance_min(isensor)
       grand_count = grand_count + 1
       !if(my_proc_id==0)  write(*,'(a,i4,i4,i6)') 'point ', x, y, grand_count
       if (mod(grand_count,nprocs) .EQ. my_proc_id) then ! calcualte if this processor should

    ! 4a3. Converting WRF data for CRTM structure
    ! --------------------------------
    !--- converting the data to CRTM structure

    !*******************************************************************************
    ! satellite information is contained in parameter arrays
    !*******************************************************************************


    ! 4a. GeometryInfo input
    ! ----------------------
    ! if this location was not designated for use by any observation, then use the average value of
    ! nearby locations. increase search "radius" until averaging over at least 2 locations.
    my_scan_angle = scan_angle(isensor,x,y)
    scan_angle_search_radius = 0
    scan_angle_search_number = 0
    if (my_scan_angle < 2*tiny(my_scan_angle)) then
      do while (scan_angle_search_number < 2)
        scan_angle_search_radius = scan_angle_search_radius + 1
        scan_angle_search_lonmin = max(1, x-scan_angle_search_radius)
        scan_angle_search_lonmax = min(ix,x+scan_angle_search_radius)
        scan_angle_search_latmin = max(1, y-scan_angle_search_radius)
        scan_angle_search_latmax = min(jx,y+scan_angle_search_radius)
        scan_angle_search_number = count(angle_sum_count(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                                 scan_angle_search_latmin:scan_angle_search_latmax) .GT. 1)
      enddo

      my_scan_angle = sum(scan_angle(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                 scan_angle_search_latmin:scan_angle_search_latmax) )
      my_zenith_angle = sum(zenith_angle(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                     scan_angle_search_latmin:scan_angle_search_latmax) )
      my_azimuth_sine = sum(azimuth_sine(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                     scan_angle_search_latmin:scan_angle_search_latmax) )
      my_azimuth_cosine = sum(azimuth_cosine(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                         scan_angle_search_latmin:scan_angle_search_latmax) )
      my_azimuth_angle_est = sum(azimuth_angle(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                           scan_angle_search_latmin:scan_angle_search_latmax) )
      my_angle_sum_count = sum(angle_sum_count(isensor,scan_angle_search_lonmin:scan_angle_search_lonmax, &
                                                           scan_angle_search_latmin:scan_angle_search_latmax) )

    else
      my_zenith_angle = zenith_angle(isensor,x,y)
      my_azimuth_sine = azimuth_sine(isensor,x,y)
      my_azimuth_cosine = azimuth_cosine(isensor,x,y)
      my_azimuth_angle_est = azimuth_angle(isensor,x,y)
      my_angle_sum_count = angle_sum_count(isensor,x,y)
    endif

    my_scan_angle = my_scan_angle / my_angle_sum_count
    my_zenith_angle = my_zenith_angle / my_angle_sum_count

    my_azimuth_angle = atan2(my_azimuth_sine / my_angle_sum_count, my_azimuth_cosine / my_angle_sum_count) / RADS_PER_DEGREE
    my_azimuth_angle = 90 - my_azimuth_angle
    if (my_azimuth_angle < 0) my_azimuth_angle = my_azimuth_angle + 360

    my_azimuth_angle_est = (my_azimuth_angle_est / my_angle_sum_count) / RADS_PER_DEGREE
    my_azimuth_angle_est = 90 - my_azimuth_angle_est
    if (my_azimuth_angle_est < 0) my_azimuth_angle_est = my_azimuth_angle_est + 360

    if(my_proc_id==0 .AND. inputfile == 'fort.80011' .AND. j .EQ. 1 ) write(*,*) "    i", i, ", my_azimuth_angle and my_azimuth_angle_est: ", my_azimuth_angle, my_azimuth_angle_est
    if ( my_scan_angle > 90 ) then
      write(*,*) "x",x," y",y, ", adjusted my_scan_angle: ", my_scan_angle
      write(*,*) "  final scan_angle_search_radius: ", scan_angle_search_radius
      write(*,*) "  my_angle_sum_count: ", my_angle_sum_count
      write(*,*) "  adjusted my_zenith_angle: ", my_zenith_angle
      write(*,*) "  scan_angle_search_lon: ", scan_angle_search_lonmin, scan_angle_search_lonmax
      write(*,*) "  scan_angle_search_lat: ", scan_angle_search_latmin, scan_angle_search_latmax
      write(*,*) "  my_proc_id: ", my_proc_id
      write(*,*)
    endif

    YEAR_STR = times(1:4)
    MONTH_STR = times(6:7)
    DAY_STR = times(9:10)
    read(YEAR_STR,*) YEAR
    read(MONTH_STR,*) MONTH
    read(DAY_STR,*) DAY

    CALL CRTM_Geometry_SetValue( Geometry, &
                                 Sensor_Zenith_Angle  = my_zenith_angle, &
                                 Sensor_Scan_Angle    = my_scan_angle, &
                                 Sensor_Azimuth_Angle = my_azimuth_angle, &
                                 Year                 = YEAR, &
                                 Month                = MONTH, &
                                 Day                  = DAY )
    ! 4b. Converting WRF data for CRTM structure
    ! --------------------------------

    if ( use_slant_path) then
      call Load_CRTM_Structures_MW_slant( Atm, Sfc, x, y, kx, &
               my_zenith_angle, my_azimuth_angle, dx, .FALSE.)
    else
      call Load_CRTM_Structures_MW( Atm, Sfc, x, y, kx, .FALSE.)
    end if

    ! 4c. Use the SOI radiative transfer algorithm
    ! --------------------------------------------
    Options%RT_Algorithm_ID = RT_SOI


    ! 4d. Specify number of streams
    IF (N_STREAMS .GT. 0) THEN
      Options%Use_N_Streams = .TRUE.
      Options%n_Streams = N_STREAMS
    END IF
    ! ============================================================================

    ! ============================================================================
    ! 5. **** CALL THE CRTM FORWARD MODEL ****

    if (Sfc(1)%Water_Coverage < 0) then
      WRITE(*,*)
      WRITE(*,*) 'x: ', x
      WRITE(*,*) 'y: ', y
      WRITE(*,*) 'xland(x,y): ', xland(x,y)
      WRITE(*,*) 'Water coverage: ', Sfc(1)%Water_Coverage
      WRITE(*,*) 'Land coverage: ', Sfc(1)%Land_Coverage
      WRITE(*,*) 'Pressure z=10: ', atm(1)%Pressure(10)
      WRITE(*,*)
    endif

    Error_Status = CRTM_Forward( Atm        , &
                                 Sfc        , &
                                 Geometry   , &
                                 ChannelInfo, &
                                 RTSolution , &
                                 Options = Options )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error in CRTM Forward Model'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF
    ! ============================================================================

    ! ============================================================================
    ! 6. **** Collecting output ****
    !
    ! User should read the user guide or the source code of the routine
    ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
    ! select the needed variables for outputs.  These variables are contained
    ! in the structure RTSolution.

    !---for file output, edited 2014.9.26
    do l = 1, n_channels_crtm_mySensor
        Tbsend(i,j,l) = real(RTSolution(l,1)%Brightness_Temperature)
        if ((i .EQ. 121 .AND. j .EQ. 1)) then
          WRITE(*,*) '  at x=121, y=1, Tbsend=',Tbsend(i,j,l)
        endif
        if (Tbsend(i,j,l) .NE. Tbsend(i,j,l)) then
          write(*,*) '  Tbsend is NaN at x=',x,' y=',y
          write(*,*) 'Pressures: ', atm(1)%Pressure
          write(*,*) 'Level_pressures: ', atm(1)%Level_Pressure
        endif
        if (Tbsend(i,j,l) .GT. 999) then
          WRITE(*,*) '  at x=',i,'y=',j,'Tbsend=',Tbsend(i,j,l)
          write(*,*) 'Pressures: ', atm(1)%Pressure
          write(*,*) 'Level_pressures: ', atm(1)%Level_Pressure
        endif
    enddo

    !--- end of x and y loops
    endif ! calcualte if this processor should
    enddo ! loop over x dimension
    enddo ! loop over y dimension

    CALL MPI_Allreduce(Tbsend,Tb_crtm,numx_crtm(isensor)*numy_crtm(isensor)*n_channels_crtm_mySensor,MPI_REAL,MPI_SUM,comm,ierr)

    ! ============================================================================

    ! ============================================================================
    !6.5  **** satellite beam convolution, writing the output ****

    if(my_proc_id==0) then

    ! write crtm BTs to file
    if (.false.) then
      do l = 1, n_channels_crtm_mySensor
        write(*,*) 'FILE_NAME should have sensor ', isensor, ' name ', trim(SENSOR_CLASS_NAMES(sensors_used(isensor)))
        write(FILE_NAME,'(a,a,a,a,i2.2,a,a)') trim(CRTM_OUT_DIR), '/Tb_', trim(SENSOR_CLASS_NAMES(sensors_used(isensor))), &
                                       '_ch', channel_numbers_crtm_mySensor(l), '_', trim(inputfile)
        write(*,*) 'saving TBs to ', trim(FILE_NAME)
        open(ens,file=trim(FILE_NAME),access='direct',recl=4)
        write(*,*) 'saving TBs file open'
        irec = 0
        do j = 1, numy_crtm(isensor)
        do i = 1, numx_crtm(isensor)
          irec= irec +1
          write( ens, rec=irec) Tb_crtm(i,j,l)
        enddo
        enddo
      enddo
       close (ens)
    endif  ! proc_0 writing simulated Tbs to file
    if(my_proc_id==0)  write(*,*) "done writing CRTM simulations to file"
    endif

    endif  ! whether the CRTM gets run for this sensor at all

    n_channels_all_mySensor = n_channels_all_eachSensor(isensor)
    ALLOCATE(channel_numbers_all_mySensor(n_channels_all_mySensor), STAT=Allocate_Status )
    channel_numbers_all_mySensor = channel_numbers_all_eachSensor(isensor,1:n_channels_all_mySensor)

    ich2 = 1
    mySensor_Tb_offset = 0
    do l = 1, isensor-1
      mySensor_Tb_offset = mySensor_Tb_offset + n_channels_all_eachSensor(l)
    enddo
    do l = 1, n_Channels_all_mySensor
      if (channel_numbers_all_mySensor(l) .EQ. channel_numbers_crtm_mySensor(ich2) ) then
        Tb(:,:,l+mySensor_Tb_offset) = Tb_crtm(:,:,ich2)
        ich2 = ich2 + 1
      endif
      if (ich2 .GT. n_channels_crtm_mySensor) exit
    enddo

    if(my_proc_id==0) then
      call cpu_time(cputime)
      WRITE(*,*) 'proc0 starting satellite beam convultion ', cputime
    end if
    if(my_proc_id==1) then
      call cpu_time(cputime)
      WRITE(*,*) 'proc1 starting satellite beam convultion ', cputime
    end if

    if(my_proc_id==0)  write(*,*) iob_radmin, nprocs, iob_radmax, isensor

    ALLOCATE(channel_search(n_Channels_all_mySensor), STAT = Allocate_Status)  ! array used to identify the channel number index from the channel number specified by the obs, for beam convolution
    do iob=iob_radmin+my_proc_id,iob_radmax,nprocs   ! parallelize across all procs by observation
      obstype = obs%type(iob)
      if ( obstype(1:9) .EQ. 'Microwave' .AND. sensor_indices(iob-iob_radmin+1) .EQ. isensor ) then
        obs_ii=obs%position(iob,1)-lon_radiance_min(isensor)+1
        obs_jj=obs%position(iob,2)-lat_radiance_min(isensor)+1
        channel_search = abs(obs%ch(iob) - channel_numbers_all_mySensor)
        ich = minloc(channel_search,1) + mySensor_Tb_offset
        xb_tb(iob) = beam_conv_gaussian_simple(numx_crtm(isensor), numy_crtm(isensor), Tb(:,:,ich), &
                         obs_ii, obs_jj, dx, dy, obs%efov_aScan(iob), obs%efov_cScan(iob) )
        if (xb_tb(iob) > 500 .OR. xb_tb(iob) .LE. 1 .OR. (xb_tb(iob) .NE. xb_tb(iob)) ) then
          WRITE(*,*) 'iob ', iob, ',  obs_ii ', obs_ii,      ',  obs_jj ', obs_jj
          WRITE(*,*) '        xb_tb(iob): ', xb_tb(iob) 
        endif
      endif
    enddo
    deallocate(channel_search)

    ! write summary of results
    if(my_proc_id==0) then
      WRITE(*,*) lon_radiance_min(isensor), lon_radiance_max(isensor), lat_radiance_min(isensor), lat_radiance_max(isensor)
      WRITE(*,'(a10,"   Tb=",f6.2,"~",f6.2)')inputfile,minval(Tb),maxval(Tb)
      WRITE(*,'(" Tb_conv=",f6.2,"~",f6.2)')minval(xb_tb(iob_radmin:iob_radmax),DIM=1,MASK=xb_tb(iob_radmin:iob_radmax) .GT. 0),maxval(xb_tb(iob_radmin:iob_radmax))
    endif
    do i = 10, nprocs, 10
      if(my_proc_id==i) then
        WRITE(*,'("proc_",i3.3," Tb_conv=",f6.2,"~",f6.2)')my_proc_id,minval(xb_tb(iob_radmin:iob_radmax),DIM=1,MASK=xb_tb(iob_radmin:iob_radmax) .GT. 0),maxval(xb_tb(iob_radmin:iob_radmax))
      endif
    enddo

    ! ============================================================================
    !  **** initializing all Tb and Tbsend fields ****
    !
    deallocate(channel_numbers_all_mySensor)

    if (n_Channels_crtm_mySensor .gt. 0) then
      Tbsend = 0.0
      Tb_crtm = 0.0
      CALL MPI_BCAST(Tbsend,numx_crtm(isensor)*numy_crtm(isensor)*max(n_channels_crtm_mySensor,1),MPI_REAL,0,comm,ierr)

      ! ============================================================================
      ! 7. **** DESTROY THE CRTM ****
      !
      deallocate(Tb_crtm)
      deallocate(Tbsend)
      deallocate(channel_numbers_crtm_mySensor)
      deallocate(RTSolution)

      Error_Status = CRTM_Destroy( ChannelInfo )
      IF ( Error_Status /= SUCCESS ) THEN
        Message = 'Error destroying CRTM'
        CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
        STOP
      END IF
    end if
    ! ============================================================================

  end do ! cycle to the next sensor

  if(my_proc_id==0) then
    call cpu_time(cputime)
    WRITE(*,*) 'finished with xb.f at ', cputime
  endif

  deallocate(Tb)
  deallocate(sensors_used)
  deallocate(lon_radiance_min)
  deallocate(lon_radiance_max)
  deallocate(lat_radiance_min)
  deallocate(lat_radiance_max)
  deallocate(scan_angle)
  deallocate(zenith_angle)
  deallocate(angle_sum_count)
  deallocate(sensor_indices)
  deallocate(Sensor_IDs_used)
  deallocate(numx_crtm)
  deallocate(numy_crtm)
  deallocate(n_channels_all_eachSensor)
  deallocate(channel_numbers_all_eachSensor)
  deallocate(channel_numbers_crtm)
  deallocate(sensor_index_crtm)
  deallocate(channel_numbers_crtm_eachSensor)
  deallocate(n_channels_crtm_eachSensor)

CONTAINS

  INCLUDE 'Load_CRTM_Structures_MW.inc'
  INCLUDE 'Load_CRTM_Structures_MW_slant.inc'

end subroutine xb_to_microwave
