!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Ensemble Kalman Filter
!
!   An EnKF data assimilation model using WRF
!   INITIATED BY FUQING ZHANG on 10/31/2001
!   Last modified  Fuqing Zhang 11/9/01
!   Last modified  Fuqing Zhang 12/31/02 for jan2000 storm enkf
!   Last modified by Altug Aksoy 04/02/2003 for sounding assimilation (output & sounding spacing parts)
!   Last modified by Altug Aksoy 04/09/2003 for parameterization of assimilated/updated variables
!   Last modified by Altug Aksoy 05/21/2003 for correlation coefficient changes in ENKF.f
!   Adapted to WRF by Zhiyong Meng  06/2005
!   Added Radar data assimilation by Yonghui Weng 03/2006
!   Added Hurricane Position and Intensity data assimilation by Yonghui Weng 01/2007
!   Added MPI by Yonghui Weng 02/2007
!   Enkf algorithm and MPI mechanism modified by Yue Ying 11/2012 to allow the use of more cpus
!   Further changes see CHANGE.log
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  program EnKF_WRF
  use constants
  use namelist_define
  use wrf_field_indices_define
  use mpi_module
  use mapinfo_define
  use netcdf
  use obs_define
  use map_utils
  use wrf_tools
  use radar
  use common_variables
  use boundary_variables
  integer            :: i_unit=80010, o_unit=90010
  integer            :: unit
  integer            :: ie, ic, it, len, i, j, k, n, ntot
  integer            :: ii, jj, kk, m, gid, sid, iid,jid, g_comm,s_comm
  real               :: member_per_cpu
  character (len=10) :: wrf_file
  double precision   :: time1, time2, time0, time_scat, time_ori
  real, allocatable, dimension(:) :: randnum
  real, allocatable, dimension(:,:,:,:,:) :: x2, x2_crnrs, x2_sn_bdy, x2_we_bdy
  type(proj_info)    :: proj
  character (len=24) :: finish_file_flg
  logical :: flag_exit
  real, allocatable, dimension(:,:,:,:,:) :: xf
  !----------------------------------------------------------------


  ! Initialize parallel stuff
  !---------------------------
  call parallel_start()



  ! Initialize file to hold time statistics
  ! ---------------------------------------
  if (my_proc_id == 0) then
    open(17, file='enkf_time.log')
    close(17)
  endif



  ! Initialize timer
  ! ----------------
  time0 = MPI_Wtime()




  ! Get WRF model information from the 1th member
  ! ----------------------------------------------
  ! ( reference to MM5_to_GrADS: module_wrf_to_grads_util.f)
  write(wrf_file,'(a5,i5.5)')'fort.',i_unit+1
  call get_wrf_info(wrf_file, ix, jx, kx, kxs, times, proj) ! saved wrf info to proj



  ! Read in namelist information
  ! ----------------------------
  call read_namelist(ix, jx, kx)





  ! Reading in observation files
  ! ----------------------------
  ! Note that the wrf file used here is NOT the mean field yet.
  write(wrf_file,'(a5,i5.5)')'fort.',i_unit+numbers_en+1
  if ( use_radar_rv == .false. .and. use_radar_rf == .false. ) then
    ! RAM efficient obs reading when no station radar is involved
    ! Read in obs data using only root
    if (my_proc_id == 0) &
      call get_all_obs(wrf_file, ix, jx, kx, times, proj)
    ! Broadcast obs data
    call MPI_BARRIER(comm, ierr)
    call bcast_obs(0)
  ! If radar is involved, all procs read in obs data
  else
    call get_all_obs(wrf_file, ix, jx, kx, times, proj)
  endif

  if ( print_detail > 1 .and. my_proc_id == 0 ) write(*,*)'All observations are loaded'
  if (my_proc_id==0) then
    write(*,'(a,x,f5.1,x,a)') 'Loaded namelist and obs in', MPI_Wtime()- time0, 'seconds.'
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30, x,f8.1,x,a,/)')'Loaded namelist and obs in', MPI_Wtime()- time0, 'seconds.'
    close(17)
  endif

  time1 = MPI_Wtime()


  !obs%num = 100


  ! Setting up parallelization and slab decomposition
  ! --------------------------------------------------
  ! Decompose domain based on n_cpus (nprocs) and n_members (numbers_en)
  ! figure out parallel strategy:
  if(.not.manual_parallel) then
    nicpu=1
    njcpu=1
    nmcpu=1
    do
      if ( nprocs/(nmcpu*nicpu*njcpu*2) < 1 ) exit
      if ( nicpu <= njcpu ) then
        nicpu=nicpu*2
      else
        njcpu=njcpu*2
      endif
    end do
  endif

  ! Enforcing situation where nmcpu is set to one.
  if ( manual_parallel .and. nmcpu .ne. 1 ) then
    if (my_proc_id == 0) write(*,*) "ERROR: Please set NMCPU=1. Exiting."
    call MPI_BARRIER(comm, ierr)
    call exit(STATUS)
  endif

  ! nicpu njcpu and nmcpu are from namelist if manual_parallel=true
  member_per_cpu=real(numbers_en+1)/nmcpu
  nm=ceiling(member_per_cpu)
  if ( my_proc_id == 0 ) then
    write(*,*) '---------------------------------------------------'
    write(*,'(a,i4,a,i5,a)') ' PARALLEL SCHEME: ',numbers_en,' ensemble members + 1 mean, ',nprocs,' cpus'
    write(*,'(a,i3,a,i3,a,i2,a)') ' domain decomposed to (ni*nj)=(',nicpu,'*',njcpu,') slabs'
    write(*,'(a,i3,a,f5.2,a,i3)') ' members distributed to ',nmcpu, ' cpu groups, member_per_cpu =', member_per_cpu, ', nm=',nm
    if ((nmcpu*nicpu*njcpu)/=nprocs) stop '*****PARALLEL ERROR: nprocs!=nmcpu*nicpu*njcpu *****'
    if (member_per_cpu<1) stop '*****PARALLEL ERROR: member_per_cpu cannot be less than 1! *****'
  endif

  ! determine which portion of grid my_proc_id is in charge of
  gid=int(my_proc_id/(nicpu*njcpu))
  sid=mod(my_proc_id,nicpu*njcpu)
  iid=mod(sid,nicpu)
  jid=int(sid/nicpu)
  call MPI_Comm_split(comm, gid, sid, s_comm, ierr)
  call MPI_Comm_split(comm, sid, gid, g_comm, ierr)

  ! Determine dimensions of decomposed slabs
  nv = 0
  do m=1,30
    if(len_trim(enkfvar(m))>=1) nv=nv+1
  enddo

  ! Two different ways to choose slab x-direction size
  if ( nicpu < int( sqrt( ix+1.) ) ) then
    ni=int((ix+1)/nicpu)+1
  else
    if ( mod( ix+1, nicpu ) == 0 ) then
      ni = int( (ix+1)/nicpu )
    else
      if (my_proc_id == 0) then
        write(*,*) 'When NICPU >= sqrt(ix+1), NICPU must be a factor of ix+1.'
        write(*,*) 'Current NICPU is not a factor of ix+1.'
        write(*,*) 'Please change NICPU to a factor of ix+1.'
        write(*,*) 'Exiting enkf.mpi.'
      endif
      call MPI_BARRIER(COMM, IERR)
      call exit(ierr)
    endif
  endif

  ! Two different ways to choose slab y-direction size
  if ( njcpu < int( sqrt( jx+1.) ) ) then
    nj=int((jx+1)/njcpu)+1
  else
    if ( mod( jx+1, njcpu ) == 0 ) then
      nj = int( (jx+1)/njcpu )
    else
      if (my_proc_id == 0) then
        write(*,*) 'When NJCPU >= sqrt(jx+1), NJCPU must be a factor of jx+1.'
        write(*,*) 'Current NJCPU is not a factor of jx+1.'
        write(*,*) 'Please change NJCPU to a factor of jx+1.'
        write(*,*) 'Exiting enkf.mpi.'
      endif
      call MPI_BARRIER(COMM, IERR)
      call exit(ierr)
    endif
  endif


  nk=kx+1
  ntot=ni*nj*nk*nv
  if(my_proc_id==0) then
     write(*,'(a,i4,a,i4,a,i4,a,i3,a,i3,a,i10)') ' On each cpu: ni*nj*nk*nv*nm=',ni,'*',nj,'*',nk,'*',nv,'*',nm,'=',ntot*nm
     if(ntot*nm>50000000) write(*,*) 'MEMORY WARNING: too much memory to be used, suggest use less cpu per node! Discard this if you already have.'
    write(*,*) '---------------------------------------------------'
    write(*,*) ' '
  endif

  ! Sanity check ni and nj
  ! Program will break if ni*nicpu < ix+1 or ni*nicpu >= ix+1+ni
  ! Likewise for the j-direction.
  flag_exit = .false.
  if ( ni*nicpu < ix+1  .or. ni*nicpu-(ix+1) >= ni ) then
    flag_exit=.true.
    if (my_proc_id == 0) then
      write(*,*)'nicpu = ',nicpu,' ni = ',ni,', ni*nicpu = ', ni*nicpu, ', ix + 1 =', ix+1
      write(*,*)'Choose a different NICPU.'
    endif
  else if (nj*njcpu < jx+1 .or. nj*njcpu-(jx+1) >= nj ) then
    flag_exit=.true.
    if (my_proc_id == 0) then
      write(*,*)'njcpu = ',njcpu,' nj = ',nj,', nj*njcpu = ', nj*njcpu, ', jx + 1 =', jx+1
      write(*,*)'Choose a different NJCPU.'
    endif
  endif

  if (flag_exit) call exit(ierr)


  ! Ensure variables needed for obs operators are present
  ! If now, extend enkfvar
  call supplement_enkfvar(nv)
  ntot=ni*nj*nk*nv

  ! Determine indices of enkfvar variables
  call find_enkfvar_indices(nv)





  ! Allocate common variables and wrf variables
  ! -------------------------------------------
  ! Constant variables
  allocate(xlat(ix,jx)); allocate(xlong(ix,jx)); allocate(lu_index(ix,jx))
  allocate(xland(ix,jx)); allocate(ind(obs%num))
  allocate(znu(kx)); allocate(znw(kx+1))
  ! Allocating slab arrays
  allocate(x(ni,nj,nk,nv,nm));  allocate(xm(ni,nj,nk,nv));
  allocate(xf(ni,nj,nk,nv,nm))
  x=0.;   xm=0.
  ! Allocating slab boundary arrays (used in compute_yf)
  allocate( x_we_bdy (2, nj, nk, nv, nm) )
  allocate( x_sn_bdy (2, ni, nk, nv, nm) )
  allocate( x_crnrs  (2,  2, nk, nv, nm) )
  allocate( xm_we_bdy(2, nj, nk, nv    ) )
  allocate( xm_sn_bdy(2, ni, nk, nv    ) )
  allocate( xm_crnrs (2,  2, nk, nv    ) )
  x_we_bdy = 0.; x_sn_bdy = 0.; x_crnrs  = 0.
  xm_we_bdy= 0.; xm_sn_bdy= 0.; xm_crnrs = 0.
  ! Allocating obs-related arrays
  allocate(ya(obs%num,numbers_en+1)); allocate(yasend(obs%num,numbers_en+1))
  allocate(yf(obs%num,numbers_en+1));
  allocate(yfm_radiance(obs%num)); allocate(yam_radiance(obs%num))
  ya = 0.; yasend = 0.; yf = 0.; yfm_radiance = 0.; yam_radiance  =0.



  ! Load constant arrays from domain
  ! --------------------------------
  call get_variable2d( wrf_file, 'XLONG     ', ix, jx, 1, xlong )
  call get_variable2d( wrf_file, 'XLAT      ', ix, jx, 1, xlat )
  call get_variable2d( wrf_file, 'LANDMASK  ', ix, jx, 1, xland )
  call get_variable2d( wrf_file, 'LU_INDEX  ', ix, jx, 1, lu_index )
  call get_variable1d( wrf_file, 'ZNW       ', kx+1, 1, znw )
  call get_variable1d( wrf_file, 'ZNU       ', kx, 1, znu )
  call get_variable0d( wrf_file, 'P_TOP     ', 1, p_top )

  if (my_proc_id==0) then
    write(*,'(a,x,f5.1,x,a)') 'Initialized slabs in', MPI_Wtime()- time1, 'seconds.'
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30, x,f8.1,x,a,/)')'Initialized slabs in', MPI_Wtime()- time1, 'seconds.'
    close(17)
  endif


  time1=MPI_Wtime()
  ! Read in initial Ensemble in x
  ! -----------------------------
  if ( my_proc_id == 0 ) write(*,*)'Read in initial Ensemble in x'
  call read_ensemble_bcast(i_unit)

  !call MPI_Barrier(comm, ierr)
  !if ( nmcpu == 1) then
  !  call read_ensemble_scat(i_unit)
  !else
  !  call read_ensemble_bcast(i_unit)
  !endif

  call MPI_Barrier( comm, ierr )
  if (my_proc_id==0) then
    write(*,'(a,x,f5.1,x,a)') 'Imported ensemble in', MPI_Wtime()- time1, 'seconds.'
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30, x,f8.1,x,a,/)')'Imported ensemble in', MPI_Wtime()- time1, 'seconds.'
    close(17)
  endif





  ! Calculate the ensemble mean
  ! --------------------------------------
  ! Slab mean
  xm = sum(x(:,:,:,:,1:numbers_en),5)/real(numbers_en)

  ! Boundary means
  call MPI_Allreduce( sum(x_we_bdy,5), xm_we_bdy, 2*nj*nk*nv, &
                      MPI_REAL, MPI_SUM, g_comm, ierr)
  call MPI_Allreduce( sum(x_sn_bdy,5), xm_sn_bdy, 2*ni*nk*nv, &
                      MPI_REAL, MPI_SUM, g_comm, ierr)
  call MPI_Allreduce( sum(x_crnrs ,5), xm_crnrs , 2*2*nk*nv, &
                      MPI_REAL, MPI_SUM, g_comm, ierr)
  xm_we_bdy = xm_we_bdy / real(numbers_en)
  xm_sn_bdy = xm_sn_bdy / real(numbers_en)
  xm_crnrs  = xm_crnrs  / real(numbers_en)

  ! Pass mean info into relevant parts of slab
  do n=1,nm
    ie=(n-1)*nmcpu+gid+1
    if(ie==numbers_en+1) then
      x(:,:,:,:,n)=xm
      x_we_bdy(:,:,:,:,n) = xm_we_bdy
      x_sn_bdy(:,:,:,:,n) = xm_sn_bdy
      x_crnrs(:,:,:,:,n)  = xm_crnrs
    endif
  enddo




  ! Outputting prior mean field
  ! ------------------------------
  call MPI_BARRIER(comm, ierr)
  if ( my_proc_id == 0 ) write(*,*)'Outputting the prior mean field'
  if(gid==0) &
    call output_RAM_efficient(i_unit+numbers_en+1,ix,jx,kx,kxs,ni,nj,nk,nv,nm,s_comm,sid,iid,jid,xm,times)

  call MPI_BARRIER( comm, ierr)





  ! Preparing observation assimilation order
  ! ---------------------------------------

  if(my_proc_id==0) then
    do it=1,obs%num
      ind(it)=it
    enddo

    ! If we want to assimilate obs with random order, random_order=.true.
    if(random_order) then
      allocate(randnum(obs%num))
      call random_seed()
      call random_number(randnum)
      call quicksort(obs%num,randnum,ind)
      deallocate(randnum)
    endif

  endif
  call MPI_Bcast(ind,obs%num,MPI_INTEGER,0,comm,ierr)


  ! If we want to assimilate nonlinear obs in batches s.t. obs in each batch
  ! are more than 1 HROI apart, set batchwise_order=.true.
  ! Written by Man-Yau Chan.
  if ( batchwise_order .or. use_nonlinear_enkf .or. use_gmm_enkf ) then

    ! Assign batch ids
    !if(my_proc_id==0) write(*,*) 'meow'
    call assign_obs_batch_id( obs )

    ! Rearrange observations by batch id
    call sort_obs_ascending_batch_id( obs )

  endif

  ! If we are using joint-space gmm_enkf, only have 1 batch
  if ( use_jointspace_gmm_enkf ) then
    obs%batch_id(:) = 1
    obs%num_batches = 1
  endif





  if (my_proc_id==0) then
    write(*,'(a,x,f5.1,x,a)') '<x> and obs order made in', MPI_Wtime()- time1, 'seconds.'
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30, x,f8.1,x,a,/)')'<x> and obs order made in', MPI_Wtime()- time1, 'seconds.'
    close(17)
  endif
  time1=MPI_Wtime()






  ! BT-based inflation if use_abei by Minamide
  ! ------------------------------------------
  perform_abei: if (use_abei .and. (raw%radiance%num.ne.0)) then
    time1 = MPI_Wtime()
    if ( my_proc_id == 0 ) write(*,*)'Running ABEI ...'
    write(wrf_file,'(a5,i5.5)')'fort.',i_unit+1
    call cal_abei(wrf_file,ix,jx,kx,kxs,ni,nj,nk,nv,nm,g_comm,s_comm,gid,sid,iid,jid,xm,x,ind,proj,times)

    ! Preparing to recalculate ensemble mean
    if ( my_proc_id == 0 ) write(*,*) 'Recalculating ensemble mean'
    ! Initialize arrays
    xm = 0.
    do n = 1, nm
      ie=(n-1)*nmcpu+gid+1
      if(ie==numbers_en+1) then
        x(:,:,:,:,n) = 0.
      endif
    enddo

    ! Calculate mean
    xm = sum( x(:,:,:,:,1:numbers_en), 5 ) / real(numbers_en)
    x(:,:,:,:,numbers_en+1) = xm

    call MPI_BARRIER( COMM, IERR)
    if ( my_proc_id == 0 ) write(*,*)' ABEI took',MPI_Wtime()-time1, 'seconds'
    if (my_proc_id==0) then
      open( 17, file='enkf_time.log', status='old',position='append', action='write')
      write(17,'(a30, x,f8.1,x,a,/)')'ABEI took', MPI_Wtime()- time1, 'seconds.'
      close(17)
    endif
    time1=MPI_Wtime()

    ! Output mean field (ABEI sometimes shift the mean field by a few bits)
    call MPI_BARRIER(comm, ierr)
    if ( my_proc_id == 0 ) write(*,*)'Outputting the prior mean field'
    if(gid==0) &
      call output_RAM_efficient(i_unit+numbers_en+1,ix,jx,kx,kxs,ni,nj,nk,nv,nm,s_comm,sid,iid,jid,xm,times)

    call MPI_BARRIER( comm, ierr)


  endif perform_ABEI




  ! Construct simulated observation ensemble
  time1=MPI_Wtime()
  flag_level_assignment=.true.
  call refresh_slab_boundaries()
  call compute_yf( wrf_file, g_comm, s_comm, gid, sid, iid, jid, proj )
  call MPI_BARRIER( COMM, IERR)
  if ( my_proc_id == 0 ) write(*,*)' Sim obs generation took ', MPI_Wtime() - time1, ' seconds.'
  if ( my_proc_id == 0) then
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30, x,f8.1,x,a,/)')'yf computation took', MPI_Wtime()- time1, 'seconds.'
    close(17)
  endif



  ! Time taken for pre-EnKF stuff
  ! -----------------------------
  time1 = MPI_Wtime()
  if ( my_proc_id == 0 ) write(*,*)' Pre-EnKF processing took ', time1-time0, ' seconds.'




  ! Apply namelist-specified inflation factor (for no inflation, inflate=1)
  ! ----------------------------------------------------------------------
  do n = 1, nm
    x(:,:,:,:,n)=inflate*(x(:,:,:,:,n)-xm)
  enddo




  ! Making copy of inflated prior and simulated obs
  ! -----------------------------------------------
  xf = x; yf = ya




  ! Run Ensemble Kalman Filter (various versions being added)
  ! ---------------------------------------------------------
  write(wrf_file,'(a5,i5.5)')'fort.',i_unit+1
  time1 = MPI_Wtime()

  ! Nonlinear batchwise EnKF (aka sequential filtering)
  if ( use_nonlinear_enkf .or. use_gmm_enkf ) then
    call batchwise_da( wrf_file, g_comm, s_comm, gid, sid, &
                         iid, jid, proj )

  ! Vanilla extended state EnKF (aka parallel EnKF filtering)
  else
    call enkf( wrf_file, g_comm, s_comm, gid, sid, iid, jid, proj)
  endif

  

  ! Report time to run the EnKF
  time2 = MPI_Wtime()
  if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' EnKF took ', time2-time1, ' seconds.'
  if ( my_proc_id == 0) then
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30, x,f8.1,x,a,/)')'EnKF took', time2- time1, 'seconds.'
    close(17)
  endif



  ! Write filter diagnostics to fort.10000 and fort.10001 files
  ! -----------------------------------------------------------
  call write_filter_diagnostics( proj )



  ! Perform covariance relaxation (e.g., RTPP)
  ! ------------------------------------------
  call covariance_relaxation( xf, wrf_file, g_comm, s_comm, gid, sid, iid, jid )


  ! Remove negative mixing ratios
  ! -----------------------------
  call minamide_neg_q_removal( g_comm, s_comm, gid, sid, iid, jid )



  ! Save the Analyses
  !------------------
  if ( my_proc_id == 0 ) write(*,*)'Output members and mean now...'
  time1 = MPI_Wtime()
  ! Using RAM efficient method when NMCPU == numbers_en+1
  if (nmcpu == numbers_en+1 ) then
    do n =1, nm
      ie = (n-1)*nmcpu+gid+1
      if ( ie <= numbers_en+1 ) then
        call output_RAM_efficient(o_unit+ie,ix,jx,kx,kxs,ni,nj,nk,nv,nm,s_comm,sid,iid,jid,x(:,:,:,:,n),times)
      endif
    enddo

  ! If nmcpu == 1, use MPI_GATHER to output WRF files
  else if ( nmcpu == 1 ) then
    call output_full_ensemble_gather( o_unit )

  ! Otherwise, if NMCPU .ne. numbers_en+1, use mpi_reduce
  else if (nmcpu .ne. numbers_en+ 1) then
    call output_full_ensemble_mpireduce( o_unit, sid, gid, iid, jid )
  endif

  time2 = MPI_Wtime()
  if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' Output took', time2-time1, ' seconds.'
  if ( my_proc_id == 0) then
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30, x,f8.1,x,a,/)')'Ens output took', time2- time1, 'seconds.'
    close(17)
  endif


  ! Safety barrier
  call MPI_BARRIER(comm, ierr)

  ! Give a finish flag
  finish_file_flg = times(1:4)//times(6:7)//times(9:10)//times(12:13)//times(15:16)//'.finish_flag'
  do k = 1, 12
    if ( finish_file_flg(k:k) .eq. " " ) finish_file_flg(k:k) ='0'
  enddo
  if ( my_proc_id == 0 ) then
    open(10,file=finish_file_flg)
    write(10,*)times
    close(10)
  endif

  ! Clean up
  deallocate(ind)
  deallocate(x)
  deallocate(xm)
  call MPI_Comm_free(g_comm,ierr)
  call MPI_Comm_free(s_comm,ierr)
  call parallel_finish()
  if ( my_proc_id == 0 ) write(*,*)' Output took ', time_end-time2, ' seconds.'
  if ( my_proc_id == 0 ) write(*,*)' All took ', time_end-time_start, ' seconds.'
  if ( my_proc_id == 0 ) write(*,'(a)')' Successful completion of EnKF.'
  if ( my_proc_id == 0) then
    open( 17, file='enkf_time.log', status='old',position='append', action='write')
    write(17,'(a30, x,f8.1,x,a,/)')'Entire enkf.mpi took', MPI_Wtime() - time_start, 'seconds.'
    close(17)
  endif

end program EnKF_WRF
