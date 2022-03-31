! ===================================================================
! Subroutines to assimilate obs in batches of non-overlapping obs.
! I.e., obs within the same batch are separated by at least 1 HROI
! in distance.
!
! Useful for obs with nonlinear obs operators.
!
! ===================================================================
! Written by Man-Yau Chan
! ===================================================================

! ===================================================================
! SUPPORTED DATA ASSIMILATION METHODS:
! ------------------------------------
! 1) Basic PSU-EnKF enkf subroutine.
!    Note: While subroutine enkf is written in joint-space form, the 
!          obs space update is effectively ignored in batchwise mode.
!
! 2) Subroutine gmm_enkf
!    Algorithm is quite similar to Chan et al 2020: An efficient
!    bi-Gaussian ensemble Kalman filter for satellite data 
!    assimilation



! ===================================================================
! IMPORTANT NOTES:
! ----------------
! 1) This system has only been designed with geostationary IR obs for
!    now. Can be modified to handle other nonlinear obs later.
! 2) Subroutine batchwise_enkf is the subroutine that evokes the full
!    batch-wise EnKF. This batchwise_enkf subroutine is called in
!    main.f.
! ===================================================================



! ==================================================================
! PROCEDURE:
! ----------
! 1) Assign batch ids to the obs (subroutine assign_obs_batch)
!    For now, all non-IR obs are assigned batch ID of -888888
!    IR obs are assigned non-zero positive batch IDs.
!
! 2) Make a master copy of variable obs
!    (subroutine copy_obs_to_master_obs)
!
! For each batch id (i.e., do b_id = 0, max_batch_id)
! 3) Reconstruct variable obs with only observations with batch id
!    b_id (subroutine select_obs_subset)
!
! 4) Reallocate yf and ya
!
! 5) Run compute_yf subroutine for subsetted obs.
!
! 6) Run enkf subroutine for subsetted obs.
!
! ==================================================================



! ==================================================================
! Subroutine to perform batchwise data assimilation
subroutine batchwise_da(wrf_file0, g_comm0,s_comm0,gid0,sid0,iid0,jid0, proj0)

  ! Load modules
  use constants
  use namelist_define
  use mpi_module
  use obs_define
  use netcdf
  use map_utils
  use common_variables
  use wrf_field_indices_define
  use gauss_mix_model_enkf

  ! Variable definitions
  implicit none
  integer, intent(in) :: sid0,gid0,iid0,jid0,g_comm0,s_comm0
  character (len=10), intent(in) :: wrf_file0
  type (proj_info)    :: proj0

  integer :: b_id, b_id0


  ! 2) Make master copy of obs
  ! --------------------------
  call copy_obs_type_variable( obs, master_obs )
  batch_start_obs_num = 0

  ! Stop-gap measures to side-step memory corruption
  !allocate( obs_mask( ix, jx) )    

  
  ! Revert x from perturbations to full state
  do n = 1, nm
    ie = (n-1)*nmcpu+gid0+1
    if ( ie <= numbers_en+1 ) then
      x(:,:,:,:,ie) = x(:,:,:,:,ie) + xm(:,:,:,:)
    endif
  enddo 



  ! Iterating over all batch ids
  ! -----------------------------
  if( my_proc_id == 0) write(*,*) 'master_obs%num_batches', master_obs%num_batches
  batch_da_loop: do b_id0 = 1, master_obs%num_batches

    !! Special handling for obs not handled by the batch assignment
    !if (b_id0 == 0) then
    !  b_id = -888888
    !else
    !  b_id = b_id0
    !endif
    b_id = b_id0



    ! 3) Subsetting observation
    ! -------------------------
    call select_obs_subset( b_id )
    call MPI_BARRIER( comm, ierr )




    ! 4) Reallocate ya and yasend
    ! ---------------------------
    deallocate( ya     ); allocate( ya    (obs%num, numbers_en+1) )
    deallocate( yasend ); allocate( yasend(obs%num, numbers_en+1) )


    ! 5) Running observation operator
    ! -------------------------------
    flag_level_assignment = .false.
    call refresh_slab_boundaries()
    call compute_yf( wrf_file0, g_comm0, s_comm0, gid0, sid0, iid0, jid0, proj0 )




    ! 6) Assigning clusters to members for GMM-EnKF
    ! ---------------------------------------------
    ! Heuristic 3-kernel clustering procedure
    call tag_clear_congestus_deep_clusters( iid0, jid0 )

    ! Remove clustering for non-IR obs
    do iob = 1, obs%num 
      if ( trim(obs%type(iob)) .ne. 'Radiance' ) kernel_id(iob,:) = 2
    enddo


    ! Special handling for two-kernel situation
    if (use_gmm_enkf .and. max_kernel_num == 2) &
      where ( kernel_id > 2 ) kernel_id = 2

    ! Special handling for one-kernel situations
    if ( use_gmm_enkf .and. max_kernel_num == 1 ) kernel_id = 1

    ! Special handling if use_nonlinear_enkf=.true. and use_gmm_enkf=.false.
    if ( use_nonlinear_enkf ) then
      kernel_id = 1
      max_kernel_num=1
    endif
    call MPI_BARRIER(comm, ierr)


    ! 7) Running DA subroutine
    ! ------------------------
    ! Gaussian mixture model EnKF subroutine
    ! (can be used for batchwise enkf as well!)
    call gmm_enkf( wrf_file0, g_comm0, s_comm0, gid0, sid0, iid0, jid0 )

    ! Recalculate ensemble mean
    xm = sum( x(:,:,:,:,1:numbers_en), 5) / real(numbers_en)
    !do n = 1, nm
    !  ie = (n-1)*nmcpu+gid0+1
    !  if ( ie == numbers_en+1 )  x(:,:,:,:,ie) = xm
    !enddo
    x(:,:,:,:,numbers_en+1) = xm

    call MPI_BARRIER( comm, ierr )

    call copy_obs_subset_to_master( b_id )


    call MPI_BARRIER( comm, ierr )
    batch_start_obs_num = batch_start_obs_num + obs%num

  enddo batch_da_loop


  ! Not doing any of the allocation and deallocation stuff for jointspace
  ! gmm-enkf
  if (.not. use_jointspace_gmm_enkf ) then
    ! Copy master_obs to variable obs
    call copy_obs_type_variable( master_obs, obs )

    ! Reallocate ya to prepare to rerun obs operator
    deallocate( ya     ); allocate( ya    (obs%num, numbers_en+1) )
    deallocate( yasend ); allocate( yasend(obs%num, numbers_en+1) )
  endif

  !if ( my_proc_id == 0) then
  !  write(*,*) x(50,50,1:3,ind_pb,1)
  !endif
  ! Run obs operator
  flag_level_assignment=.false.
  call refresh_slab_boundaries()
  if ( .not. use_jointspace_gmm_enkf ) &
    call compute_yf( wrf_file0, g_comm0, s_comm0, gid0, sid0, iid0, jid0, proj0 )


  ! Converting x from array of full states to array of perturbations
  do n = 1, nm
    ie = (n-1)*nmcpu+gid0+1
    if ( ie <= numbers_en+1 ) then
      x(:,:,:,:,ie)= x(:,:,:,:,ie) - xm(:,:,:,:)
    endif
  enddo 



end subroutine batchwise_da
! =================================================================








! ==================================================================
! Subroutine to select observations with a specific batch id (b_id)
subroutine select_obs_subset( b_id )

  ! Load modules
  use constants
  use namelist_define
  use mpi_module
  use obs_define
  use netcdf
  use map_utils


  ! Variable definitions
  implicit none
  integer, intent(in) :: b_id
  integer :: iob, iob2


  ! 1) Deallocate obs variable
  ! --------------------------
  call deallocate_obs_type( obs )
  obs%num = 0


  ! 2) Count number of obs with desired batch id
  ! --------------------------------------------
  count_batch_obs: do iob = 1, master_obs%num
    if ( master_obs%batch_id(iob) == b_id ) obs%num = obs%num + 1
  enddo count_batch_obs


  ! 3) Reallocate obs variable
  ! --------------------------
  call allocate_obs_type( obs, obs%num )


  ! 4) Copy obs with desired batch id into variable obs
  iob2 = 1
  copy_batch_obs: do iob = 1, master_obs%num

    if ( master_obs%batch_id(iob) == b_id ) then

      ! Make copy
      call copy_obs_type_entry( master_obs, iob, obs, iob2 )

      ! increment obs relative counter
      iob2 = iob2 + 1

    endif

  enddo copy_batch_obs


  ! 5) Counting the number of radiance observations
  ! -----------------------------------------------
  raw%radiance%num = 0
  do iob = 1, obs%num
    if ( trim( obs%type(iob) ) == 'Radiance' ) &
      raw%radiance%num = raw%radiance%num + 1
  enddo
  

end subroutine select_obs_subset
! ==================================================================





! ==================================================================
! Subroutine to copy from assimilated obs to master_obs
subroutine copy_obs_subset_to_master( b_id )

  ! Load modules
  use constants
  use namelist_define
  use mpi_module
  use obs_define
  use netcdf
  use map_utils


  ! Variable definitions
  implicit none
  integer, intent(in) :: b_id
  integer :: iob, iob2, start_iob


  ! 1) Find index of batch b_id's first obs in master_obs
  ! -----------------------------------------------------
  seek_batch_first_obs: do iob = 1, master_obs%num
    if ( master_obs%batch_id(iob) == b_id ) exit seek_batch_first_obs
  enddo seek_batch_first_obs
  start_iob = iob

  ! 2) Copy over subsets
  ! ----------------------
  copy_subset_to_master: do iob2 = 1, obs%num

    call copy_obs_type_entry( obs, iob2, master_obs, start_iob + iob2-1 )

  enddo copy_subset_to_master

 

end subroutine copy_obs_subset_to_master
! ==================================================================

! ==================================================================



! =====================================================================================
! Subroutine to refresh boundary variables (e.g., x_we_bdy) 
! Notes:
! 1) Required for batchwise DA and GMM-EnKF. B'cos boundary variables are needed for 
!    compute_yf, but batchwise DA and GMM-EnKF does not update boundary variables
! 2) THIS SUBROUTINE CALLS MPI-USING SUBROUTINES. HANDLE WITH CARE!
! =====================================================================================
subroutine refresh_slab_boundaries()


  ! Load derived types
  use namelist_define
  use wrf_tools
  use mpi_module
  use common_variables
  use boundary_variables

  implicit none
  integer :: sid, iid, jid ! ids to indicate slab position
  integer :: target_id    ! Process id of receiving process

  ! Pass information eastwards
  call communicate_boundary_eastward()

  ! Pass information westwards
  call communicate_boundary_westward()

  ! Pass information northwards
  call communicate_boundary_northward()

  ! Pass information southwards
  call communicate_boundary_southward()

  ! Pass information regarding corners
  call communicate_boundary_corners()

  ! Recalculate the boundary mean values
  xm_we_bdy = sum( x_we_bdy(:,:,:,:,1:numbers_en), 5 ) / numbers_en
  xm_sn_bdy = sum( x_sn_bdy(:,:,:,:,1:numbers_en), 5 ) / numbers_en
  xm_crnrs  = sum( x_crnrs (:,:,:,:,1:numbers_en), 5 ) / numbers_en
  x_we_bdy(:,:,:,:,numbers_en+1) = xm_we_bdy
  x_sn_bdy(:,:,:,:,numbers_en+1) = xm_sn_bdy
  x_crnrs (:,:,:,:,numbers_en+1) = xm_crnrs

end subroutine refresh_slab_boundaries
! =====================================================================================




! =====================================================================================
! Subroutine to pass boundary information eastwards
! THIS IS AN MPI SUBROUTINE. HANDLE WITH CARE!
! 1) Figure out sid, gid, iid and jid for each process
! 2) Pass information eastwards by:
!    a) Slabs with iid = 0, 2, 4, ... transmit x(ni,:,:,:,:) to x_we_bdy(1,:,:,:,:) of
!       slabs with iid = 1, 3, 5, ...
!    b) Slabs with iid = 1, 3, 5, ... transmit x(ni,:,:,:,:) to x_we_bdy(1,:,:,:,:) of
!       slabs with iid = 2, 4, 6, ...
subroutine communicate_boundary_eastward()

  ! Load derived types
  use namelist_define
  use wrf_tools
  use mpi_module
  use common_variables
  use boundary_variables

  implicit none
  integer :: sid, iid, jid ! ids to indicate slab position
  integer :: target_id    ! Process id of receiving process
  real, dimension(nj, nk, nv, nm) :: buffer_we ! Buffer used to pass edges eastward


  ! Step 1: Figure out process sid, gid, iid, jid 
  ! ---------------------------------------------
  sid=mod(my_proc_id,nicpu*njcpu)
  iid=mod(sid,nicpu)  ! Zonal index for the slab
  jid=int(sid/nicpu)  ! Meridional index for the slab

  ! Step 2a: Eastward information passing (phase a)
  ! -----------------------------------------------
  ! Processes with iid = 0, 2, 4, ... send out information
  eastward_a_send: if ( mod( iid, 2 ) == 0 ) then

    ! Figure out targetted process
    ! Note that because nmcpu = 1, sid = my_proc_id
    target_id = jid * nicpu + iid + 1

    ! Pack info into buffer
    buffer_we = x( ni,:,:,:,: )
    
    ! Send information to target process
    ! Note that iid = nicpu-1 is at the edge of domain. So no sending from iid=nicpu-1
    if ( iid < nicpu-1) &
      call MPI_SEND(buffer_we, nj*nk*nv*nm, MPI_REAL, target_id, my_proc_id, comm, ierr)

  endif eastward_a_send


  ! Processes with iid = 1, 3, 5, .... receive info
  eastward_a_recv: if ( mod( iid, 2 ) == 1 ) then

    ! Figure out the process that is sending info
    target_id = jid * nicpu + iid -1


    ! Receive information
    call MPI_RECV( buffer_we, nj*nk*nv*nm, MPI_REAL, target_id, target_id, comm, &
                   MPI_STATUS_IGNORE, ierr)
    
    ! Put info into boundary variable
    x_we_bdy(1,:,:,:,:) = buffer_we

  endif eastward_a_recv


  ! Wait for all processes to reach this point
  call MPI_Barrier( comm, ierr )


  ! Step 2b: Eastward information passing (phase b)
  ! -----------------------------------------------
  ! Processes with iid = 1, 3, 5 .... send out info to iid = 2, 4, 6
  eastward_b_send: if ( mod(iid, 2) == 1) then

    ! Figure out the receiving process id
    target_id = jid * nicpu + iid + 1

    ! Pack info into buffer
    buffer_we = x(ni,:,:,:,:)

    ! Send info to target process
    ! Note that iid = nicpu-1 is at the edge of domain. So no sending from iid=nicpu-1
    if ( iid < nicpu - 1 ) &
      call MPI_SEND(buffer_we, nj*nk*nv*nm, MPI_REAL, target_id, my_proc_id, comm, ierr)

  endif eastward_b_send 

  ! Processes with iid = 2, 4, 6, .... receive info
  eastward_b_recv: if ( mod( iid, 2 ) == 0 .and. iid >= 2 ) then

    ! Figure out the process that is sending info
    target_id = jid * nicpu + iid-1

    ! Receive information
    call MPI_RECV( buffer_we, nj*nk*nv*nm, MPI_REAL, target_id, target_id, comm, &
                   MPI_STATUS_IGNORE, ierr)

    ! Put info into boundary variable
    x_we_bdy(1,:,:,:,:) = buffer_we

  endif eastward_b_recv

  call MPI_BARRIER( comm, ierr )


end subroutine communicate_boundary_eastward
! =====================================================================================




! =====================================================================================
! Subroutine to pass boundary information westwards
! THIS IS AN MPI SUBROUTINE. HANDLE WITH CARE!
! 1) Figure out sid, gid, iid and jid for each process
! 2) Pass information westwards by:
!    a) Slabs with iid = 1, 3, 5, ... transmit x(1 ,:,:,:,:) to x_we_bdy(2,:,:,:,:) of
!       slabs with iid = 0, 2, 4, ...
!    b) Slabs with iid = 2, 4, 6, ... transmit x(1 ,:,:,:,:) to x_we_bdy(2,:,:,:,:) of
!       slabs with iid = 1, 3, 5, ...
subroutine communicate_boundary_westward()

  ! Load derived types
  use namelist_define
  use wrf_tools
  use mpi_module
  use common_variables
  use boundary_variables

  implicit none
  integer :: sid, iid, jid ! ids to indicate slab position
  integer :: target_id    ! Process id of receiving process
  real, dimension(nj, nk, nv, nm) :: buffer_we ! Buffer used to pass edges eastward


  ! Step 1: Figure out process sid, gid, iid, jid 
  ! ---------------------------------------------
  sid=mod(my_proc_id,nicpu*njcpu)
  iid=mod(sid,nicpu)  ! Zonal index for the slab
  jid=int(sid/nicpu)  ! Meridional index for the slab

  ! Step 2a: Westward information passing (phase a)
  ! -----------------------------------------------
  ! Processes with iid = 1, 3, 5, ... send out information
  westward_a_send: if ( mod( iid, 2 ) == 1 ) then

    ! Figure out targetted process
    ! Note that because nmcpu = 1, sid = my_proc_id
    target_id = jid * nicpu + iid - 1

    ! Pack info into buffer
    buffer_we = x( 1,:,:,:,: )
    
    ! Send information to target process
    call MPI_SEND(buffer_we, nj*nk*nv*nm, MPI_REAL, target_id, my_proc_id, comm, ierr)

  endif westward_a_send


  ! Processes with iid = 0, 2, 4, .... receive info
  westward_a_recv: if ( mod( iid, 2 ) == 0 ) then

    ! Figure out the process that is sending info
    target_id = jid * nicpu + iid +1

    ! Receive information
    if ( iid < nicpu-1 ) then
      call MPI_RECV( buffer_we, nj*nk*nv*nm, MPI_REAL, target_id, target_id, comm, &
                     MPI_STATUS_IGNORE, ierr)
    
      ! Put info into boundary variable
      x_we_bdy(2,:,:,:,:) = buffer_we
    endif

  endif westward_a_recv


  ! Wait for all processes to reach this point
  call MPI_Barrier( comm, ierr )


  ! Step 2b: Westward information passing (phase b)
  ! -----------------------------------------------
  ! Processes with iid = 2, 4, 6 .... send out info to iid = 1, 3, 5
  westward_b_send: if ( mod(iid, 2) == 0 .and. iid >= 2) then

    ! Figure out the receiving process id
    target_id = jid * nicpu + iid - 1

    ! Pack info into buffer
    buffer_we = x(1,:,:,:,:)

    ! Send info to target process
    call MPI_SEND(buffer_we, nj*nk*nv*nm, MPI_REAL, target_id, my_proc_id, comm, ierr)

  endif westward_b_send 

  ! Processes with iid = 1, 3, 5, .... receive info
  westward_b_recv: if ( mod( iid, 2 ) == 1 ) then

    ! Figure out the process that is sending info
    target_id = jid * nicpu + iid + 1

    ! Receive information
    if ( iid < nicpu-1 ) then
      call MPI_RECV( buffer_we, nj*nk*nv*nm, MPI_REAL, target_id, target_id, comm, &
                     MPI_STATUS_IGNORE, ierr)

      ! Put info into boundary variable
      x_we_bdy(2,:,:,:,:) = buffer_we
    endif

  endif westward_b_recv

  call MPI_BARRIER( comm, ierr )


end subroutine communicate_boundary_westward
! =====================================================================================




! =====================================================================================
! Subroutine to pass boundary information northwards
! THIS IS AN MPI SUBROUTINE. HANDLE WITH CARE!
! 1) Figure out sid, gid, iid and jid for each process
! 2) Pass information northwards by:
!    a) Slabs with jid = 0, 2, 4, ... transmit x(:,nj,:,:,:) to x_sn_bdy(1,:,:,:,:) of
!       slabs with jid = 1, 3, 5, ...
!    b) Slabs with jid = 1, 3, 5, ... transmit x(:,nj,:,:,:) to x_sn_bdy(1,:,:,:,:) of
!       slabs with jid = 2, 4, 6, ...
subroutine communicate_boundary_northward()

  ! Load derived types
  use namelist_define
  use wrf_tools
  use mpi_module
  use common_variables
  use boundary_variables

  implicit none
  integer :: sid, iid, jid ! ids to indicate slab position
  integer :: target_id    ! Process id of receiving process
  real, dimension(ni, nk, nv, nm) :: buffer_sn ! Buffer used to pass edges northward


  ! Step 1: Figure out process sid, gid, iid, jid 
  ! ---------------------------------------------
  sid=mod(my_proc_id,nicpu*njcpu)
  iid=mod(sid,nicpu)  ! Zonal index for the slab
  jid=int(sid/nicpu)  ! Meridional index for the slab

  ! Step 2a: Northward information passing (phase a)
  ! -----------------------------------------------
  ! Processes with jid = 0, 2, 4, ... send out information
  northward_a_send: if ( mod( jid, 2 ) == 0 ) then

    ! Figure out targetted process
    ! Note that because nmcpu = 1, sid = my_proc_id
    target_id = (jid+1) * nicpu + iid

    ! Pack info into buffer
    buffer_sn = x( :,nj,:,:,: )
    
    ! Send information to target process
    ! Note that jid = njcpu-1 is at the edge of domain. So no sending from jid=njcpu-1
    if ( jid < njcpu-1) &
      call MPI_SEND(buffer_sn, ni*nk*nv*nm, MPI_REAL, target_id, my_proc_id, comm, ierr)

  endif northward_a_send


  ! Processes with gid = 1, 3, 5, .... receive info
  northward_a_recv: if ( mod( jid, 2 ) == 1 ) then

    ! Figure out the process that is sending info
    target_id = (jid-1) * nicpu + iid 

    ! Receive information
    call MPI_RECV( buffer_sn, ni*nk*nv*nm, MPI_REAL, target_id, target_id, comm, &
                   MPI_STATUS_IGNORE, ierr)
    
    ! Put info into boundary variable
    x_sn_bdy(1,:,:,:,:) = buffer_sn

  endif northward_a_recv


  ! Wait for all processes to reach this point
  call MPI_Barrier( comm, ierr )


  ! Step 2b: Northward information passing (phase b)
  ! -----------------------------------------------
  ! Processes with jid = 1, 3, 5 .... send out info to jid = 2, 4, 6
  northward_b_send: if ( mod(jid, 2) == 1) then

    ! Figure out the receiving process id
    target_id = (jid+1) * nicpu + iid 

    ! Pack info into buffer
    buffer_sn = x(:,nj,:,:,:)

    ! Send info to target process
    ! Note that iid = nicpu-1 is at the edge of domain. So no sending from iid=nicpu-1
    if ( jid < njcpu - 1 ) &
      call MPI_SEND(buffer_sn, ni*nk*nv*nm, MPI_REAL, target_id, my_proc_id, comm, ierr)

  endif northward_b_send 

  ! Processes with jid = 2, 4, 6, .... receive info
  northward_b_recv: if ( mod( jid, 2 ) == 0 .and. jid >= 2 ) then

    ! Figure out the process that is sending info
    target_id = (jid-1) * nicpu + iid

    ! Receive information
    call MPI_RECV( buffer_sn, ni*nk*nv*nm, MPI_REAL, target_id, target_id, comm, &
                   MPI_STATUS_IGNORE, ierr)

    ! Put info into boundary variable
    x_sn_bdy(1,:,:,:,:) = buffer_sn

  endif northward_b_recv

  call MPI_BARRIER( comm, ierr )


end subroutine communicate_boundary_northward
! =====================================================================================




! =====================================================================================
! Subroutine to pass boundary information southwards
! THIS IS AN MPI SUBROUTINE. HANDLE WITH CARE!
! 1) Figure out sid, gid, iid and jid for each process
! 2) Pass information southwards by:
!    a) Slabs with jid = 1, 3, 5, ... transmit x(:,1,:,:,:) to x_sn_bdy(2,:,:,:,:) of
!       slabs with jid = 0, 2, 4, ...
!    b) Slabs with jid = 2, 4, 6, ... transmit x(:,1,:,:,:) to x_sn_bdy(2,:,:,:,:) of
!       slabs with jid = 1, 3, 5, ...
subroutine communicate_boundary_southward()

  ! Load derived types
  use namelist_define
  use wrf_tools
  use mpi_module
  use common_variables
  use boundary_variables

  implicit none
  integer :: sid, iid, jid ! ids to indicate slab position
  integer :: target_id    ! Process id of receiving process
  real, dimension(ni, nk, nv, nm) :: buffer_sn ! Buffer used to pass edges eastward


  ! Step 1: Figure out process sid, gid, iid, jid 
  ! ---------------------------------------------
  sid=mod(my_proc_id,nicpu*njcpu)
  iid=mod(sid,nicpu)  ! Zonal index for the slab
  jid=int(sid/nicpu)  ! Meridional index for the slab

  ! Step 2a: Southward information passing (phase a)
  ! -----------------------------------------------
  ! Processes with jid = 1, 3, 5, ... send out information
  southward_a_send: if ( mod( jid, 2 ) == 1 ) then

    ! Figure out targetted process
    ! Note that because nmcpu = 1, sid = my_proc_id
    target_id = (jid-1) * nicpu + iid 

    ! Pack info into buffer
    buffer_sn = x( :,1,:,:,: )
    
    ! Send information to target process
    call MPI_SEND(buffer_sn, ni*nk*nv*nm, MPI_REAL, target_id, my_proc_id, comm, ierr)

  endif southward_a_send


  ! Processes with jid = 0, 2, 4, .... receive info
  southward_a_recv: if ( mod( jid, 2 ) == 0 ) then

    ! Figure out the process that is sending info
    target_id = (jid+1) * nicpu + iid 

    ! Receive information
    if ( jid < njcpu -1 ) then
      call MPI_RECV( buffer_sn, ni*nk*nv*nm, MPI_REAL, target_id, target_id, comm, &
                     MPI_STATUS_IGNORE, ierr)
    
      ! Put info into boundary variable
      x_sn_bdy(2,:,:,:,:) = buffer_sn
    endif

  endif southward_a_recv


  ! Wait for all processes to reach this point
  call MPI_Barrier( comm, ierr )


  ! Step 2b: Southward information passing (phase b)
  ! -----------------------------------------------
  ! Processes with jid = 2, 4, 6 .... send out info to jid = 1, 3, 5
  southward_b_send: if ( mod(jid, 2) == 0 .and. jid >= 2) then

    ! Figure out the receiving process id
    target_id = (jid-1) * nicpu + iid 

    ! Pack info into buffer
    buffer_sn = x(:,1,:,:,:)

    ! Send info to target process
    call MPI_SEND(buffer_sn, ni*nk*nv*nm, MPI_REAL, target_id, my_proc_id, comm, ierr)

  endif southward_b_send 

  ! Processes with jid = 1, 3, 5, .... receive info
  southward_b_recv: if ( mod( jid, 2 ) == 1 ) then

    ! Figure out the process that is sending info
    target_id = (jid+1) * nicpu + iid

    ! Receive information
    if ( jid < njcpu -1 ) then
      call MPI_RECV( buffer_sn, ni*nk*nv*nm, MPI_REAL, target_id, target_id, comm, &
                     MPI_STATUS_IGNORE, ierr)

      ! Put info into boundary variable
      x_sn_bdy(2,:,:,:,:) = buffer_sn
    endif

  endif southward_b_recv

  call MPI_BARRIER( comm, ierr )


end subroutine communicate_boundary_southward
! ========================================================================================




! ========================================================================================
subroutine communicate_boundary_corners()

  ! Load derived types
  use namelist_define
  use wrf_tools
  use mpi_module
  use common_variables
  use boundary_variables

  implicit none
  integer :: sid, iid, jid ! ids to indicate slab position
  integer :: target_id    ! Process id of receiving process
  real, dimension( nicpu, njcpu, 2, 2, nk, nv, nm) :: buffer_send ! Buffer used to send corners
  real, dimension( nicpu, njcpu, 2, 2, nk, nv, nm) :: buffer_recv ! Buffer used to recv corners


  ! Step 1: Figure out process sid, gid, iid, jid 
  ! ---------------------------------------------
  sid=mod(my_proc_id,nicpu*njcpu)
  iid=mod(sid,nicpu)  ! Zonal index for the slab
  jid=int(sid/nicpu)  ! Meridional index for the slab


  ! Step 2: File all corner info into buffer
  ! ----------------------------------------
  buffer_send = 0.

  ! Pass info towards the northeast
  if ( ( iid < nicpu - 1 ) .and. ( jid < njcpu - 1 ) ) &
    buffer_send( iid+1, jid+1, 1,1,:,:,:) = x(ni,nj,:,:,:)

  ! Pass info towards the southeast
  if ( ( iid < nicpu - 1 ) .and. ( jid > 0 ) ) &
    buffer_send( iid+1, jid-1, 1,2,:,:,:) = x(ni,1,:,:,:)

  ! Pass info towards the northwest
  if ( ( iid > 0 ) .and. ( jid < njcpu - 1 ) ) &
    buffer_send( iid-1, jid+1, 2,1,:,:,:) = x(1,nj,:,:,:)

  ! Pass info towards the southwest
  if ( ( iid > 0 ) .and. ( jid > 0 ) ) &
    buffer_send( iid-1, jid-1, 2,2,:,:,:) = x(1,1,:,:,:)


  ! Step 3: Pass corner info across all procs
  ! -----------------------------------------
  buffer_recv = 0.
  call MPI_ALLREDUCE( buffer_send, buffer_recv, nicpu*njcpu*2*2*nk*nv*nm, &
                      MPI_REAL, MPI_SUM, comm, ierr )
  call MPI_BARRIER( comm, ierr )


  ! Step 4: Load corner info into boundary
  ! --------------------------------------
  x_crnrs(:,:,:,:,:) = buffer_recv(iid, jid, :,:,:,:,:)


end subroutine communicate_boundary_corners
! =====================================================================================

