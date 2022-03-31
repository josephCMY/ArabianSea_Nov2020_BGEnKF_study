! ================================================================================
! Module to perform Bi-Gaussian Ensemble Kalman Filter
! Original Publication:
! 1) Chan, M.-Y., Anderson J. L. and Chen X. (2020): An Efficient Bi-Gaussian
!    Ensemble Kalman Filter for Satellite Infrared Radiance Data Assimilation.
!    Monthly Weather Review. doi:10.1175/MWR-D-20-0142.1
!
! Notes
! 1) The cluster transfer process is now completely deterministic
!    Different from Chan et al (2020), which used a stochastic approach to
!    generate the resampling matrix.
! 2) BGEnKF is designed with state space in mind. Can be extended to joint space
!    in the future.
! 3) Unlike the PSU-EnKF, the variable x used in this program MUST BE the full
!    state, not the perturbation state!
! 4) For ease of coding, the joint-space GMM-EnKF redefines variable x so that
!    observations and aux variables are included.
! ================================================================================

module gauss_mix_model_enkf

  ! Importing useful modules
  ! ------------------------
  use constants
  use namelist_define
  use wrf_field_indices_define
  use obs_define
  use mpi_module
  !use netcdf
  use common_variables



  ! ==============================================================================

  ! Variable definitions
  ! --------------------
  implicit none


  ! Gaussian Mixture Model variables
  ! --------------------------------
  ! Integers to indicate which kernel each member belong to
  ! Has dimensions of (obs%num,numbers_en)
  integer :: cid
  integer, allocatable, dimension(:,:) :: kernel_id

  ! Parameters to include ya and aux into model slab variable x
  integer :: nv_increase

  ! Gaussian kernel parameters
  real, allocatable, dimension(:) :: prior_cluster_yvariance
  real, allocatable, dimension(:) :: prior_cluster_ymean
  real, allocatable, dimension(:) :: prior_weights, poste_weights
  real, allocatable, dimension(:,:,:) :: prior_x_mean, poste_x_mean

  ! Variables for sanity checking cluster transfer procedures
  real :: origin_xvar, origin_xavg      ! Pre-transfer values
  real :: finish_xvar, finish_xavg      ! Post-transfer values
  integer :: check_i, check_j, check_k  ! Slab indices to perform checks for


  integer, parameter :: min_mems_per_kernel = 5      ! Minimum number of members needed in each kernel
                                                     ! Otherwise, revert from GMM-EnKF to EnKF.



  ! Number of obs assimilated prior to the current batch
  integer :: batch_start_obs_num

  ! Auxiliary variable used for joint-space clustering.
  real,    allocatable, dimension(:,:      ) :: aux, aux_send

  ! Kernel IDs when assimilating a specific observation
  integer, allocatable, dimension(:) :: prior_id, poste_id

  ! Messenger variable to communicate kernel ids into the
  ! remove_overly_small_clusters subroutine
  integer, allocatable, dimension(:) :: msg_id

  ! Flags to indicate if an expanding cluster member was inflated
  integer, allocatable, dimension(:) :: inds_infl, inds_uninfl, inds_incoming, flag_infl

  ! Kalman gains and sqrt factors
  real, allocatable, dimension(:,:,:) :: km
  real :: alpha ! Square-root modification factor
  real, allocatable, dimension(:,:,:) :: mean_x_increment

  ! Variable to store backup copy of slab
  real, allocatable, dimension(:,:,:,:) :: backup_slab, ens_increment

  ! Variable to store backup copy of ensemble obs
  real, allocatable, dimension(:) :: backup_ya, ya_increment

  ! Variable to store backup copy of auxiliary variable
  real, allocatable, dimension(:) :: backup_aux, aux_increment

  ! Variables relating to expanding cluster
  integer :: n_shift, n_src_mems, prior_cluster_size
  real :: exp_infl_fac


  ! MPI and WRF file related variables
  ! ----------------------------------
  integer              :: sid,gid,iid,jid,g_comm,s_comm
  character (len=10)   :: wrf_file
  !type (proj_info)     :: proj
  !integer :: my_proc_id ! Remove this in the future


  ! Obs-related variables
  ! ---------------------
  character (len=10)   :: obstype
  integer              :: iunit, ounit  ! Used for writing fort.1000* files


  ! Variables for counting
  ! ----------------------
  integer   :: i, j, k, m, n, iob, iiob, nob, ie, ii, jj, kk, is, it
  integer   :: ig, iv, i1,j1, itot, iit


  ! Slab-related and update zone variables
  ! --------------------------------------
  integer   :: ist,ied,jst,jed,kst,ked,istart,iend,jstart,jend,uist,uied,ujst,ujed
  integer   :: sist,sied,sjst,sjed, sistart,siend,sjstart,sjend,sis,sie,sjs,sje


  ! Variables needed to do joint-space updates
  ! ------------------------------------------
  ! Index positions where the auxilaries and observables are slotted into the
  ! variable x. Needed for performing joint-space updates
  integer, allocatable, dimension(:,:) :: obs_mapping_array
  integer, allocatable, dimension(:,:) :: aux_mapping_array
  integer, dimension(3) :: roi_holder

  ! Misc variables
  ! --------------
  ! Other handy variable definitions
  ! fac = 1/(n-1) used for averaging in variance calculation
  ! loc_fac = localization factor
  ! ngx, ngz = roi in horizontall and vertical
  ! var,cov = error variance,covariance of something
  ! y_hxm = y-hxm (the innovation vector, mean), hxa is the perturbation (of ensemble members).
  real      :: fac,d,var,cov,y_hxm, loc_fac
  real      :: ir_bias, ir_bias_error
  real, allocatable, dimension(:) :: hxa
  integer   :: ngx, ngz
  integer, dimension(8)  :: values
  character (len=10) :: filename, date, time, zone, varname
  character (len=20) :: format1
  real ::  timer,t0,t_update_x, t_update_y
  real :: t1,t2, t_ab, t_c, t_d, t_e, t_f, t_g, t_h, t_buffer
  integer :: grid_id, fid, rcode, num_update_var, assimilated_obs_num, update_flag
  real      :: gaussdev, error, xb
  real, dimension(3)      :: center_xb
  ! state vectors (x): xm=ensemble mean; x=perturbation from ensemble mean
  ! _f=first guess; _t=truth
  ! ya=Hxa,Hxm: the state vector translated to observation space. ya(number of obs, number of members + mean)
  ! km        : kalman gain in updating x
  !real, dimension(ni,nj,nk,nv)                   :: x_t, std_x,std_xf
  !real, dimension(3,3,nk,nv,nm,nicpu*njcpu) :: xobsend, xob
  !real, dimension(3,3,nk,nv) :: tmp, tmpsend
  !real, dimension(obs%num,numbers_en+1) :: ya, yasend, yf
  real  :: ya_mean
  real :: dx
  ! for satellite radiance
  integer :: iob_radmin,iob_radmax
  !real, dimension(obs%num) :: ya_ca
  !real, dimension(obs%num,nm) :: yasend_tb
  !real, dimension(ni,nj,nk)     :: xq_n,xq_p
  !real, dimension(ni,nj,nk,nm)  :: xq_nsend,xq_psend
  ! Empirical Localization Functions
  integer :: n_r=101
  integer :: ich, iob_satch,num_satch
  ! parameter estimation
  integer :: iost,update_file
  ! Variables for slab tagging
  integer :: flag_slab_update, flag_run_xvar
  integer, allocatable, dimension(:,:) :: slab_mask, obs_mask

  ! ==============================================================================




  ! Subroutine section to follow
  CONTAINS



  ! ==============================================================================
  ! Main subroutine to perform GMM-EnKF
  subroutine gmm_enkf(wrf_file0, g_comm0, s_comm0, gid0, sid0, iid0, jid0 )!, proj0)

    ! Arguments to ingest during subroutine call
    implicit none
    integer, intent(in) :: sid0, gid0, iid0, jid0, g_comm0, s_comm0
    character (len=10), intent(in) :: wrf_file0
    !type( proj_info ) :: proj0


    ! ------------------------------------------------------
    ! GMM-EnKF preamble
    ! ------------------------------------------------------

    ! Pass arguments into the module variables
    wrf_file = wrf_file0; !proj = proj0
    g_comm = g_comm0; s_comm = s_comm0
    gid = gid0; sid = sid0; iid = iid0; jid = jid0

    ! Pretty print out to indicate start of GMM-EnKF
    !read(wrf_file(6:10),'(i5)')iunit
    if ( my_proc_id==0 ) then
       write(*, *)'   '
       write(*, *)'---------------------------------------------------'
       write(*, *)'.... Begin GMM-EnKF Assimilation ....'
       write(*, *)'   '
    endif

    ! Initialize counter for the number of observations assimilated
    assimilated_obs_num = 0

    ! Immediate termination of subroutine if no obs are present
    if ( obs%num == 0 ) then
      if ( my_proc_id == 0) write(*,*) 'No obs found. Exiting gmm_enkf subroutine.'
      return
    endif


    ! Quick check on IR biases
    if ( trim(obs%type(1)) == 'Radiance' .and. my_proc_id == 0) then
      ir_bias = sum( ya(:,numbers_en+1) - obs%dat(:) ) / obs%num
      ir_bias_error = sum( (ya(:,numbers_en+1) - obs%dat(:))**2 )/obs%num &
                      - ir_bias**2
      ir_bias_error = sqrt( ir_bias_error ) / sqrt( obs%num * 1. )

      write(*,*) 'Mean innovation of', obs%num,' IR obs:' , &
                  ir_bias, '+-', ir_bias_error
    endif

    ! Allocate arrays for GMM-EnKF
    call allocate_gmm_enkf_arrays()

    ! Compute auxiliary variables for joint-space GMM-EnKF
    nv_increase=0
    if ( use_jointspace_gmm_enkf ) then
      call construct_obs_aux_to_x_mapping_array()
      call map_modelspace_x_to_jointspace_x()      
      call MPI_BARRIER( comm, ierr )
    endif
    if (my_proc_id == 0) write(*,*)size(x)

    ! Load slab and domain info, and define mask that marks slab zone
    call get_slab_info()

    ! Initialize Kalman gain arrays
    km =0.

    ! Initialize timer variables
    !timer=MPI_Wtime()
    t_update_x = 0.; t_update_y = 0.
    t_ab = 0.; t_c = 0.; t_d = 0.; t_d = 0.; t_e = 0.;
    t_f = 0.; t_g = 0.; t_h = 0.


    ! Loop over every obs in the obs object
    if(my_proc_id==0) write(*,*) 'Assimilating obs...'




    ! -------------------------------------------------------------
    ! Serially assimilate observations via GMM-EnKF
    ! -------------------------------------------------------------
    obs_assimilate_cycle : do it = 1,obs%num

      ! Load obs info and compute innovation
      call load_obs_and_innovation()

      ! Quality control obs (current approach can be refined in the future)
      call crude_obs_quality_control()

      ! If QC rejected obs, skip over to the next observation
      if ( obs%kick_flag(iob) == 1 ) cycle obs_assimilate_cycle

      ! Perform Adaptive Observation Error Inflation (if needed)
      if (use_aoei) call evoke_aoei()

      ! If using joint-space GMM-EnKF, use auxiliary variable to tag kernel
      jointspace_kernel_tagging: if ( use_jointspace_gmm_enkf ) then

        ! Tag using auxiliary
        kernel_id=1
        where (aux > 1e-3) kernel_id =2
        if (my_proc_id ==0) &
          write(*,'(2(a,x,i4,3x))') 'kernel 1 count:',count( kernel_id(iob,:) == 1 ), &
                                    'kernel 2 count:',count( kernel_id(iob,:) == 2 )

        ! Special measure to protect against weight sampling errors
        if ( ( count( kernel_id(iob,:) == 1 ) < min_mems_per_kernel ) &
               .or. ( count( kernel_id(iob,:) == 2 ) < min_mems_per_kernel ) )  kernel_id(iob,:) = 1

      endif jointspace_kernel_tagging

      ! Check if the obs will update the slab
      call check_obs_slab_overlap()

      ! Estimate prior kernel parameters (weight, mean and covariances)
      call estimate_prior_kernel_parameters()

      ! Compute prior and posterior kernel weights for this observation
      call compute_posterior_weights()

      
      ! Heuristic check to see if expanding cluster has enough members
      limit_growing_cluster_sampling_errors: do cid = 1, max_kernel_num
        if (  ( prior_weights(cid) < poste_weights(cid) ) &
              .and. ( count(kernel_id(iob,:) == cid) < min_expanding_kernel_prior_size) &
           )  then
          kernel_id(iob,:) = 1
          prior_weights = 0.; poste_weights = 0.
          prior_weights(1) = 1.; poste_weights(1) = 1.
          if (my_proc_id == 0) write(*,*)'Heuristic check triggered.'
          
          ! Re-estimate prior kernel parameters (weight, mean and covariances)
          call estimate_prior_kernel_parameters()

          !if (my_proc_id == 0) write(*,'(a,x,2f4.2)')'New prior weights', prior_weights
          !if (my_proc_id == 0) write(*,'(a,x,2(f4.2,x))')'New poste weights', poste_weights

        endif
      enddo limit_growing_cluster_sampling_errors


      ! If kernel weight updates are unphysical, default from BGEnKF to EnKF 
      ! Case 1: clear obs DA caused clear cluster shrinkage
      ! Case 2: cloudy obs DA caused cloud cluster shrinkage
      default_BGEnKF_to_EnKF: if ( trim( obs%type(iob) ) == 'Radiance' ) then

        ! Case 1: clear obs DA caused clear cluster shrinkage
        default_case1: if ( (obs%dat(iob) > 290. ) .and. ( prior_weights(1) > poste_weights(1) ) ) then
          if(my_proc_id==0) write(*,*) 'Obs',iob,'clear obs caused clear cluster shrinkage'
          kernel_id(iob,:) = 1
          prior_weights = 0.; poste_weights = 0.
          prior_weights(1) = 1.; poste_weights(1) = 1.
          ! Re-estimate prior kernel parameters (weight, mean and covariances)
          call estimate_prior_kernel_parameters()

        endif default_case1
        
       ! Case 2: cloudy obs DA caused cloud cluster shrinkage
       default_case2: if ( (obs%dat(iob) < 280. ) .and. ( prior_weights(2) > poste_weights(2) ) ) then
          if(my_proc_id==0) write(*,*) 'Obs',iob,'cloudy obs caused cloudy cluster shrinkage'
          kernel_id(iob,:) = 1
          prior_weights = 0.; poste_weights = 0.
          prior_weights(1) = 1.; poste_weights(1) = 1.
          ! Re-estimate prior kernel parameters (weight, mean and covariances)
          call estimate_prior_kernel_parameters()

        endif default_case2

      endif default_BGEnKF_to_EnKF

      !call MPI_BARRIER(comm, ierr)

      ! Determine how members will migrate between clusters
      call setup_member_migration_paths()

      ! Print out the migration pattern
      if ( my_proc_id == 0 .and. print_detail  > 10 ) then
        write(*,*) 'Prior weights', prior_weights
        write(*,*) 'Prior means', prior_cluster_ymean
        write(*,*) 'Prior variance', prior_cluster_yvariance
        write(*,*) 'Poste weights', poste_weights
        write(*,*) 'Printing out the migration paths'
        do ie = 1, numbers_en
          write(*,'(a,x,i4,x,2(a,x,f6.1,x),2(a,x,i4,x))') &
            'member',ie, 'ya=',ya(iob, ie),'yo=',obs%dat(iob), &
            'prior cluster', prior_id(ie),'poste cluster', poste_id(ie)
        enddo
      endif




      ! -------------------------------------------------------------------------
      ! Applying GMM-EnKF updates to the slab held by process
      ! -------------------------------------------------------------------------
      ! Iterate thru model field variables
      update_x_var : do m = 1, nv+nv_increase

        ! Check if variable needs updating. If not, skip to next variable
        check_updatevar_list: if (m<=nv) then
          varname=enkfvar(m)
          update_flag = 0
          do iv = 1, num_update_var
            if ( varname .eq. updatevar(iv) ) update_flag = 1
          enddo
          if ( update_flag==0 ) cycle update_x_var
        endif check_updatevar_list


        ! Set update zone to entire slab when handling observables and auxiliary
        ! ----------------------------------------------------------------------
        change_aux_obs_roi: if (m>nv) then
          ! Saving a copy of the actual observation roi
          roi_holder(:) = obs%roi(iob,:)
          ! Hijack the obs roi
          obs%roi(iob,:) = 999999
        endif change_aux_obs_roi

        ! -----------------------------------------------------------------------
        ! Updating model variable fields using observation iob
        !
        ! Procedure:
        ! ----------
        ! 1) Make a backup copy of model field slab (backup_slab)
        !
        ! 2) Kalman filter step on individual cluster.
        !    For each cluster:
        !    a. Compute localized Kalman gain and EnSRF sqrt factor
        !    b. Apply localized Kalman Filter update on all cluster members
        !
        ! 3) For each shrinking cluster, recenter the non-outgoing members on
        !    posterior kernel mean.
        !
        ! 4) Generate new states in each expanding clusters:
        !    a. Generate the resampling matrix E described in Chan et al (2020)
        !    b. Overwrite incoming members' states with new resampled states
        !    c. Inflate pre-existing states (if needed)
        !
        ! 5) Applying localization to model-space updates
        !    a. Compute update increment using x - backup_slab
        !    b. Localize update increment
        !    c. Apply localized update increment:
        !       x = backup_slab + loc * ( x - backup_slab )
        ! -----------------------------------------------------------------------

        ! Step 1: Make backup slab of all members
        ! ---------------------------------------
        backup_slab = x(:,:,:,m,1:numbers_en)


        ! Step 2: Apply Kalman filter on each cluster
        ! -------------------------------------------
        update_cluster_x_var: do cid = 1, max_kernel_num

          ! Identify observation update zone within slab
          ! In the future, each kernel can potentially have diff ROI
          call identify_obs_update_zone()

          ! Skip update to x update zone has no overlap with slab
          if ( flag_slab_update == 0 ) cycle update_cluster_x_var

          ! Compute prior slab mean, Kalman gain and sqrt factor
          call compute_kalman_update_parameters()

          ! Perform Kalman filter update on cluster members
          call apply_enkf_update()

        enddo update_cluster_x_var



        ! Step 3: Recenter shrinking cluster means
        ! ---------------------------------------
        recenter_shrinking_cluster_means: do cid = 1, max_kernel_num

          ! Skip over non-shrinking clusters
          if ( prior_weights(cid) <= poste_weights(cid) ) &
            cycle recenter_shrinking_cluster_means

          ! Recenter residual members in shrinking cluster
          call recenter_residual_members()

        enddo recenter_shrinking_cluster_means



        ! Step 4: Apply expansion to expanding cluster
        ! --------------------------------------------
        grow_expanding_clusters: do cid = 1, max_kernel_num

          ! Skip over non-expanding clusters
          if ( prior_weights(cid) >= poste_weights(cid) ) &
            cycle grow_expanding_clusters

          ! Generate parameters relating to resampling
          call expansion_resampling_params()

          ! Preparing pre-expansion statistics. Will be later used for sanity
          ! checking the expansion procedure.
          call sanity_check_cluster_expansion_preamble()

          ! Estimate expanding cluster posterior mean state
          call estimate_cluster_xa_mean()

          ! Convert all cluster states to perturbations
          call cluster_states_to_perturbations()

          ! Inflate (prior_cluster_size - n_src_mems) perturbations
          call inflate_cluster_perts()

          ! Construct new perturbations using matrix E
          call construct_new_perturbations()

          ! Apply mean to all perturbations in cluster
          call poste_cluster_perturbations_to_states()

          ! Sanity check the cluster expansion procedure
          call sanity_check_cluster_expansion()

        enddo grow_expanding_clusters


        ! Special sanity check subroutine
        ! ------------------------------
        if ( nint( obs%position(iob,1) ) >= istart               &
             .and. nint( obs%position(iob,1) ) <= iend           &
             .and. nint( obs%position(iob,2) ) >= jstart         &
             .and. nint( obs%position(iob,2) ) <= jend           &
             .and. count( prior_weights > 0. ) == max_kernel_num &
             .and. count( poste_weights > 0. ) == max_kernel_num &
             .and. print_detail > 0                              &
             .and. trim(varname) == 'QVAPOR'                     ) then
           call gmm_enkf_check()
         endif


        ! Undo ROI hijacking for observables and auxiliary variables
        ! ----------------------------------------------------------
        if (m>nv) obs%roi(iob,:) = roi_holder(:)



        ! Step 5: Apply localization
        ! --------------------------
        ! Compute increment for each member
        ens_increment(uist:uied, ujst:ujed, kst:ked,:) &
          = x(uist:uied, ujst:ujed, kst:ked, m,1:numbers_en) &
            - backup_slab(uist:uied, ujst:ujed, kst:ked,:)


        ! Apply localization factor on model-space increment
        ! --------------------------------------------------
        localize_modelspace_increment: if (m<=nv) then
          do k = kst, ked
            do j = ujst, ujed
              do i = uist, uied
                ! Compute Gaspari-Cohn localization factor
                call corr( (  i + istart - 1 ) - obs%position(iob,1), &
                           (  j + jstart - 1 ) - obs%position(iob,2), &
                           k - obs%position(iob,3), ngx, ngz, loc_fac )
                ! Apply GC localization factor
                ens_increment(i,j,k,:) = ens_increment(i,j,k,:) * loc_fac
              enddo
            enddo
          enddo
        endif localize_modelspace_increment



        ! Apply localization factor on observable increments
        ! --------------------------------------------------
        localize_ya_increment: if ( (m>nv) .and. (m<=nv+nv_increase/2) ) then
          localize_ya_increment_loop: do iit = 1, obs%num
            iiob = ind( iit )

            ! Skip over observables that have different m.
            ! This is to prevent performing localization more than once.
            if ( obs_mapping_array(iiob,4) /= m ) &
              cycle localize_ya_increment_loop

            ! skip updates to observables outside the update zone (not slab)
            if ( obs%position(iiob,1)<ist .or. obs%position(iiob,1)>ied .or. &
                 obs%position(iiob,2)<jst .or. obs%position(iiob,2)>jed .or. &
                 obs%position(iiob,3)<kst .or. obs%position(iiob,3)>ked ) &
               cycle localize_ya_increment_loop


            ! Compute Gaspari-Cohn correlation factor
            call corr( real(obs%position(iiob,1)-obs%position(iob,1)), &
                       real(obs%position(iiob,2)-obs%position(iob,2)), &
                       real(obs%position(iiob,3)-obs%position(iob,3)), &
                       obs%roi(iob,1), obs%roi(iob,2), loc_fac )

            ! Apply localization factor
            i = obs_mapping_array(iiob,1); j = obs_mapping_array(iiob,2)
            k = obs_mapping_array(iiob,3)
            ens_increment(i,j,k,:) = ens_increment(i,j,k,:) * loc_fac

          enddo localize_ya_increment_loop
        endif localize_ya_increment



        ! Apply localization factor on auxiliary increments
        ! --------------------------------------------------
        localize_aux_increment: if ( (m>nv+nv_increase/2) ) then
          localize_aux_increment_loop: do iit = 1, obs%num
            iiob = ind( iit )

            ! Skip over observables that have different m.
            ! This is to prevent performing localization more than once.
            if ( aux_mapping_array(iiob,4) /= m ) &
              cycle localize_aux_increment_loop

            ! skip updates to observables outside the update zone (not slab)
            if ( obs%position(iiob,1)<ist .or. obs%position(iiob,1)>ied .or. &
                 obs%position(iiob,2)<jst .or. obs%position(iiob,2)>jed .or. &
                 obs%position(iiob,3)<kst .or. obs%position(iiob,3)>ked ) &
               cycle localize_aux_increment_loop

            ! Compute Gaspari-Cohn correlation factor
            call corr( real(obs%position(iiob,1)-obs%position(iob,1)), &
                       real(obs%position(iiob,2)-obs%position(iob,2)), &
                       real(obs%position(iiob,3)-obs%position(iob,3)), &
                       obs%roi(iob,1), obs%roi(iob,2), loc_fac )

            ! Apply localization factor
            i = aux_mapping_array(iiob,1); j = aux_mapping_array(iiob,2)
            k = aux_mapping_array(iiob,3)
            ens_increment(i,j,k,:) = ens_increment(i,j,k,:) * loc_fac

          enddo localize_aux_increment_loop
        endif localize_aux_increment


        ! Apply the localized increment
        x(uist:uied, ujst:ujed, kst:ked, m, 1:numbers_en) &
          = backup_slab(uist:uied, ujst:ujed, kst:ked, :) &
            + ens_increment(uist:uied, ujst:ujed, kst:ked, :)

        ! Recalculate mean state
        x(uist:uied, ujst:ujed, kst:ked, m, numbers_en+1) &
        = sum( x(uist:uied, ujst:ujed, kst:ked, m, 1:numbers_en), 4) / numbers_en

      enddo update_x_var


      ! Refresh observables and auxiliary variables
      ! -------------------------------------------
      extract_obs_aux_from_jointspace: if (use_jointspace_gmm_enkf) then
        do iiob = 1, obs%num
          ya(iiob, :) = x( obs_mapping_array(iiob,1), obs_mapping_array(iiob,2), &
                           obs_mapping_array(iiob,3), obs_mapping_array(iiob,4), : ) *1.
          aux(iiob, :)= x( aux_mapping_array(iiob,1), aux_mapping_array(iiob,2), &
                           aux_mapping_array(iiob,3), aux_mapping_array(iiob,4), : )*1.
        enddo
      endif extract_obs_aux_from_jointspace
      !ya(:,numbers_en+1) = sum( ya(:,1:numbers_en), 2)/numbers_en
      !aux(:,numbers_en+1) = sum( aux(:,1:numbers_en), 2)/numbers_en
   

    enddo obs_assimilate_cycle


    ! Revert from joint-space to model-space
    ! --------------------------------------
    if ( use_jointspace_gmm_enkf ) call map_jointspace_x_to_modelspace_x()


    if (my_proc_id == 0 .and. print_detail > 10) then
      do cid = 1, 3
        write(*,'(2(a,x,i5,x))') 'targetted cluster',cid, 'size:', &
                                 count(poste_id == cid)
      enddo
    endif

    call MPI_BARRIER( comm, ierr )
    ! Deallocate dynamically allocated arrays used by gmm-enkf
    call deallocate_gmm_enkf_arrays()

    end subroutine gmm_enkf
    ! =================================================================================




    ! =========================================================================
    ! Subroutine to generate the index positions of obs and aux on the
    ! extended x variables.
    subroutine construct_obs_aux_to_x_mapping_array()

      implicit none
      integer :: vv

      ! Figure out how many "faux" variables we need to add to x
      nv_increase =  ceiling( (obs%num)*1. / (ni*nj*nk) ) * 2

      ! Set up mapping rules for observables
      iob = 1
      do vv = nv+1, nv + nv_increase/2
        do kk = 1, nk
          do jj = 1, nj
            do ii = 1, ni

              ! Store mapping rules
              obs_mapping_array(iob,1) = ii
              obs_mapping_array(iob,2) = jj
              obs_mapping_array(iob,3) = kk
              obs_mapping_array(iob,4) = vv

              ! Increment observation being addressed
              iob = iob + 1

              ! Terminate loop if all obs have been assigned positions on x
              if ( iob > obs%num ) exit
            enddo

            ! Terminate loop if all obs have been assigned positions on x
            if ( iob > obs%num ) exit
          enddo

          ! Terminate loop if all obs have been assigned positions on x
          if ( iob > obs%num ) exit
        enddo

        ! Terminate loop if all obs have been assigned positions on x
        if ( iob > obs%num ) exit
      enddo


      ! Set up mapping rules for auxiliary variables
      iob = 1
      do vv = nv+nv_increase/2+1, nv+nv_increase
        do kk = 1, nk
          do jj = 1, nj
            do ii = 1, ni

              ! Store mapping rulez
              aux_mapping_array(iob,1) = ii
              aux_mapping_array(iob,2) = jj
              aux_mapping_array(iob,3) = kk
              aux_mapping_array(iob,4) = vv

              ! Increment observation being addressed
              iob = iob + 1

              ! Terminate loop if all obs have been assigned positions on x
              if ( iob > obs%num ) exit
            enddo

            ! Terminate loop if all obs have been assigned positions on x
            if ( iob > obs%num ) exit
          enddo

          ! Terminate loop if all obs have been assigned positions on x
          if ( iob > obs%num ) exit
        enddo

        ! Terminate loop if all obs have been assigned positions on x
        if ( iob > obs%num ) exit
      enddo


      ! Sanity check: were all obs positions assigned?
      if (iob <= obs%num) then
        write(*,'(a,x,i,a,x,i,x,a)') &
          'BUG: Only', iob, "out of", obs%num, &
          "observations were assigned mapping rules for jointspace filtering"
        write(*,*)'Terminating process', my_proc_id
        call exit()
      endif
      call MPI_BARRIER( comm, ierr)

      ! Print statement to indicate all is well
      if (my_proc_id == 0) then
        write(*,*) "Mapping rules for jointpsace filtering assigned"
      endif

    end subroutine construct_obs_aux_to_x_mapping_array
    ! =========================================================================




    ! =================================================================================
    ! Subroutine to redefine x s.t. it contains the auxiliary and observables
    subroutine map_modelspace_x_to_jointspace_x()

      implicit none

      real, allocatable, dimension(:,:,:,:,:) :: x_hold

      ! Allocate temporary holding array and store model space elements
      allocate( x_hold ( ni, nj, nk, nv, nm) )
      x_hold = x(:,:,:,:,:) * 1.


      ! Construct the auxiliary variable
      if (my_proc_id == 0) write(*,*) 'entering compute_auxiliary_variable'
      call compute_auxiliary_variable( iid, jid )

      ! Deallocate x
      deallocate(x)

      ! Reallocate x to include aux and obs variables
      allocate( x( ni, nj, nk, nv + nv_increase, nm ) )

      ! Insert placeholder
      x(:,:,:,:,:) = -888888

      ! Stick in model variables
      x(:,:,:,1:nv,:) = x_hold *1.

      ! Stick in observables using the mapping rules
      iob=1
      do iob = 1, obs%num
          x( obs_mapping_array(iob,1), obs_mapping_array(iob,2), &
             obs_mapping_array(iob,3), obs_mapping_array(iob,4), : ) &
          = ya(iob,:)
      enddo

      ! Stick in auxiliary using the mapping rules
      do iob = 1, obs%num
          x( aux_mapping_array(iob,1), aux_mapping_array(iob,2), &
             aux_mapping_array(iob,3), aux_mapping_array(iob,4), : ) &
          = aux(iob,:)
      enddo

      ! Freeing up memory
      deallocate(x_hold)

    end subroutine map_modelspace_x_to_jointspace_x
    ! =================================================================================




    ! =================================================================================
    ! Subroutine to remove observables and auxiliary from x
    subroutine map_jointspace_x_to_modelspace_x()

      implicit none
      real, allocatable, dimension(:,:,:,:,:) :: x_hold

      ! Allocate temporary holding array and store model space elements
      allocate( x_hold ( ni, nj, nk, nv, nm) )
      x_hold = x(:,:,:,1:nv,:)

      ! Reallocate x and put stuff back in
      deallocate(x); allocate( x(ni, nj, nk, nv, nm) )
      x(:,:,:,:,:) = x_hold(:,:,:,:,:)

      ! Free up memory
      deallocate( x_hold )

    end subroutine map_jointspace_x_to_modelspace_x
    ! =================================================================================




    ! =================================================================================
    ! Subroutine to dynamically allocate arrays for gmm enkf
    subroutine allocate_gmm_enkf_arrays()

      implicit none

      ! Gaussian kernel parameters
      allocate( prior_cluster_yvariance    ( max_kernel_num         ) )
      allocate( prior_cluster_ymean        ( max_kernel_num         ) )
      allocate( prior_weights              ( max_kernel_num         ) )
      allocate( poste_weights              ( max_kernel_num         ) )
      allocate( prior_x_mean               ( ni, nj, nk             ) )
      allocate( poste_x_mean               ( ni, nj, nk             ) )

      ! Kernel IDs
      !allocate( kernel_id                  ( obs%num, numbers_en    ) )
      allocate( prior_id                   ( numbers_en             ) )
      allocate( poste_id                   ( numbers_en             ) )

      ! Flags used to handle resampling of expanding cluster
      allocate( inds_infl                  ( numbers_en             ) )
      allocate( inds_uninfl                ( numbers_en             ) )
      allocate( inds_incoming              ( numbers_en             ) )
      allocate( flag_infl                  ( numbers_en             ) )

      ! Kalman filter update related arrays
      allocate( km                         ( ni, nj, nk             ) )
      allocate( mean_x_increment           ( ni, nj, nk             ) )


      ! Arrays used for localization
      allocate( backup_slab                ( ni, nj, nk, numbers_en ) )
      allocate( ens_increment              ( ni, nj, nk, numbers_en ) )
      allocate( backup_ya                  ( numbers_en             ) )
      allocate( ya_increment               ( numbers_en             ) )
      allocate( backup_aux                 ( numbers_en             ) )
      allocate( aux_increment              ( numbers_en             ) )

      ! Array to hold full ens obs perturbation
      allocate( hxa                        ( numbers_en             ) )

      ! Arrays to handle the masks
      allocate( slab_mask                  ( ix+1, jx+1             ) )
      allocate( obs_mask                   ( ix+1, jx+1             ) )

      ! Arraya to handle joint-space updates
      if ( use_jointspace_gmm_enkf ) then
        allocate( obs_mapping_array        ( obs%num, 4            ) )
        allocate( aux_mapping_array        ( obs%num, 4            ) )
      endif


    end subroutine allocate_gmm_enkf_arrays




    ! =================================================================================
    ! Subroutine to deallocate arrays for gmm enkf
    subroutine deallocate_gmm_enkf_arrays()

      implicit none

      ! Gaussian kernel parameters
      deallocate( prior_cluster_yvariance  )
      deallocate( prior_cluster_ymean      )
      deallocate( prior_weights            )
      deallocate( poste_weights            )
      deallocate( prior_x_mean             )
      deallocate( poste_x_mean             )

      ! Kernel IDs
      deallocate( prior_id                 )
      deallocate( poste_id                 )

      ! Flags used to handle resampling of expanding cluster
      deallocate( inds_infl                )
      deallocate( inds_uninfl              )
      deallocate( inds_incoming            )
      deallocate( flag_infl                )

      ! Kalman filter update related arrays
      deallocate( km                       )
      deallocate( mean_x_increment         )

      ! Arrays used for localization
      deallocate( backup_slab              )
      deallocate( ens_increment            )
      deallocate( backup_ya                )
      deallocate( ya_increment             )
      deallocate( backup_aux               )
      deallocate( aux_increment            )

      ! Array to hold full ens obs perturbation
      deallocate( hxa                      )

      ! Arrays to handle the masks
      deallocate( slab_mask                )
      deallocate( obs_mask                 )

      ! Arraya to handle joint-space updates
      if ( use_jointspace_gmm_enkf ) then
        deallocate( obs_mapping_array  )
        deallocate( aux_mapping_array  )
      endif


  end subroutine deallocate_gmm_enkf_arrays



  ! ==============================================================================
  ! Subroutine to obtain basic information about the slab held by process, and the
  ! domain itself. Also defines the mask that marks out the slab.
  subroutine get_slab_info()

    implicit none

    ! Get domain info
    !call open_file(wrf_file, nf_nowrite, fid)
    !rcode = nf_get_att_int(fid, nf_global, 'GRID_ID', grid_id)
    !rcode = nf_get_att_real(fid, nf_global, 'DX', dx)
    !dx=dx/1000.
    num_update_var = 0
    do m = 1, 20
      if ( len_trim(updatevar(m))>=1 ) num_update_var=num_update_var+1
    enddo
    !call close_file(fid)


    !subdomain start and end indices in full domain
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


    ! Define mask to determine slab location
    slab_mask = 0
    slab_mask( istart:iend, jstart:jend  ) =1


  end subroutine get_slab_info





  ! =============================================================================
  ! Subroutine to get obs info and compute innovation
  subroutine load_obs_and_innovation()

    implicit none

    ! Get basic information about selected observation
    iob=ind(it)
    obs%kick_flag(iob)=0
    obstype = obs%type(iob)
    error = obs%err(iob)

    ! Print out information about the observation
    y_hxm = obs%dat(iob) - ya(iob,numbers_en+1)
    if ( my_proc_id==0 ) write(*,'(/,a,i6,a,f10.2,a,f10.2,a,f10.4,a,f8.2,a,i4,a,i4,a,i4,a,i2)') &
                               'No.',iob+batch_start_obs_num, &
                               ' '//obstype//' =', obs%dat(iob),  &
                               ' ya=', ya(iob,numbers_en+1), ' y-ya=', y_hxm, &
                               ' err=',error,' hroi=',obs%roi(iob,1),'(', &
                               obs%roi(iob,3),') vroi=',obs%roi(iob,2),'at z', &
                               int(obs%position(iob,3))

    ! Load info about localization
    ngx = obs%roi(iob,1)
    ngz = obs%roi(iob,2)


  end subroutine load_obs_and_innovation
  ! =============================================================================





  ! =============================================================================
  ! Subroutine to do quality control on obs.
  ! Note that this QC is very very crude.
  subroutine crude_obs_quality_control()

    implicit none

    ! Ultra crude quality control
    if( abs(y_hxm)>(error*5.) .and. &
        .not. ( obstype=='min_slp   ' .or. obstype=='longitude ' .or. &
                obstype=='latitude  ' .or. obstype=='slp       '&
                .or. obstype=='Radiance  ' .or. obstype=='Microwave ' &
              ) &
       ) then
       !.or. obstype=='Radiance  ' .or. obstype(1:9)=='Psounding') ) then
      if ( my_proc_id==0 ) write(*,*)' ...kicked off for large error'
      obs%kick_flag(iob)=1
      return
    endif

    !! If AOEI is disabled, then we will also QC radiance obs
    !if ( (abs(y_hxm) > error*5  .and. trim(obstype) == 'Radiance' ) &
    !     .and. use_aoei == .false. ) then
    !  if ( my_proc_id==0)  write(*,*)' ...kicked off for large error'
    !  obs%kick_flag(iob)=1
    !  return
    !endif


    ! Special handling for Microwave observaitons
    if ( obstype=='Microwave ' .and. (ya(iob,numbers_en+1) .NE. ya(iob,numbers_en+1)) ) then
      if ( my_proc_id==0 ) write(*,*)' ...kicked off for NaN'
      obs%kick_flag(iob)=1
      return
    endif

    ! If quality control passes, increment number of assimilated obs
    assimilated_obs_num=assimilated_obs_num+1

  end subroutine crude_obs_quality_control
  ! ==============================================================================





  ! =============================================================================
  ! Subroutine to perform Adaptive Observation Error Inflation
  ! Publication: Minamide and Zhang (2017)
  subroutine evoke_aoei()

    implicit none

    ! Compute ensemble variance of simulated obs (used for AOEI)
    var = 0.
    do ie = 1, numbers_en
       hxa(ie) = ya(iob,ie) - ya(iob,numbers_en+1)
       var     = var + hxa(ie)*hxa(ie)
    enddo
    fac  = 1./real(numbers_en-1)


    ! Apply Minamide and Zhang (2017) Adaptive Obs Error Inflation
    if ( ( trim( obstype ) == 'Radiance'  .and. use_aoei ) .or. &
         ( trim( obstype ) == 'Microwave' .and. aoei_microwave ) ) then
      error = max( error*error,  y_hxm*y_hxm - fac*var )
      error = sqrt(error)
      if ( my_proc_id == 0 .and. error > obs%err(iob) ) &
          write(*,*) 'observation-error inflated to ',error
    endif

  end subroutine evoke_aoei
  ! =============================================================================




  ! =============================================================================
  ! Subroutine to check if obs updates the process' slab
  ! Uses the slab_mask defined earlier in get_slab_info to check
  subroutine check_obs_slab_overlap()

    implicit none

    ! Define zone spanned by selected observation
    ngx = max( obs%roi( iob, 1), obs%roi(iob,3) )
    ist = max( update_is, int(obs%position(iob,1))-ngx )
    ied = min( update_ie, int(obs%position(iob,1))+ngx )
    jst = max( update_js, int(obs%position(iob,2))-ngx )
    jed = min( update_je, int(obs%position(iob,2))+ngx )

    ! Check if obs update zone overlaps with slab using masks
    flag_slab_update = 0
    obs_mask = 0
    obs_mask(ist:ied,jst:jed) = 1
    flag_run_xvar = sum( obs_mask(ist:ied,jst:jed) *slab_mask(ist:ied, jst:jed) )

  end subroutine check_obs_slab_overlap
  ! =============================================================================





  ! =============================================================================
  ! Subroutine to estimate the prior kernels' weight, mean and covariance
  subroutine estimate_prior_kernel_parameters()

    implicit none
    ! Variable definitions
    integer, dimension( max_kernel_num ) :: kernel_counts
    integer :: kid

    ! Initialize values
    prior_cluster_yvariance = 0.; prior_cluster_ymean = 0.; kernel_counts = 0;

    ! Compute prior kernel mean, variance and weight in obs space
    compute_kernel_stats_loop: do kid = 1, max_kernel_num

      ! Compute kernel mean
      compute_prior_cluster_ymean: do n = 1, numbers_en
        if( kernel_id(iob,n) == kid ) then
          prior_cluster_ymean(kid) = prior_cluster_ymean(kid) + ya( iob, n )
          kernel_counts(kid) = kernel_counts(kid) + 1
        endif
      enddo compute_prior_cluster_ymean
      prior_cluster_ymean(kid) = prior_cluster_ymean(kid) / kernel_counts(kid)

      ! Compute kernel variance
      compute_prior_cluster_yvariance: do n = 1, numbers_en
        if( kernel_id(iob,n) == kid ) then
          prior_cluster_yvariance(kid) = prior_cluster_yvariance(kid) &
                                         + ( ya(iob, n) - prior_cluster_ymean(kid) )**2.
        endif
      enddo compute_prior_cluster_yvariance
      prior_cluster_yvariance(kid) = prior_cluster_yvariance(kid) / (kernel_counts(kid)-1.)

      ! Compute kernel weight
      prior_weights(kid) = (kernel_counts(kid)*1.) / (numbers_en*1.)

    enddo compute_kernel_stats_loop


  end subroutine estimate_prior_kernel_parameters
  ! =============================================================================




  ! =============================================================================
  ! Subroutine to determine prior and posterior weights of the kernels
  ! See Eqn (3) of Chan et al (2020)
  subroutine compute_posterior_weights()

    implicit none

    ! Allocate useful variables
    real, dimension( max_kernel_num ) :: poste_kernel_counts
    integer :: kid, biggest_cluster
    real :: expo, norm


    !if ( my_proc_id == 0) write(*,*) 'starting prior weights',prior_weights
    ! Evaluate posterior weight equation for GMM EnKF
    eval_poste_weight_loop: do kid = 1, max_kernel_num

      ! Compute exponent term
      expo = -0.5 * ( obs%dat(iob) - prior_cluster_ymean(kid) ) **2
      expo = expo / ( prior_cluster_yvariance(kid) + error**2 )

      ! Compute normalization of the distribution
      norm = 2. * pi * ( prior_cluster_yvariance(kid) + error**2 )
      norm = sqrt(norm)

      ! Evaluate posterior weight
      poste_weights(kid) = prior_weights(kid) * exp( expo ) / norm

      ! Exception for zero prior weight
      if (prior_weights(kid) == 0) poste_weights(kid) = 0


    enddo eval_poste_weight_loop
    !if ( my_proc_id == 0) write(*,*) 'starting poste weights',poste_weights


    ! If float precision can't handle the weight computation, then
    ! nullify weight updates and disable cluster transfer.
    if (sum( poste_weights ) < 1e-8) poste_weights = prior_weights


    ! Normalize posterior weight
    norm = sum( poste_weights )
    poste_weights = poste_weights/norm



    ! Ensuring that the posterior weight values are multiples of 1/numbers_en
    ! -----------------------------------------------------------------------
    ! Compute naive posterior kernel count and identify largest cluster
    poste_kernel_counts = poste_weights * numbers_en
    do kid = 1, max_kernel_num
      poste_kernel_counts(kid) = 1.*nint( poste_kernel_counts(kid) )
    enddo


    ! Identify biggest cluster
    biggest_cluster = 1
    do kid=1, max_kernel_num
      if ( poste_kernel_counts(kid) > poste_kernel_counts(biggest_cluster) ) &
        biggest_cluster = kid
    enddo

    ! Adjust largest cluster to ensure that the total kernel count is correct
    poste_kernel_counts(biggest_cluster) &
      = poste_kernel_counts(biggest_cluster) &
        + ( numbers_en - sum(poste_kernel_counts) )

    ! Check if total poste kernel count is correct
    if ( nint(sum( poste_kernel_counts ) - numbers_en) .ne. 0 ) then
      ! If not correct, exit
      if ( my_proc_id == 0 ) then
        write(*,*) 'Computed posterior cluster sizes do not sum to numbers_en.'
        write(*,*) 'Sum of cluster sizes = ', nint(sum( poste_kernel_counts ))
        write(*,*) '          numbers_en = ', numbers_en
        write(*,*) 'Exiting.'
      endif
      call exit()
    endif

    ! Update the posterior kernel weights
    poste_weights = poste_kernel_counts / (numbers_en*1.)


  end subroutine compute_posterior_weights
  ! ==============================================================================




  ! ==============================================================================
  ! Subroutine to determine how members are migrated between clusters
  ! General idea is to first identify the members that are migrating out of their
  ! prior cluster. Then, we absorb these migrating members into the expanding
  ! posterior cluster.
  subroutine setup_member_migration_paths()

    implicit none
    integer :: kid, prior_size, poste_size, n_remove, n_absorb
    integer :: selected_mem, mem_count
    real :: abs_pert, min_pert
    integer, dimension( numbers_en ) :: tmp_id
    logical :: flag_purgatory, flag_wrong_size

    ! Load prior cluster ids
    prior_id = kernel_id(iob, :)
    tmp_id = kernel_id(iob,:)
    poste_id = kernel_id(iob,:)


    ! Identify members migrating from shrinking clusters
    ! -------------------------------------------------
    cluster_shrink_tagging_loop: do kid = 1, max_kernel_num

      ! Skipping over non-shrinking clusters
      if ( poste_weights(kid) >= prior_weights(kid) ) &
        cycle cluster_shrink_tagging_loop

      ! Determine number of members in the prior cluster
      prior_size = nint( prior_weights(kid) * numbers_en )

      ! Determine number of members in the posterior cluster
      poste_size = nint( poste_weights(kid) * numbers_en )

      ! Determine number of members to remove
      n_remove = prior_size - poste_size

      if (my_proc_id == 0 .and. print_detail > 10) &
        write(*,'(4(a,x,i5,x))') 'Cluster',cid,'prior size',prior_size,&
                                  'poste size', poste_size, &
                                  'n_remove', n_remove

      ! Selecting members that have the smallest observation perturbations
      ! These members will be removed from the cluster after the Kalman gain
      ! computation step
      identify_outgoing_mems: do ii = 1, n_remove

        ! Iteratively select member with smallest sim obs perturbation
        min_pert = 99999; selected_mem = -9999
        select_smallest_pert: do n = 1, numbers_en

          ! Skipping over members belonging to other clusters or have been
          ! tagged for removal
          if ( tmp_id(n) .ne. kid ) cycle select_smallest_pert

          ! Check perturbation length
          abs_pert = abs( ya(iob,n) - prior_cluster_ymean(kid) )

          ! If abs pert length < min_pert, update selected mem
          if ( abs_pert < min_pert ) then
            selected_mem = n
            min_pert = abs_pert
          endif

        enddo select_smallest_pert

        ! Apply tag to the identified outgoing member
        tmp_id( selected_mem ) = -999
        poste_id( selected_mem ) = -999

      enddo identify_outgoing_mems

    enddo cluster_shrink_tagging_loop



    ! Now absorb the outgoing members into expanding clusters
    ! -------------------------------------------------------
    cluster_expand_tagging_loop: do kid = 1, max_kernel_num

      ! Ignore clusters that are not expanding
      if ( poste_weights(kid) <= prior_weights(kid) ) &
        cycle cluster_expand_tagging_loop

      ! Determine number of members to absorb
      prior_size = nint( prior_weights(kid) * numbers_en )
      poste_size = nint( poste_weights(kid) * numbers_en )
      n_absorb = poste_size - prior_size

      ! Select members to absorb into expanding cluster
      ii = 0  ! Variable for counting
      mem_absorb_loop: do n =1, numbers_en

        ! Skip over members without -999 tag
        if ( poste_id(n) .ne. -999 ) cycle mem_absorb_loop

        ! Absorb member
        poste_id( n ) = kid

        ! Increment counter
        ii = ii + 1

        ! If we have absorbed the requisite number of members, terminate loop
        if ( ii == n_absorb ) exit mem_absorb_loop

      enddo mem_absorb_loop

    enddo cluster_expand_tagging_loop



    ! Checking for unclustered members
    ! --------------------------------
    flag_purgatory = .false.
    check_purgatory: do n = 1, numbers_en
      if ( poste_id(n) == -999 ) then
        flag_purgatory = .true.
        if (my_proc_id == 0) then
          write(*,'(a,i4,a,i9)') &
            'Member ',n,' was not absorbed into any cluster for obs', iob
        endif
      endif
    enddo check_purgatory

    ! If there are unclustered members, exit code
    if (flag_purgatory) then
      if (my_proc_id==0) &
        write(*,*) 'Unclustered posterior members detected! Exiting.'
      !call MPI_BARRIER( comm, ierr )
      call exit()
    endif



    ! Checking if poste_id is consistent with poste_size
    ! --------------------------------------------------
    flag_wrong_size = .false.
    check_poste_size: do kid = 1, max_kernel_num

      ! Targetted posterior size
      poste_size = nint(poste_weights(kid) * numbers_en)

      ! Initialize counting variable
      mem_count = 0

      ! Iterate thru ensemble and identify posterior members in cluster
      do n = 1, numbers_en
        if ( poste_id(n) == kid ) mem_count = mem_count + 1
      enddo

      ! Check if correct number of poste_id's relating to cluster
      if ( poste_size .ne. mem_count ) then
        flag_wrong_size = .true.
        if ( my_proc_id==0 ) then
          write(*,'(a,i9,a,i3,a,i4,a,i4,a)') &
          'For obs ', iob,' cluster ',kid,', targetted poste size is', &
          poste_size,'. However, ',mem_count,' members counted!'
        endif
      endif

    enddo check_poste_size

    ! Exit if cluster has wrong number of members
    if (flag_wrong_size) then
      if (my_proc_id==0) &
        write(*,*) 'Incorrect posterior cluster sizes. Exiting.'
      !call MPI_BARRIER( comm, ierr )
      call exit()
    endif



    !! Remove overly small posterior clusters
    !! --------------------------------------
    !allocate( msg_id( numbers_en ) )
    !! Transmit the posterior id into the messenger variable
    !msg_id(:) = poste_id(:)
    !! Remove overly small posterior clusters
    !call remove_overly_small_clusters()
    !! Save the processed ids into poste_id
    !poste_id(:) = msg_id(:)
    !deallocate( msg_id )

    ! Update posterior cluster weights to reflect the ad-hoc adjusted sizes
    ! ---------------------------------------------------------------------
    do i = 1, max_kernel_num
      poste_weights(i) = count( poste_id == i ) * 1.0 / numbers_en
    enddo



    if (my_proc_id == 0) then
      write(*,'(a,x,2f)') 'prior weights', prior_weights
      write(*,'(a,x,2f)')'poste weights', poste_weights
    endif



  end subroutine setup_member_migration_paths
  ! ==============================================================================




  ! ==============================================================================
  ! Subroutine to identify the observation update zone within the slab
  subroutine identify_obs_update_zone()

    implicit none

    ! a. Identify ROI (special handling for Radiance obs)
    ! ---------------------------------------------------
    ngx = obs%roi(iob,1)
    ngz = obs%roi(iob,2)

    ! Special ROI handling for Radiance obs (ROI varies with obs)
    if (obstype=='Radiance  ') then
      ! Use Successive Covariance Localization (SCL) by Zhang et al. (2009)
      ! MWR for BT assimilation
      if(varname=='QCLOUD    ' .or. varname=='QRAIN     ' .or. varname=='QICE      ' .or. & !varname=='QVAPOR    ' .or.
                varname=='QGRAUP    ' .or. varname=='QSNOW     ') then
        if(obs%roi(iob,1) == 0) then
          update_flag = 0
        else
          ngx = obs%roi(iob,1)
        endif
      else
        if(obs%roi(iob,3) == 0) then
          update_flag = 0
        else
          ngx = obs%roi(iob,3)
        endif
      endif
      ! --- SCL end
      !
      !! --- Use Empirical Localization Functions (ELFs)
      !if (use_elf) then
      !  call corr_elf(varname, obs%sat(iob), obs%ch(iob), ya_ca(iob), ngx, ngz, obs%position(iob,3))
      !  if ( my_proc_id==0 .and. (varname=='U         '.or.varname=='V         '.or.varname=='QVAPOR    ')) &
      !      write(*,*) varname, ya_ca(iob), ngx, ngz, obs%position(iob,3)
      !endif
      !if (ngx == 0 .or. ngz == 0) update_flag = 0
      !! --- ELFs end
    endif


    ! b. Check if obs update zone overlaps with slab
    ! -----------------------------------------------
    ! Start and end indices of obs update zone
    ist = max( update_is, int(obs%position(iob,1))-ngx )
    ied = min( update_ie, int(obs%position(iob,1))+ngx )
    jst = max( update_js, int(obs%position(iob,2))-ngx )
    jed = min( update_je, int(obs%position(iob,2))+ngx )

    ! Check if obs update zone overlaps with slab using masks
    flag_slab_update = 0
    obs_mask = 0
    obs_mask(ist:ied,jst:jed) = 1
    flag_slab_update = sum( obs_mask(ist:ied,jst:jed) *  slab_mask(ist:ied,jst:jed) )


    ! c. Identify update zone in slab
    ! -------------------------------
    ! Note that these indices are in terms of slab's indices, not domain indices
    uist = ist - istart + 1
    uied = ied - istart + 1
    ujst = jst - jstart + 1
    ujed = jed - jstart + 1
    if ( uist < 1                  ) uist = 1
    if ( ujst < 1                  ) ujst = 1
    if ( uist >= iend - istart + 1 ) uist = iend - istart + 1
    if ( ujst >= jend - jstart + 1 ) ujst = jend - jstart + 1

    if ( uied < 1                  ) uied = 1
    if ( ujed < 1                  ) ujed = 1
    if ( uied >= iend - istart + 1 ) uied = iend - istart + 1
    if ( ujed >= jend - jstart + 1 ) ujed = jend - jstart + 1


    ! Check if variable is 2d or 3d
    !call wrf_var_dimension ( wrf_file, varname, ix, jx, kx, kxs, ii, jj, kk )
    ii = var_dims(m,1); jj = var_dims(m,2); kk = var_dims(m,3)
    if ( kk == 1 ) then
      kst  = 1
      ked = 1
    else
      kst = max( update_ks, int(obs%position(iob,3))-ngz )
      ked = min( update_ke, int(obs%position(iob,3))+ngz )
    endif

    ! Set up indices for sanity checks (needed later)
    check_i = nint( (uist + uied)/2. )
    check_j = nint( (ujst + ujed)/2. )
    check_k = nint( ( kst +  ked)/2. )

  end subroutine identify_obs_update_zone
  ! ==============================================================================



  ! ==============================================================================
  ! Subroutine to compute cluster cid's mean state slab, Kalman gain and
  ! square-root factor
  subroutine compute_kalman_update_parameters()

    implicit none
    integer :: cluster_size

    ! a. Compute prior kernel slab mean
    ! ---------------------------------
    prior_id = kernel_id(iob,:)
    prior_x_mean = 0.; cluster_size = 0
    estimate_prior_slab_mean: do n = 1, numbers_en

      ! Skipping over irrelevant members
      if ( prior_id(n) .ne. cid ) cycle estimate_prior_slab_mean

      ! Aggregate sample info
      prior_x_mean(uist:uied, ujst:ujed, kst:ked) &
        = prior_x_mean(uist:uied, ujst:ujed, kst:ked) &
          + x(uist:uied, ujst:ujed, kst:ked, m, n)
      cluster_size = cluster_size + 1

    enddo estimate_prior_slab_mean

    ! Estimate average
    prior_x_mean(uist:uied, ujst:ujed, kst:ked) &
      = prior_x_mean(uist:uied, ujst:ujed, kst:ked) / cluster_size

    ! b. Compute Kalman gain matrix
    ! -----------------------------
    ! Compute covariance for cluster cid
    km = 0.
    estimate_cov_xy: do n = 1, numbers_en

      ! Skipping over irrelevant members
      if ( prior_id(n) .ne. cid ) cycle estimate_cov_xy

      ! Aggregate sample info
      km( uist:uied, ujst:ujed, kst:ked ) &
        = km( uist:uied, ujst:ujed, kst:ked ) &
          + (  ( x(uist:uied, ujst:ujed, kst:ked, m, n)  &
                 - prior_x_mean(uist:uied, ujst:ujed, kst:ked) ) &
             * ( ya(iob,n) - prior_cluster_ymean(cid) ) &
            )
    enddo estimate_cov_xy

    ! Normalize the covariance
    km(uist:uied, ujst:ujed, kst:ked) &
      = km(uist:uied, ujst:ujed, kst:ked) / (cluster_size-1)


    ! c. Compute Kalman gain cluster cid
    ! -----------------------------------
    km(uist:uied, ujst:ujed, kst:ked) &
      = km(uist:uied, ujst:ujed, kst:ked) &
        / ( prior_cluster_yvariance(cid) + error*error )
    !if (my_proc_id == 0) then
    !  write(*,*) 'obs ',iob, ' yf variance ', prior_cluster_yvariance(cid)
    !  write(*,*) 'yf mean',  prior_cluster_ymean(cid), error*error
    !endif


    ! d. Compute square-root factor for cluster cid
    ! ---------------------------------------------
    alpha = error*error
    alpha = alpha / ( prior_cluster_yvariance(cid) + error*error )
    alpha = sqrt(alpha)
    alpha = (1.+ alpha)**-1

  end subroutine compute_kalman_update_parameters
  ! =============================================================================



  ! =============================================================================
  ! Subroutine to apply Kalman Filter update on cluster
  subroutine apply_enkf_update()

    implicit none
    real :: innovation

    ! Compute innovation
    ! ------------------
    innovation = obs%dat(iob) - prior_cluster_ymean(cid)

    ! Compute cluster cid posterior x mean
    ! ------------------------------------------
    poste_x_mean = 0.
    poste_x_mean( uist:uied, ujst:ujed, kst:ked ) &
      = prior_x_mean( uist:uied, ujst:ujed, kst:ked ) &
         + km( uist:uied, ujst:ujed, kst:ked ) * innovation

    ! Apply Kalman filter update to all members in cluster cid
    ! -----------------------------------------------------------
    kalman_filter_x_loop: do n = 1, numbers_en

      ! Skip over irrelevant memebrs
      if ( prior_id(n) .ne. cid ) cycle kalman_filter_x_loop

      ! Remove prior mean from member
      x(uist:uied, ujst:ujed, kst:ked, m, n) &
        = x(uist:uied, ujst:ujed, kst:ked, m, n) &
          - prior_x_mean( uist:uied, ujst:ujed, kst:ked )

      ! Apply perturbation update
      x(uist:uied, ujst:ujed, kst:ked, m, n) &
        = x(uist:uied, ujst:ujed, kst:ked, m, n) &
          - km( uist:uied, ujst:ujed, kst:ked) &
             * ( ya(iob,n) - prior_cluster_ymean(cid) ) &
             * alpha

      ! Add in posterior mean
      x(uist:uied, ujst:ujed, kst:ked, m, n) &
        = x(uist:uied, ujst:ujed, kst:ked, m, n) &
          + poste_x_mean( uist:uied, ujst:ujed, kst:ked )

    enddo kalman_filter_x_loop

  end subroutine apply_enkf_update
  ! =============================================================================




  ! =============================================================================
  ! Subroutine to recenter remaining members in shrinking cluster
  subroutine recenter_residual_members()

    implicit none
    integer :: prior_cluster_size, poste_cluster_size
    real, dimension( ni, nj, nk ) :: residual_mean

    ! a. Compute posterior mean
    ! --------------------------
    poste_x_mean = 0.
    prior_cluster_size = nint(prior_weights(cid) * numbers_en)
    estimate_poste_xmean: do n = 1, numbers_en

      ! SKip over irrelevant members
      if ( prior_id(n) .ne. cid ) cycle estimate_poste_xmean

      ! Aggregate stats
      poste_x_mean(uist:uied, ujst:ujed, kst:ked) &
        = poste_x_mean(uist:uied, ujst:ujed, kst:ked) &
          + x(uist:uied, ujst:ujed, kst:ked, m, n)

    enddo estimate_poste_xmean

    poste_x_mean(uist:uied, ujst:ujed, kst:ked) &
      = poste_x_mean(uist:uied, ujst:ujed, kst:ked) / prior_cluster_size


    ! b. Compute residual member poste mean
    ! -------------------------------------
    residual_mean = 0.
    poste_cluster_size = nint(poste_weights(cid) * numbers_en)
    estimate_residual_xmean: do n = 1, numbers_en

      ! SKip over irrelevant members
      if ( poste_id(n) .ne. cid ) cycle estimate_residual_xmean

      ! Aggregate stats
      residual_mean(uist:uied, ujst:ujed, kst:ked) &
        = residual_mean(uist:uied, ujst:ujed, kst:ked) &
          + x(uist:uied, ujst:ujed, kst:ked,m,n)

    enddo estimate_residual_xmean


    residual_mean(uist:uied, ujst:ujed, kst:ked) &
      = residual_mean(uist:uied, ujst:ujed, kst:ked) / poste_cluster_size


    ! c. Recenter residual members
    ! ----------------------------
    apply_residual_recentering: do n = 1, numbers_en

      ! SKip over irrelevant members
      if ( poste_id(n) .ne. cid ) cycle apply_residual_recentering

      ! Recenter mean
      x(uist:uied, ujst:ujed, kst:ked,m,n) &
        = x(uist:uied, ujst:ujed, kst:ked,m,n) &
          - residual_mean(uist:uied, ujst:ujed, kst:ked) &
          + poste_x_mean(uist:uied, ujst:ujed, kst:ked)

    enddo apply_residual_recentering


    ! d. Sanity check the residual members' mean
    ! ------------------------------------------
    finish_xavg = 0.; origin_xavg = poste_x_mean( check_i, check_j, check_k )

    ! Compute post-transfer shrinking cluster mean
    sanity_check_residual_mem_mean: do n=1, numbers_en

      ! SKip over irrelevant members
      if ( poste_id(n) .ne. cid ) cycle sanity_check_residual_mem_mean

      finish_xavg = finish_xavg + x( check_i, check_j, check_k, m, n )

    enddo sanity_check_residual_mem_mean
    finish_xavg = finish_xavg / poste_cluster_size

    ! Compare the post-transfer mean against pre-transfer mean
    if ( abs( finish_xavg - origin_xavg )/abs(origin_xavg) > 1e-3 &
         .and. abs(origin_xavg) > 1e-3  ) then
       write(*,*) 'ERROR: Shrinking cluster mean shifted!'
       write(*,*) 'Variable: ', varname
       write(*,"(2(a,x,f,x))") &
         'Pre-transfer mean:', origin_xavg, &
         'Post-transfer mean:', finish_xavg
       write(*,*) 'Exiting.'
       call exit(ierr)
     endif

  end subroutine recenter_residual_members
  ! =============================================================================



  ! =============================================================================
  ! Subroutine to generate parameters controlling resampling of expanding cluster
  subroutine expansion_resampling_params()

    implicit none
    integer :: nrows, ncols


    ! Figure out how many members to shift into cluster
    ! -------------------------------------------------
    n_shift = nint( poste_weights(cid) * numbers_en ) &
              - nint( prior_weights(cid) * numbers_en )


    ! Figure out how many members to use for resampling n_shift mems
    ! --------------------------------------------------------------
    ! See Eqn (18) of Chan et al (2020)
    prior_cluster_size = nint( prior_weights(cid) * numbers_en )
    if (n_shift .le. prior_cluster_size) then
      n_src_mems = n_shift - 1
    else
      n_src_mems = prior_cluster_size
    endif

    ! Compute expansion inflation factor
    ! ----------------------------------
    ! See definition after Eqn (13)
    exp_infl_fac = n_shift + prior_cluster_size - 1.
    exp_infl_fac = exp_infl_fac / (prior_cluster_size - 1.)
    exp_infl_fac = sqrt( exp_infl_fac )

    nrows = n_src_mems; ncols = n_shift


  end subroutine expansion_resampling_params
  ! =============================================================================





  ! =============================================================================
  ! Subroutine to generate matrix E referenced in Chan et al (2020)
  subroutine gen_resampling_matrixE( E )

    implicit none

    ! Output matrix E
    real, dimension( n_src_mems, n_shift ), intent(inout) :: E

    ! Convenient variables
    integer :: lapack_info
    integer :: nrows, ncols
    ! Note nrows=m, ncols = n_shift, prior_cluster_size = N^F_{E,e}
    ! in the Chan et al (2020) BGEnKF paper

    ! Important matrices
    ! orthoW = (L_W)^-1 * W in Chan et al.
    ! E and E' follow Chan et al.
    real, dimension(n_src_mems, n_shift) :: orthoW, primeE

    ! Variables for sanity checking
    real, dimension( n_src_mems, n_src_mems ) :: EET


    ! Translate variables
    ! -------------------
    nrows = n_src_mems; ncols = n_shift



    ! Generate orthoW = (L_W)^-1 * W.
    ! -------------------------------
    call gen_orthoW_matrix( nrows, ncols, orthoW)


    ! Generate E' = L_E * (L_W)^-1 * W = L_E * orthoW
    ! -----------------------------------------------
    call gen_primeE_matrix( nrows, ncols, orthoW, primeE)


    ! E = barE + primeE
    ! -----------------
    E = primeE + (exp_infl_fac -1) / ncols


    ! Sanity check: is Eqn (B1) satisfied?
    ! ------------------------------------
    ! Eqn (B1): E * E^T = nshift/(N^f_E,e -1) I
    call sgemm( 'n','t', nrows, nrows, ncols, 1., E, nrows, &
                E, nrows, 0., EET, nrows )
    do ii = 1, nrows
      EET(ii,ii) = EET(ii,ii) - (ncols*1.)/(prior_cluster_size-1.)
    enddo

    if ( sum(abs(EET)) > 1e-3) then
      write(*,*) 'ERROR: E * E^T != nshift/(N^f_E,e -1) I'
      write(*,*) 'Exiting'
      call exit()
    endif


    ! Sanity check: is Eqn (B2) satisfied?
    ! ------------------------------------
    ! Eqn (B2): Sum of columns of E = k-1
    EET = 0.
    do ii = 1, ncols
      EET(:,1) = EET(:,1) + E(:,ii)
    enddo

    EET(:,1) = EET(:,1) - (exp_infl_fac-1.)

    if ( sum(abs(EET(:,1))) > 1e-3 ) then
      write(*,*) 'ERROR: sum of cols in E != k-1'
      write(*,*) 'Exiting'
      call exit()
    endif

  end subroutine gen_resampling_matrixE
  ! =============================================================================




  ! ============================================================================
  ! Subroutine to generate W s.t. WWT is identity.
  ! I.e., WT is the right inverse of W
  subroutine gen_orthoW_matrix( nrows, ncols, orthoW )

    implicit none

    ! Input arguments
    integer, intent(in) :: nrows, ncols
    real, dimension(nrows, ncols), intent(inout):: orthoW

    ! Matrices needed to generate orthogonal W
    real, dimension(nrows, ncols) :: W
    real, dimension(nrows, nrows) :: WWT, L_W, WWT2, inv_L_W

    ! Variables for counting
    integer :: lapack_info

    ! Variable for sanity checkin
    real :: check



    ! Construct W and make sure that columns sum to 0
    ! -----------------------------------------------
    ! Using a pseudo truncated identity matrix
    W = 0.
    do ii = 1, nrows
        W(ii,ii) = 1.
    enddo
    ! Ensuring sum of all columns is zero
    W = W - 1./ncols



    ! Cholesky decomp WWT
    ! -------------------
    WWT=0.
    call sgemm( 'n', 't', nrows, nrows, ncols, 1., W, nrows, &
                 W, nrows, 0., WWT, nrows )

    L_W = WWT
    call spotrf( 'L', nrows, L_W, nrows, lapack_info )

    ! Zeroing unneeded part of matrix
    do ii =1,nrows
      do jj = ii+1, nrows
        L_W(ii,jj) = 0.
      enddo
    enddo



    ! Adjust W s.t. WWT = I
    ! ---------------------
    ! Lower triangle invert L_W
    inv_L_W = L_W
    call strtri( 'L', 'N', nrows, inv_L_W, nrows, lapack_info )

    ! Apply inv_L_W onto W matrix to generate orthoW
    call sgemm( 'n','n',nrows, ncols, nrows, 1., inv_L_W, nrows, &
                 W, nrows, 0., orthoW, nrows)



    ! Sanity check1: is orthoW * orthoW^T = I?
    ! ----------------------------------------
    WWT2 = 0.
    call sgemm( 'n', 't', nrows, nrows, ncols, 1., orthoW, nrows, &
                 orthoW, nrows, 0., WWT2, nrows )

    ! Check if WWT2 is identity
    do ii = 1, nrows
      WWT2(ii,ii) = WWT2(ii,ii) -1.
    enddo
    check = sum(abs(WWT2))
    if ( check > 1e-3 ) then
      write(*,*) 'ERROR: orthoW *orthoW^T is not identity!'
      write(*,*) 'Exiting'
      call exit()
    endif

    ! Sanity check 2: Is the sum of columns in orthoW equals zero?
    ! ------------------------------------------------------------
    W = orthoW
    do jj =2, ncols
      W(:,1) = W(:,1) + orthoW(:,jj)
    enddo
    check = sum( abs( W(:,1) ) )
    if ( check > 1e-3 ) then
      write(*,*) 'ERROR: cols of orthoW do not sum to zero!'
      write(*,*) 'Exiting'
      call exit()
    endif

  end subroutine gen_orthoW_matrix




  ! ============================================================================
  ! Subroutine to constuct E' discussed in Chan et al (2020) appendix
  ! B. Will use orthoW generated from gen_orthoW_matrix
  subroutine gen_primeE_matrix( nrows, ncols, orthoW, primeE )

    implicit none

    ! Input arguments
    integer, intent(in) :: nrows, ncols
    real, dimension(nrows, ncols), intent(in)    :: orthoW
    real, dimension(nrows, ncols), intent(inout) :: primeE

    ! Variables relating to C_E and L_E
    real, dimension(nrows, nrows) :: C_E, L_E, check

    ! Misc variables
    integer :: lapack_info


    ! Generate the C_E matrix from Eqn (B4).
    ! -------------------------------------
    C_E = 0
    do ii = 1, nrows
      C_E(ii,ii) = ncols*1. / (prior_cluster_size - 1.0)
    enddo
    C_E = C_E - ( (exp_infl_fac - 1.)**2 ) / (ncols*1.)


    ! Cholesky decomp C_E to generate L_E
    ! -----------------------------------
    L_E = C_E
    call spotrf( 'L', nrows, L_E, nrows, lapack_info )
    ! Zeroing unneeded things
    do ii =1,nrows
      do jj = ii+1, nrows
        L_E(ii,jj) = 0.
      enddo
    enddo


    ! Construct E' matrix using E' = L_E * orthoW
    ! -------------------------------------------
    call sgemm( 'n','n', nrows, ncols, nrows, 1., L_E, nrows, &
                orthoW, nrows, 0., primeE, nrows )


    ! Sanity check: is E' * E'^T = C_E?
    ! -------------------------------
    call sgemm( 'n','t', nrows, nrows, ncols, 1., primeE, nrows, &
                primeE, nrows, 0., check, nrows )

    if ( sum( abs( check - C_E) ) > 1e-3 ) then
      write(*,*) "ERROR: E' * E'^T != C_E"
      write(*,*) "Exiting"
      call exit()
    endif


  end subroutine gen_primeE_matrix
  ! ============================================================================



  ! ============================================================================
  ! Subroutine to estimate mean state of cluster
  ! Stores mean state in poste_x_mean
  subroutine estimate_cluster_xa_mean()

    implicit none

    ! Compute posterior ensemble mean
    ! -------------------------------
    poste_x_mean = 0.
    prior_cluster_size = 0
    estimate_xa_mean: do n = 1, numbers_en

      ! Ignore irrelevant mems
      if( prior_id(n) .ne. cid ) cycle estimate_xa_mean

      ! Aggregate stats
      poste_x_mean( uist:uied, ujst:ujed, kst:ked ) &
        = poste_x_mean( uist:uied, ujst:ujed, kst:ked ) &
          + x( uist:uied, ujst:ujed, kst:ked, m, n )

      prior_cluster_size = prior_cluster_size + 1

    enddo estimate_xa_mean
    poste_x_mean(uist:uied, ujst:ujed, kst:ked) &
      = poste_x_mean(uist:uied, ujst:ujed, kst:ked) / prior_cluster_size


  end subroutine estimate_cluster_xa_mean
  ! ===========================================================================



  ! ===========================================================================
  ! Subroutine to compute expanding cluster mean and variance before cluster
  ! expansion. Will be used for sanity checks later.
  subroutine sanity_check_cluster_expansion_preamble()

    implicit none

    ! Compute x mean before cluster transfer
    origin_xavg = sum( x( check_i, check_j, check_k, m, 1:numbers_en), &
                       mask = prior_id == cid )
    origin_xavg = origin_xavg / ( numbers_en * prior_weights(cid) )


    ! Compute x var before cluster transfer
    origin_xvar = sum( ( x( check_i, check_j, check_k, m, 1:numbers_en) &
                         - origin_xavg ) ** 2, &
                       mask = prior_id == cid )
    origin_xvar = origin_xvar / ( numbers_en * prior_weights(cid) - 1 )


  end subroutine sanity_check_cluster_expansion_preamble
  ! ===========================================================================



  ! ============================================================================
  ! Subroutine to convert cluster states to perturbation states
  ! Stores mean state in poste_x_mean
  subroutine cluster_states_to_perturbations()

    implicit none

    remove_mean_loop: do n = 1, numbers_en

      ! Ignore irrelevant mems
      if( prior_id(n) .ne. cid ) cycle remove_mean_loop

      x( uist:uied, ujst:ujed, kst:ked, m, n ) &
        = x( uist:uied, ujst:ujed, kst:ked, m, n ) &
          - poste_x_mean( uist:uied, ujst:ujed, kst:ked )

    enddo remove_mean_loop


  end subroutine cluster_states_to_perturbations
  ! ===========================================================================




  ! ============================================================================
  ! Subroutine to inflate (prior_cluster_size - n_src_mems) perturbations of the
  ! expanding cluster
  subroutine inflate_cluster_perts()

    implicit none

    integer :: infl_count, infl_max

    ! Number of perts to inflate
    infl_max = prior_cluster_size - n_src_mems
    if (print_detail > 10) &
      write(*,'(a,x,i5)') 'Number of members to inflate', infl_max

    ! Initialize counter
    infl_count = 0;

    ! Initialize indices where the perts were inflated
    inds_infl = -999



    ! Inflate first infl_max perturbations
    ! ---------------------------------------
    flag_infl = 0
    inflate_cluster_perts_loop: do n = 1, numbers_en

      ! Ignore irrelevant mems
      if( prior_id(n) .ne. cid ) cycle inflate_cluster_perts_loop

      ! Check if we have inflated all members
      if ( infl_count >= infl_max ) exit inflate_cluster_perts_loop

      ! Apply inflation factor
      x(uist:uied, ujst:ujed, kst:ked, m, n) &
        = x(uist:uied, ujst:ujed, kst:ked, m, n)  * exp_infl_fac

      ! Increment counter
      infl_count = infl_count + 1

      ! Mark this member as an inflated member
      inds_infl(infl_count) = n
      flag_infl(n) = 1

      ! Check if we have inflated all members
      if ( infl_count >= infl_max ) exit inflate_cluster_perts_loop


    enddo inflate_cluster_perts_loop

    if( infl_count .ne. infl_max) then
      write(*,*) 'infl_count',infl_count
      write(*,*) 'infl_max', infl_max
      call MPI_ABORT(comm, i, ierr )
    endif


  end subroutine inflate_cluster_perts
  ! ============================================================================




  ! =============================================================================
  ! Subroutine to construct new members from uninflated members
  subroutine construct_new_perturbations()

    implicit none

    integer :: new_ind, src_ind
    real, dimension( n_src_mems, n_shift ) :: E


    ! Exception case: if only 1 incoming member
    ! -----------------------------------------
    solo_incoming: if (n_shift == 1) then

      solo_incoming_loop: do n = 1, numbers_en
        if ( (prior_id(n) .ne. cid) .and. (poste_id(n) == cid) ) then
          x(uist:uied, ujst:ujed, kst:ked, m, n) = 0.
          exit solo_incoming_loop
        endif
      enddo solo_incoming_loop

      return ! Return to main program
    endif solo_incoming



    ! a. Generate flags to id uninflated original members
    ! ---------------------------------------------------
    inds_uninfl = -999
    ii = 0
    identify_uninflated_mems: do n = 1, numbers_en

      ! Skip over irrelevant mems
      if ( prior_id(n) .ne. cid ) cycle identify_uninflated_mems

      ! If member is inflated, skip over
      if ( flag_infl(n) == 1 ) cycle identify_uninflated_mems

      ! Save the index of uninflated member
      ii = ii +1
      inds_uninfl(ii) = n

    enddo identify_uninflated_mems




    ! b. Generate flags to mark incoming members and zero them
    ! --------------------------------------------------------
    inds_incoming = -999
    ii = 0
    flag_incoming_members: do n = 1, numbers_en
      if ( (prior_id(n) .ne. cid) .and. (poste_id(n) == cid) ) then
        ii=ii+1
        inds_incoming(ii) = n
        x(uist:uied, ujst:ujed, kst:ked, m, n) = 0.
      endif
    enddo flag_incoming_members
    if (print_detail > 10) then
      write(*,'(a,i5)') 'Number of incoming members:', &
                        count( inds_incoming > 0 )
      write(*,'(a,i5)') 'Number of uninflated members:', &
                        count( inds_uninfl > 0 )

    endif



    ! c. Construct new perturbation by resampling
    ! --------------------------------------------
    ! Set up resampling matrix
    call gen_resampling_matrixE( E )


    ! Apply resampling matrix
    make_new_perts: do ii = 1, n_shift

      ! Index of new member
      new_ind = inds_incoming(ii)

      ! Construct new perturbations
      ! Each perturbation is a linear combination of uninflated members
      ! times coefficients from a column of E
      linear_sum_perts: do jj = 1, n_src_mems

        ! Contribution from member
        src_ind = inds_uninfl(jj)

        ! Apply contribution from member
        x(uist:uied, ujst:ujed, kst:ked, m, new_ind) &
          = x(uist:uied, ujst:ujed, kst:ked, m, new_ind) &
            + E(jj, ii) * x(uist:uied, ujst:ujed, kst:ked, m, src_ind)

      enddo linear_sum_perts

    enddo make_new_perts

  end subroutine construct_new_perturbations
  ! =============================================================================




  ! =============================================================================
  ! Subroutine to convert cluster perturbations into states
  subroutine poste_cluster_perturbations_to_states()

    implicit none

    add_xa_mean: do n = 1, numbers_en

      ! ignore irrelevant members
      if ( poste_id(n) .ne. cid ) cycle add_xa_mean

      ! Add posterior mean
      x(uist:uied, ujst:ujed, kst:ked, m, n) &
        = x(uist:uied, ujst:ujed, kst:ked, m, n) &
          + poste_x_mean(uist:uied, ujst:ujed, kst:ked)

    enddo add_xa_mean

  end subroutine poste_cluster_perturbations_to_states
  ! =============================================================================




  ! ===========================================================================
  ! Subroutine to check if the cluster expansion procedure preserved the
  ! expanding cluster's mean and covariance.
  subroutine sanity_check_cluster_expansion()

    implicit none

    logical :: exit_flag

    ! Compute x mean after cluster transfer
    finish_xavg = sum( x( check_i, check_j, check_k, m, 1:numbers_en), &
                       mask = poste_id == cid )
    finish_xavg = finish_xavg / ( numbers_en * poste_weights(cid) )


    ! Compute x var before cluster transfer
    finish_xvar = sum( ( x( check_i, check_j, check_k, m, 1:numbers_en) &
                         - finish_xavg ) ** 2, &
                       mask = poste_id == cid )
    finish_xvar = finish_xvar / ( numbers_en * poste_weights(cid) - 1 )


    ! Check if the cluster mean are changed by the cluster expansion
    exit_flag=.false.
    if ( abs( finish_xavg - origin_xavg )/abs(origin_xavg) > 1e-3 &
         .and. abs( origin_xavg ) > 1e-3 ) then
      write(*,'(a,x,a)') 'Cluster expansion changed mean value for', varname
      write(*,'(2(a,x,f,x))') 'Old mean', origin_xavg, &
                              'New mean', finish_xavg
      write(*,'(2(a,x,i4,4x))') 'Old size', int(prior_weights(cid) * numbers_en), &
                                'New size', int(poste_weights(cid) * numbers_en)
      exit_flag = .true.
    endif

    ! Check if the cluster variance are changed by the cluster expansion
    if ( abs( finish_xvar - origin_xvar )/abs(origin_xvar) > 1e-3 &
         .and. abs(origin_xvar) > 1e-3 ) then
      write(*,'(a,x,a)') 'Cluster expansion changed variance value for', varname
      write(*,'(3(a,x,f,x))') 'Old variance', origin_xvar, &
                              'New variance', finish_xvar, &
                              'Old-to-new variance ratio', finish_xvar / origin_xvar
      write(*,*) 'perturbations', &
                 pack( x( check_i, check_j, check_k, m, 1:numbers_en), &
                       poste_id == cid )

      write(*,'(2(a,x,i4,4x))') 'Old size', int(prior_weights(cid) * numbers_en), &
                                'New size', int(poste_weights(cid) * numbers_en)
      write(*,'(2(a,x,i4,4x))') 'Old size', count( prior_id==cid ), &
                                'New size', count( poste_id==cid )

      exit_flag = .true.
    endif

    ! Exiting if sanity check failed
    if ( exit_flag ) then
      write(*,*) 'Exiting program'
      call MPI_ABORT( comm,i, ierr )
      call exit(ierr)
    endif


  end subroutine sanity_check_cluster_expansion
  ! ===========================================================================









  ! ===============================================================================
  ! Subroutine to compute auxiliary variable used in clustering
  ! Currently using column-integrated ice.
  subroutine compute_auxiliary_variable( iid0, jid0 )

    implicit none

    ! Input argument to indicate  process' slab
    integer, intent(in) :: iid0, jid0
    ! Array to hold kernel ids within a slab. Will be sent to all procs
    integer, allocatable, dimension(:,:) :: kernel_id_send
    ! Useful arrays
    real, allocatable, dimension(:,:) :: pres, tk, tv, qhydro, qvapor, delz
    real, allocatable, dimension(:,:) :: air_density !, liq_ice_content
    ! Counting variables
    integer :: iob, kk, ie
    ! Position related variables
    integer :: istart, iend, jstart, jend, obs_ii, obs_jj, slab_ii, slab_jj
    ! Useful constants
    real, parameter :: P1000MB=100000.D0
    real, parameter :: R_D=287.D0
    real, parameter :: Cpd=7.D0*R_D/2.D0
    real, parameter :: Re=6378000.0
    ! Useful array for flagging
    real, allocatable, dimension(:,:) :: ones_array
    real, allocatable, dimension(:) :: flag_deep, flag_congestus, flag_cloud


    ! Store iid and jid values
    iid = iid0; jid = jid0

    ! Allocate auxiliary variables
    if( allocated( aux      ) )deallocate(aux)
    if( allocated( aux_send ) )deallocate(aux_send)
    allocate(      aux(obs%num,numbers_en+1) );      aux = 0.
    allocate( aux_send(obs%num,numbers_en+1) ); aux_send = 0.

    ! Allocate useful variable arrays
    allocate( pres            ( kx, numbers_en+1 ) ); pres            = 0
    allocate( tv              ( kx, numbers_en+1 ) ); tv              = 0
    allocate( tk              ( kx, numbers_en+1 ) ); tk              = 0
    allocate( qhydro          ( kx, numbers_en+1 ) ); qhydro          = 0
    allocate( qvapor          ( kx, numbers_en+1 ) ); qvapor          = 0
    !allocate( liq_ice_content ( kx, numbers_en+1 ) ); liq_ice_content = 0
    allocate( air_density     ( kx, numbers_en+1 ) ); air_density     = 0
    allocate( delz            ( kx, numbers_en+1 ) ); delz            = 0
    allocate( ones_array      ( kx, numbers_en   ) ); ones_array      = 0
    allocate( flag_deep       ( numbers_en       ) ); flag_deep       = 0
    allocate( flag_congestus  ( numbers_en       ) ); flag_congestus  = 0
    allocate( flag_cloud      ( numbers_en       ) ); flag_cloud      = 0
    ones_array = 1.

    ! Figure out the slab's start and end indices (in terms of full domain
    ! indices)
    istart=iid*ni+1
    iend=(iid+1)*ni
    jstart=jid*nj+1
    jend=(jid+1)*nj
    if(iid==(nicpu-1)) iend=ix+1
    if(jid==(njcpu-1)) jend=jx+1


    ! Compute auxiliary variable (col-integrated frozen water)
    ! --------------------------------------------------------
    loop_compute_aux_variable: do iob = 1, obs%num

      !! Skip over non-radiance obs
      !if ( trim(obs%type(iob)) .ne. 'Radiance' ) &
      !  cycle obs_col_clustering_loop

      ! Check if the observation lives within the slab
      ! If obs is outside slab, skip over
      obs_ii = int( obs%position(iob, 1) )
      obs_jj = int( obs%position(iob, 2) )
      if (  ( ( istart > obs_ii ) .or. ( iend < obs_ii ) ) &
          .or. ( ( jstart > obs_jj ) .or. ( jend < obs_jj ) ) ) then
        cycle loop_compute_aux_variable
      endif


      ! Compute slab-relative coordinates
      slab_ii = obs_ii - istart + 1
      slab_jj = obs_jj - jstart + 1


      ! Setup useful variables
      ! ----------------------
      ! Load pressure
      pres = x( slab_ii, slab_jj, 1:kx, ind_p, : )   &
             + x( slab_ii, slab_jj, 1:kx, ind_pb, : )

      ! Load temperature
      tk = ( x( slab_ii, slab_jj, 1:kx, ind_t, : ) + 300.0 ) &
           * ( (pres / P1000MB) ** (R_D/Cpd) )


      ! Load QVAPOR
      qvapor = x( slab_ii, slab_jj, 1:kx, ind_qvapor, : )


      ! Load all ice mixing ratios
      qhydro = x( slab_ii, slab_jj, 1:kx, ind_qice  , : ) &
               + x( slab_ii, slab_jj, 1:kx, ind_qsnow , : ) &
               + x( slab_ii, slab_jj, 1:kx, ind_qgraup, : )


      ! Compute model layer thicknesses
      delz(2:kx,:) = (  ( x(slab_ii, slab_jj, 3:(kx+1), ind_ph, :) &
                          + x(slab_ii, slab_jj, 3:(kx+1), ind_phb, :) ) &
                      - ( x(slab_ii, slab_jj, 2:kx, ind_ph, :) &
                          + x(slab_ii, slab_jj, 2:kx, ind_phb, :) ) &
                     ) / 9.806
      delz(1,:) = ( x(slab_ii, slab_jj, 1, ind_ph, :) &
                    + x(slab_ii, slab_jj, 1, ind_phb, :) )/9.806 &
                  - x( slab_ii, slab_jj, 1, ind_hgt, :)


      ! Compute virtual temperature with no approximations (see wikipedia)
      tv = tk * (qvapor + 0.622) / ( 0.622 * ( 1.+qvapor) )

      ! Compute air density using modified ideal gas law
      air_density = pres / (287*tv)


      ! Compute auxiliary variable (column frozen mass content)
      aux_send(iob, :) = sum( qhydro(:,:) * air_density(:,:)  &
                              * delz(:,:)  , DIM=1 )


    enddo loop_compute_aux_variable



    ! Now pass auxiliary info to all processes.
    call MPI_Allreduce( aux_send, aux, obs%num*numbers_en,  &
                        MPI_REAL, MPI_SUM, comm, ierr )


    ! Release memory
    deallocate( pres            )
    deallocate( tv              )
    deallocate( tk              )
    deallocate( qhydro          )
    deallocate( qvapor          )
    !deallocate( liq_ice_content )
    deallocate( air_density     )
    deallocate( delz            )
    deallocate( ones_array      )
    deallocate( flag_deep       )
    deallocate( flag_congestus  )
    deallocate( flag_cloud      )

    call MPI_BARRIER( comm, ierr)
  end subroutine compute_auxiliary_variable
  ! ===============================================================================







  ! ==============================================================================
  ! SEQUENTIAL FILTER SPECIFIC FUNCTION
  ! -----------------------------------
  ! Subroutine to separate members into clear sky, congestus and deep cloud
  ! clusters for each observation in structured variable obs.
  ! Results are stored in allocatable array kernel_id
  ! 1 -- clear, 2 -- congestus, 3 -- deep cloud
  subroutine tag_clear_congestus_deep_clusters( iid0, jid0 )

    implicit none

    ! Input argument to indicate  process' slab
    integer, intent(in) :: iid0, jid0
    ! Array to hold kernel ids within a slab. Will be sent to all procs
    integer, allocatable, dimension(:,:) :: kernel_id_send
    ! Useful arrays
    real, allocatable, dimension(:,:) :: pres, tk, tv, qhydro, qvapor, delz
    real, allocatable, dimension(:,:) :: air_density, liq_ice_content
    ! Counting variables
    integer :: iob, kk, ie
    ! Position related variables
    integer :: istart, iend, jstart, jend, obs_ii, obs_jj, slab_ii, slab_jj
    ! Useful constants
    real, parameter :: P1000MB=100000.D0
    real, parameter :: R_D=287.D0
    real, parameter :: Cpd=7.D0*R_D/2.D0
    real, parameter :: Re=6378000.0
    ! Useful array for flagging
    real, allocatable, dimension(:,:) :: ones_array
    real, allocatable, dimension(:) :: flag_deep, flag_congestus, flag_cloud


    ! Store iid and jid values
    iid = iid0; jid = jid0

    ! Allocate kernel_id array and array for MPI transmission
    if( allocated(kernel_id) ) deallocate(kernel_id)
    allocate( kernel_id      ( obs%num, numbers_en ) )
    allocate( kernel_id_send ( obs%num, numbers_en ) )
    kernel_id = 0; kernel_id_send = 0

    ! Allocate useful variable arrays
    allocate( pres            ( kx, numbers_en+1 ) ); pres            = 0
    allocate( tv              ( kx, numbers_en+1 ) ); tv              = 0
    allocate( tk              ( kx, numbers_en+1 ) ); tk              = 0
    allocate( qhydro          ( kx, numbers_en+1 ) ); qhydro          = 0
    allocate( qvapor          ( kx, numbers_en+1 ) ); qvapor          = 0
    allocate( liq_ice_content ( kx, numbers_en+1 ) ); liq_ice_content = 0
    allocate( air_density     ( kx, numbers_en+1 ) ); air_density     = 0
    allocate( delz            ( kx, numbers_en+1 ) ); delz            = 0
    allocate( ones_array      ( kx, numbers_en   ) ); ones_array      = 0
    allocate( flag_deep       ( numbers_en       ) ); flag_deep       = 0
    allocate( flag_congestus  ( numbers_en       ) ); flag_congestus  = 0
    allocate( flag_cloud      ( numbers_en       ) ); flag_cloud      = 0
    ones_array = 1.

    ! Figure out the slab's start and end indices (in terms of full domain
    ! indices)
    istart=iid*ni+1
    iend=(iid+1)*ni
    jstart=jid*nj+1
    jend=(jid+1)*nj
    if(iid==(nicpu-1)) iend=ix+1
    if(jid==(njcpu-1)) jend=jx+1


    ! Assign cluster id to members for each observation
    ! -------------------------------------------------
    obs_col_clustering_loop: do iob = 1, obs%num

      !! Skip over non-radiance obs
      !if ( trim(obs%type(iob)) .ne. 'Radiance' ) &
      !  cycle obs_col_clustering_loop

      ! Check if the observation lives within the slab
      ! If obs is outside slab, skip over
      obs_ii = int( obs%position(iob, 1) )
      obs_jj = int( obs%position(iob, 2) )
      if (  ( ( istart > obs_ii ) .or. ( iend < obs_ii ) ) &
          .or. ( ( jstart > obs_jj ) .or. ( jend < obs_jj ) ) ) then
        cycle obs_col_clustering_loop
      endif


      ! Compute slab-relative coordinates
      slab_ii = obs_ii - istart + 1
      slab_jj = obs_jj - jstart + 1


      ! Setup useful variables
      ! ----------------------
      ! Load pressure
      pres = x( slab_ii, slab_jj, 1:kx, ind_p, : )   &
             + x( slab_ii, slab_jj, 1:kx, ind_pb, : )

      ! Load temperature
      tk = ( x( slab_ii, slab_jj, 1:kx, ind_t, : ) + 300.0 ) &
           * ( (pres / P1000MB) ** (R_D/Cpd) )

      ! Load QVAPOR
      qvapor = x( slab_ii, slab_jj, 1:kx, ind_qvapor, : )

      ! Load all liq+ice mixing ratios
      qhydro = x( slab_ii, slab_jj, 1:kx, ind_qcloud  , : ) &
               + x( slab_ii, slab_jj, 1:kx, ind_qrain , : ) &
               + x( slab_ii, slab_jj, 1:kx, ind_qice  , : ) &
               + x( slab_ii, slab_jj, 1:kx, ind_qsnow , : ) &
               + x( slab_ii, slab_jj, 1:kx, ind_qgraup, : )

      ! Compute model layer thicknesses
      delz(2:kx,:) = (  ( x(slab_ii, slab_jj, 3:(kx+1), ind_ph, :) &
                          + x(slab_ii, slab_jj, 3:(kx+1), ind_phb, :) ) &
                      - ( x(slab_ii, slab_jj, 2:kx, ind_ph, :) &
                          + x(slab_ii, slab_jj, 2:kx, ind_phb, :) ) &
                     ) / 9.806
      delz(1,:) = ( x(slab_ii, slab_jj, 2, ind_ph, :) &
                    + x(slab_ii, slab_jj, 2, ind_phb, :) )/9.806 &
                  - x( slab_ii, slab_jj, 1, ind_hgt, :)

      ! Compute virtual temperature with no approximations (see wikipedia)
      tv = tk * (qvapor + 0.622) / ( 0.622 * ( 1.+qvapor) )

      ! Compute air density using modified ideal gas law
      air_density = pres / (287*tv)


      ! Compute liq+ice water content
      liq_ice_content = qhydro * air_density * delz


      ! Separate members into clusters for this obs
      ! -------------------------------------------
      ! Iterate over members
      mem_cluster_loop: do ie = 1, numbers_en
  
        ! Seek deep and dense clouds
        seek_deep_dense_clouds: do kk = 1, kx
  
          ! If deep and dense clouds detected, note and exit loop
          ! DEPRECATED!
          if ( (pres(kk,ie) < 700*100.) &
                 .and. ( liq_ice_content(kk,ie) > 0.01 ) ) then
            kernel_id_send( iob, ie ) = 2
            exit seek_deep_dense_clouds
          endif
  
        enddo seek_deep_dense_clouds

        ! If member has been assigned deep cloud, move on to next member
        if ( kernel_id_send( iob, ie ) > 0 ) cycle mem_cluster_loop

  
        ! Seek congestus clouds
        seek_congestus_clouds: do kk = 1, kx

          ! If congested clouds detected, note and exit loop
          ! Written for window-BT
          if ( ( pres(kk,ie) < 800*100. ) &
                 .and. ( liq_ice_content(kk,ie) > 0.001 ) ) then
            kernel_id_send( iob, ie ) = 1
            exit seek_congestus_clouds
          endif

        enddo seek_congestus_clouds

        ! If member has neither deep nor congestus clouds, the kernel_id_send
        ! is 0

      enddo mem_cluster_loop

    enddo obs_col_clustering_loop


    ! Now pass the kernel id information to all processes
    call MPI_Allreduce( kernel_id_send, kernel_id, obs%num*numbers_en,  &
                        MPI_INTEGER, MPI_SUM, comm, ierr )

    ! Shifting values
    kernel_id = kernel_id + 1


    ! Ensuring that all clusters have at least min_mem_per_kernel members or
    ! have no members at all
    allocate( msg_id( numbers_en ) )
    do iob = 1, obs%num
      ! Save kernel id into the messenger variable msg_id
      msg_id(:) = kernel_id(iob,:)
      ! Remove overly small clusters
      call remove_overly_small_clusters()
      ! Save messenger variable into kernel_id
      kernel_id(iob,:) = msg_id(:)
    enddo
    deallocate( msg_id )
  


  end subroutine tag_clear_congestus_deep_clusters
  ! ===============================================================================




  ! ===============================================================================
  ! Subroutine to eliminate overly-small clusters before performing GMM-EnKF
  ! Designed for up to three clusters
  ! Requires kernel ids to be inserted into array msg_id to work
  subroutine remove_overly_small_clusters()
  
    implicit none

    integer, allocatable, dimension(:) :: cluster_size
    logical :: exit_flag

    ! Allocate cluster size
    allocate( cluster_size( 3 ) )

    ! Compute kernel sizes using IDs stored in msg_id
    do i = 1, 3
      cluster_size(i) = count( msg_id == i )
    enddo

    if ( my_proc_id == 0 .and. print_detail > 10 ) write(*,*) 'Raw cluster sizes:', cluster_size

    ! If clear cluster is too small, combine clear and shallow cloud clusters
    if ( cluster_size(1) < min_mems_per_kernel ) &
      where ( msg_id == 1 ) msg_id = 2

    ! If deep or shallow cloud clusters are too small, combine deep and shallow cloud
    ! clusters
    if ( ( cluster_size(3) < min_mems_per_kernel ) &
         .or. ( cluster_size(2) < min_mems_per_kernel ) ) &
      where ( msg_id == 3) msg_id = 2

    ! Recalculate the cluster sizes
    do i = 1, 3
      cluster_size(i) = count( msg_id == i )
    enddo

    ! At this point, if cluster 2 was originally too small, we have absorbed
    ! deep clouds (cluster 3) in. If cluster 2 is still too small, we will
    ! absorb the clear-sky cluster (cluster 2)  into cluster 2.
    if ( cluster_size(2) < min_mems_per_kernel ) &
      where( msg_id == 1 ) msg_id = 2

    ! Recalculate cluster sizes and check if cluster sizes are sane
    exit_flag=.false.
    do i = 1, 3
      ! measuring cluster size
      cluster_size(i) = count( msg_id == i )
      ! working on cluster checker.
      if (  ( cluster_size(i) < min_mems_per_kernel ) &
            .and. ( cluster_size(i) > 0 )  ) then
        if ( my_proc_id==0 ) &
          write(*,*) 'Cluster', i,' has inappropriate size:', cluster_size
        exit_flag=.true.
      endif
    enddo

    ! If cluster size was inappropriate, exit
    if (exit_flag) then
      write(*,*) 'Exiting'
      call exit(ierr)
    endif

    if ( my_proc_id == 0 .and. print_detail > 10) write(*,*) 'Adjusted cluster sizes:', cluster_size

  end subroutine remove_overly_small_clusters
  ! ====================================================================================




  ! ====================================================================================
  ! Subroutine run checks on GMM-EnKF statistics. These checks are:
  ! 1) Are the posterior weights within +- 1./numbers_en from theoretical values?
  ! 2) Are posterior cluster means nigh identical to the Kalman Filter values?
  ! 3) Is the posterior expanding cluster state identical to Kalman filter values?
  subroutine gmm_enkf_check()

    implicit none

    real, allocatable, dimension(:):: prior_state_means, poste_state_means
    real, allocatable, dimension(:):: prior_state_vars, poste_state_vars, c_xy
    real, allocatable, dimension(:):: thry_state_means, thry_state_vars
    real, allocatable, dimension(:):: thry_poste_weights
    logical :: weight_check, mean_check, var_check
    real :: weight_tolerance
    integer :: obs_i, obs_j, obs_k

    ! Allocate variables for sanity check
    allocate( prior_state_means ( max_kernel_num ) ); prior_state_means  = 0.
    allocate( poste_state_means ( max_kernel_num ) ); poste_state_means  = 0.
    allocate( prior_state_vars  ( max_kernel_num ) ); prior_state_vars   = 0.
    allocate( poste_state_vars  ( max_kernel_num ) ); poste_state_vars   = 0.
    allocate( c_xy              ( max_kernel_num ) ); c_xy               = 0.
    allocate( thry_state_means  ( max_kernel_num ) ); thry_state_means   = 0.
    allocate( thry_state_vars   ( max_kernel_num ) ); thry_state_vars    = 0.
    allocate( thry_poste_weights( max_kernel_num ) ); thry_poste_weights = 0.

    ! Tolerance in weights
    weight_tolerance = 1./numbers_en

    ! Figure out obs index position on the slab
    obs_i = nint( obs%position(iob,1) ) - istart+1
    obs_j = nint( obs%position(iob,2) ) - jstart+1
    obs_k = 10 !nint( obs%position(iob,3) )


    ! PART 1: Estimate ensemble statistics and calculate theoretical posterior values
    ! --------------------------------------------------------------------------------
    ! Compute prior state means, variances and covariances using backup_slab
    do cid = 1, max_kernel_num
      prior_state_means(cid) = sum( backup_slab(obs_i,obs_j,obs_k,1:numbers_en),    &
                                    mask = (prior_id == cid) )                      &
                                  / ( prior_weights(cid) * numbers_en )

      prior_state_vars(cid)   = sum( ( backup_slab(obs_i,obs_j,obs_k,1:numbers_en)   &
                                       - prior_state_means(cid) )**2,                &
                                     mask = (prior_id == cid) )                      &
                                   / ( prior_weights(cid) * numbers_en - 1 )

      c_xy(cid) = sum( ( backup_slab(obs_i,obs_j,obs_k,1:numbers_en)                &
                         - prior_state_means(cid) )                                 &
                       * ( ya(iob,1:numbers_en) - prior_cluster_ymean(cid)),        &
                       mask = (prior_id == cid) )                                   &
                     / ( prior_weights(cid) * numbers_en - 1 )
    enddo


    ! Compute theoretical posterior weight values
    thry_poste_weights = exp( (-0.5)*( prior_cluster_ymean - obs%dat(iob) )**2    &
                              / ( prior_cluster_yvariance + error**2 ) )          &
                         /sqrt( 2*3.1415*(prior_cluster_yvariance + error**2) )
    thry_poste_weights = thry_poste_weights * prior_weights
    thry_poste_weights = thry_poste_weights / sum( thry_poste_weights )


    ! Compute theoretical Kalman filter posterior mean for state.
    ! <x^a> = <x^f> + (C_xy/(var_f + var_o) ) * ( y - <Hx^f>)
    thry_state_means = prior_state_means                                          &
                       +( c_xy / ( prior_cluster_yvariance + error**2 ) )         &
                        *( obs%dat(iob) - prior_cluster_ymean )


    ! Compute theoretical Kalman filter posterior variance for state.
    ! var_a = var_f - (C_xy/(var_f + var_o) ) * C_xy
    thry_state_vars  = prior_state_vars                                           &
                       - c_xy**2 / ( prior_cluster_yvariance + error**2 )


    ! Compute ensemble posterior state means and variances
    do cid = 1, max_kernel_num
      poste_state_means(cid) = sum( x(obs_i,obs_j,obs_k,m,1:numbers_en),    &
                                    mask = (poste_id == cid) )              &
                                  / ( poste_weights(cid) * numbers_en )

      poste_state_vars(cid)   = sum( ( x(obs_i,obs_j,obs_k,m,1:numbers_en)   &
                                       - poste_state_means(cid) )**2,        &
                                     mask = (poste_id == cid) )              &
                                   / ( poste_weights(cid) * numbers_en - 1 )
    enddo



    ! PART 2: Check theoretical posterior values against posterior ensemble values
    ! ----------------------------------------------------------------------------

    ! Check 1: Are the posterior weights within tolerable range from theory?
    weight_check = .true.
    do cid = 1, max_kernel_num
      if ( abs( poste_weights(cid) - thry_poste_weights(cid) ) > weight_tolerance ) &
        weight_check = .false.
    enddo

    ! Check 2: Are the posterior means within tolerable range from theory?
    mean_check = .true.
    do cid = 1, max_kernel_num
      if ( abs( (poste_state_means(cid) - thry_state_means(cid)) &
                 / thry_state_means(cid) ) > 1e-3 ) &
        mean_check = .false.
    enddo

    ! Check 3: Is the expanding cluster standard deviation within tolerable range from theory?
    var_check = .true.
    do cid = 1, max_kernel_num
      ! Ignore non-expanding clusters
      if ( prior_weights(cid) >= poste_weights(cid) ) cycle

      ! Check standard deviations
      if ( abs( sqrt(poste_state_vars(cid)) - sqrt(thry_state_vars(cid)) ) &
            / sqrt(thry_state_vars(cid) ) > 1e-3  ) &
        var_check = .false.
    enddo



    ! PART 3: Print out sanity check information
    ! ------------------------------------------
    ! Set up file to output stuff
    open( unit=19, file = 'gmm_enkf_check.log', status='old', position='append', &
          action='write')

    ! Write observation information
    write(19, '(/,/,a,x,i7,x,a,2(a,x,f8.2,x),13(/,a,x,f,x,f),/,3(a,x,L,x))') &
      'Obs', iob+batch_start_obs_num, obstype, ', obs val=', obs%dat(iob), &
      ', obs err var=', error**2, &
      '    Prior count:         ', prior_weights * numbers_en, &
      '    Prior weight:        ', prior_weights, &
      '    Prior obs means:     ', prior_cluster_ymean, &
      '    Prior obs variances: ', prior_cluster_yvariance, &
      '    Poste weights (ensemble):   ', poste_weights, &
      '    Poste weights (theoretical):', thry_poste_weights, &
      '    Prior state means:              ', prior_state_means, &
      '    Prior state vars:               ' , prior_state_vars, &
      '    Prior state covs:               ' , c_xy,             &
      '    Poste state means (ensemble):   ', poste_state_means, &
      '    Poste state means (theoretical):', thry_state_means, &
      '    Poste state std dev (ensemble):   ', sqrt(poste_state_vars), &
      '    Poste state std dev (theoretical):', sqrt(thry_state_vars), &
      '    Weights check: ', weight_check, &
      '    Mean check:    ', mean_check, &
      '    Std dev check: ', var_check

    ! Close file
    close(19)

    ! Deallocate variables
    deallocate( prior_state_means  )
    deallocate( poste_state_means  )
    deallocate( prior_state_vars   )
    deallocate( poste_state_vars   )
    deallocate( c_xy               )



  end subroutine gmm_enkf_check
  ! ====================================================================================




end module gauss_mix_model_enkf
