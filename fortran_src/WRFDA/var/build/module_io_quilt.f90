



































MODULE module_wrf_quilt



















































  USE module_internal_header_util
  USE module_timing

  INTEGER, PARAMETER :: int_num_handles = 99
  INTEGER, PARAMETER :: max_servers = int_num_handles+1  
  LOGICAL, DIMENSION(0:int_num_handles) :: okay_to_write, int_handle_in_use, okay_to_commit
  INTEGER, DIMENSION(0:int_num_handles) :: int_num_bytes_to_write, io_form
  REAL, POINTER,SAVE :: int_local_output_buffer(:)
  INTEGER,      SAVE :: int_local_output_cursor
  LOGICAL          :: quilting_enabled
  LOGICAL          :: disable_quilt = .FALSE.
  INTEGER          :: prev_server_for_handle = -1
  INTEGER          :: server_for_handle(int_num_handles)
  INTEGER          :: reduced(2), reduced_dummy(2)
  LOGICAL, EXTERNAL :: wrf_dm_on_monitor

  INTEGER :: mpi_comm_avail,availrank
  LOGICAL :: in_avail=.false., poll_servers=.false.

  INTEGER nio_groups
  INTEGER :: mpi_comm_local
  LOGICAL :: compute_node
  LOGICAL :: compute_group_master(max_servers)
  INTEGER :: mpi_comm_io_groups(max_servers)
  INTEGER :: nio_tasks_in_group
  INTEGER :: nio_tasks_per_group
  INTEGER :: ncompute_tasks
  INTEGER :: ntasks
  INTEGER :: mytask

  INTEGER, PARAMETER           :: onebyte = 1
  INTEGER comm_io_servers, iserver, hdrbufsize, obufsize
  INTEGER, DIMENSION(4096)     :: hdrbuf
  INTEGER, DIMENSION(int_num_handles)     :: handle


  CONTAINS

    INTEGER FUNCTION get_server_id ( dhandle )









      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dhandle
      IF ( dhandle .GE. 1 .AND. dhandle .LE. int_num_handles ) THEN
        IF ( server_for_handle ( dhandle ) .GE. 1 ) THEN
          get_server_id = server_for_handle ( dhandle )
        ELSE
           IF(poll_servers) THEN
              
              call wrf_quilt_find_server(server_for_handle(dhandle))
           ELSE
              
              prev_server_for_handle = mod ( prev_server_for_handle + 1 , nio_groups )
              server_for_handle( dhandle ) = prev_server_for_handle+1
           ENDIF
           get_server_id=server_for_handle(dhandle)
        ENDIF
      ELSE
         CALL wrf_message('module_io_quilt: get_server_id bad dhandle' )
      ENDIF
    END FUNCTION get_server_id

    SUBROUTINE set_server_id ( dhandle, value )
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: dhandle, value
       IF ( dhandle .GE. 1 .AND. dhandle .LE. int_num_handles ) THEN
         server_for_handle(dhandle) = value
       ELSE
         CALL wrf_message('module_io_quilt: set_server_id bad dhandle' )
       ENDIF
    END SUBROUTINE set_server_id

    LOGICAL FUNCTION get_poll_servers() 
      implicit none
      get_poll_servers=poll_servers
    end FUNCTION get_poll_servers

    SUBROUTINE int_get_fresh_handle( retval )










      INTEGER i, retval
      retval = -1
      DO i = 1, int_num_handles
        IF ( .NOT. int_handle_in_use(i) )  THEN
          retval = i
          GOTO 33
        ENDIF
      ENDDO
33    CONTINUE
      IF ( retval < 0 )  THEN
        CALL wrf_error_fatal3("<stdin>",181,&
"frame/module_io_quilt.F: int_get_fresh_handle() can not")
      ENDIF
      int_handle_in_use(i) = .TRUE.
      NULLIFY ( int_local_output_buffer )
    END SUBROUTINE int_get_fresh_handle

    SUBROUTINE setup_quilt_servers ( nio_tasks_per_group,     &
                                     mytask,                  &
                                     ntasks,                  &
                                     nproc_x,                 &
                                     nproc_y,                 &
                                     n_groups_arg,            &
                                     nio,                     &
                                     mpi_comm_wrld,           &
                                     mpi_comm_local,          &
                                     mpi_comm_io_groups)




























































      USE module_configure
      USE module_dm, ONLY : compute_mesh
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER,                      INTENT(IN)  :: nio_tasks_per_group, mytask, ntasks, &
                                                   n_groups_arg, mpi_comm_wrld
      INTEGER,                      INTENT(IN)  :: nproc_x, nproc_y
      INTEGER,  INTENT(OUT)                     :: mpi_comm_local, nio
      INTEGER, DIMENSION(100),      INTENT(OUT) :: mpi_comm_io_groups

      INTEGER                     :: i, j, ii, comdup, ierr, niotasks, n_groups, iisize
      INTEGER, DIMENSION(ntasks)  :: icolor
      CHARACTER*128 mess
      INTEGER :: io_form_setting
      INTEGER :: me
      INTEGER :: k, m, nprocx, nprocy
      LOGICAL :: reorder_mesh



      CALL nl_get_io_form_history(1,   io_form_setting) ; call sokay( 'history', io_form_setting )
      CALL nl_get_io_form_restart(1,   io_form_setting) ; call sokay( 'restart', io_form_setting )
      CALL nl_get_io_form_auxhist1(1,  io_form_setting) ; call sokay( 'auxhist1', io_form_setting )
      CALL nl_get_io_form_auxhist2(1,  io_form_setting) ; call sokay( 'auxhist2', io_form_setting )
      CALL nl_get_io_form_auxhist3(1,  io_form_setting) ; call sokay( 'auxhist3', io_form_setting )
      CALL nl_get_io_form_auxhist4(1,  io_form_setting) ; call sokay( 'auxhist4', io_form_setting )
      CALL nl_get_io_form_auxhist5(1,  io_form_setting) ; call sokay( 'auxhist5', io_form_setting )
      CALL nl_get_io_form_auxhist6(1,  io_form_setting) ; call sokay( 'auxhist6', io_form_setting )
      CALL nl_get_io_form_auxhist7(1,  io_form_setting) ; call sokay( 'auxhist7', io_form_setting )
      CALL nl_get_io_form_auxhist8(1,  io_form_setting) ; call sokay( 'auxhist8', io_form_setting )
      CALL nl_get_io_form_auxhist9(1,  io_form_setting) ; call sokay( 'auxhist9', io_form_setting )
      CALL nl_get_io_form_auxhist10(1, io_form_setting) ; call sokay( 'auxhist10', io_form_setting )
      CALL nl_get_io_form_auxhist11(1, io_form_setting) ; call sokay( 'auxhist11', io_form_setting )

      n_groups = n_groups_arg
      IF ( n_groups .LT. 1 ) n_groups = 1

      compute_node = .TRUE.







      nio = nio_tasks_per_group
      ncompute_tasks = ntasks - (nio * n_groups)
      IF ( ncompute_tasks .LT. nio ) THEN 
        WRITE(mess,'("Not enough tasks to have ",I3," groups of ",I3," I/O tasks. No quilting.")')n_groups,nio
        nio            = 0
        ncompute_tasks = ntasks
      ELSE                                   
        WRITE(mess,'("Quilting with ",I3," groups of ",I3," I/O tasks.")')n_groups,nio
      ENDIF                                   
      CALL wrf_message(mess)

      IF ( nio .LT. 0 ) THEN
        nio = 0
      ENDIF
      IF ( nio .EQ. 0 ) THEN
        quilting_enabled = .FALSE.
        mpi_comm_local = mpi_comm_wrld
        mpi_comm_io_groups = mpi_comm_wrld
        RETURN
      ENDIF
      quilting_enabled = .TRUE.



      DO i = 1, ncompute_tasks
        icolor(i) = 0
      ENDDO
      ii = 1

      DO i = ncompute_tasks+1, ntasks, nio
        DO j = i, i+nio-1
          icolor(j) = ii
        ENDDO
        ii = ii+1
      ENDDO
      CALL MPI_Comm_dup(mpi_comm_wrld,comdup,ierr)
      CALL MPI_Comm_split(comdup,icolor(mytask+1),mytask,mpi_comm_local,ierr)


      CALL nl_get_reorder_mesh(1,reorder_mesh)
      IF ( reorder_mesh ) THEN
        reorder_mesh = .FALSE.
        CALL nl_set_reorder_mesh(1,reorder_mesh)
        CALL wrf_message('Warning: reorder_mesh does not work with quilting. Disabled reorder_mesh.')
      ENDIF
      
      IF ( nproc_x .NE. -1 .AND. nproc_y .NE. -1 ) THEN
       nprocx=nproc_x
       nprocy=nproc_y
      ELSE
       CALL compute_mesh( ncompute_tasks, nprocx, nprocy )
      ENDIF

      nio = min(nio,nprocy)
      m = mod(nprocy,nio)  
      ii = 1
      DO j = 1, nio, 1
         DO k = 1,nprocy/nio+min(m,1)
           DO i = 1, nprocx
             icolor(ii) = j - 1
             ii = ii + 1
           ENDDO
         ENDDO
         m = max(m-1,0)
      ENDDO


      DO j = 1, n_groups
        
        DO i = ncompute_tasks+1,ntasks
          icolor(i) = MPI_UNDEFINED
        ENDDO
        ii = 0
        DO i = ncompute_tasks+(j-1)*nio+1,ncompute_tasks+j*nio
          icolor(i) = ii
          ii = ii+1
        ENDDO
        CALL MPI_Comm_dup(mpi_comm_wrld,comdup,ierr)
        CALL MPI_Comm_split(comdup,icolor(mytask+1),mytask, &
                            mpi_comm_io_groups(j),ierr)
      ENDDO

         if(nio_groups==1) then
            poll_servers=.false.
            call wrf_message('Server polling is does not work with one io group.  Disabled poll_servers.')
         endif

      if(poll_servers) then
         
         
         
         
         

         

         call mpi_comm_rank(mpi_comm_wrld,me,ierr)

         icolor=MPI_UNDEFINED
         in_avail=.false.

         if(wrf_dm_on_monitor()) then
            in_avail=.true. 
         endif
         icolor(1)=1

         do j=1,n_groups
            i=ncompute_tasks+j*nio-1
            if(me+1==i) then
               in_avail=.true. 
            endif
            icolor(i)=1
         enddo

         CALL MPI_Comm_dup(mpi_comm_wrld,comdup,ierr)
         CALL MPI_Comm_split(comdup,icolor(me+1),me, &
                             mpi_comm_avail,ierr)

         availrank=MPI_UNDEFINED
         if(in_avail) then
            call mpi_comm_rank(mpi_comm_avail,availrank,ierr)
         endif

      endif

      compute_group_master = .FALSE.
      compute_node         = .FALSE.

      DO j = 1, n_groups

         IF ( mytask .LT. ncompute_tasks .OR.                                                  &    
              (ncompute_tasks+(j-1)*nio .LE. mytask .AND. mytask .LT. ncompute_tasks+j*nio) &    
            ) THEN

         CALL MPI_Comm_Size( mpi_comm_io_groups(j) , iisize, ierr )
         
         
         CALL MPI_Comm_Rank( mpi_comm_io_groups(j) , me , ierr )

         
         
         
         IF (ncompute_tasks+(j-1)*nio .LE. mytask .AND. mytask .LT. ncompute_tasks+j*nio) THEN
            mpi_comm_io_groups(1) = mpi_comm_io_groups(j)
         ELSE
            compute_node = .TRUE.
            
            
            
            compute_group_master(j) = (me .EQ. 0)


         ENDIF
         ENDIF
      ENDDO

    END SUBROUTINE setup_quilt_servers

    SUBROUTINE sokay ( stream, io_form )
    USE module_state_description
    CHARACTER*(*) stream
    CHARACTER*256 mess
    INTEGER io_form

    SELECT CASE (io_form)
      CASE ( IO_NETCDF   )
         RETURN
      CASE ( IO_INTIO   )
         RETURN
      CASE ( IO_GRIB1 )
         RETURN
      CASE (0)
         RETURN
      CASE DEFAULT
         WRITE(mess,*)' An output format has been specified that is incompatible with quilting: io_form: ',io_form,' ',TRIM(stream)
         CALL wrf_error_fatal3("<stdin>",478,&
mess)
    END SELECT
    END SUBROUTINE sokay

    SUBROUTINE quilt













      USE module_state_description
      USE module_quilt_outbuf_ops
      USE module_configure, only : grid_config_rec_type, model_config_rec, model_to_grid_config_rec
      IMPLICIT NONE
      INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101
      integer, parameter  :: WRF_FILE_OPENED_FOR_WRITE            = 102
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103
      integer, parameter  :: WRF_REAL                             = 105
      integer, parameter  :: WRF_DOUBLE                           = 105
      integer, parameter  :: WRF_FLOAT=WRF_REAL
      integer, parameter  :: WRF_INTEGER                          = 106
      integer, parameter  :: WRF_LOGICAL                          = 107
      integer, parameter  :: WRF_COMPLEX                          = 108
      integer, parameter  :: WRF_DOUBLE_COMPLEX                   = 109
      integer, parameter  :: WRF_FILE_OPENED_FOR_UPDATE           = 110


      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102
      TYPE (grid_config_rec_type)  :: config_flags
      INTEGER itag, ninbuf, ntasks_io_group, ntasks_local_group, mytask_local, ierr
      INTEGER istat
      INTEGER mytask_io_group
      INTEGER   :: nout_set = 0
      INTEGER   :: obufsize, bigbufsize, chunksize, sz
      REAL, DIMENSION(1)      :: dummy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: obuf, bigbuf
      REAL,    ALLOCATABLE, DIMENSION(:) :: RDATA
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDATA
      CHARACTER (LEN=512) :: CDATA
      CHARACTER (LEN=80) :: fname
      INTEGER icurs, hdrbufsize, itypesize, ftypesize, rtypesize, Status, fstat, io_form_arg
      INTEGER :: DataHandle, FieldType, Comm, IOComm, DomainDesc, code, Count
      INTEGER, DIMENSION(3) :: DomainStart , DomainEnd , MemoryStart , MemoryEnd , PatchStart , PatchEnd
      INTEGER :: dummybuf(1)
      INTEGER :: num_noops, num_commit_messages, num_field_training_msgs, hdr_tag
      CHARACTER (len=256) :: DateStr , Element, VarName, MemoryOrder , Stagger , DimNames(3), FileName, SysDepInfo, mess
      INTEGER, EXTERNAL :: use_package
      LOGICAL           :: stored_write_record, retval
      INTEGER iii, jjj, vid, CC, DD, dom_id
      LOGICAL           :: call_server_ready

logical okay_to_w
character*120 sysline

      dom_id = 1 
      CALL model_to_grid_config_rec ( dom_id , model_config_rec , config_flags )









      SysDepInfo = " "
      if ( config_flags%use_netcdf_classic ) SysDepInfo="use_netcdf_classic"
      CALL ext_ncd_ioinit( SysDepInfo, ierr )
      SysDepInfo = " "
      CALL ext_int_ioinit( SysDepInfo, ierr )
      CALL ext_gr1_ioinit( SysDepInfo, ierr)

      call_server_ready = .true. 

      okay_to_commit = .false.
      stored_write_record = .false.
      ninbuf = 0
      
      
      
      
      
      CALL Mpi_Comm_Size (  mpi_comm_io_groups(1), ntasks_io_group,    ierr  )
      CALL MPI_COMM_RANK( mpi_comm_io_groups(1), mytask_io_group,    ierr )
      CALL Mpi_Comm_Size (  mpi_comm_local,        ntasks_local_group, ierr  )
      CALL MPI_COMM_RANK( mpi_comm_local,        mytask_local,       ierr )

      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
      IF ( itypesize <= 0 ) THEN
        CALL wrf_error_fatal3("<stdin>",612,&
"external/RSL/module_dm.F: quilt: type size <= 0 invalid")
      ENDIF







       CC = ntasks_io_group - 1

       DD = ncompute_tasks / ntasks_local_group








okay_to_w = .false.
      DO WHILE (.TRUE.)  














         if(poll_servers .and. call_server_ready) then
            call_server_ready=.false.
            
            
            call wrf_quilt_server_ready()
         endif

        
        

        
        ! if needed (currently needed only for ioclose).
        reduced_dummy = 0
        CALL MPI_Reduce( reduced_dummy, reduced, 2, MPI_INTEGER, MPI_SUM, mytask_io_group, mpi_comm_io_groups(1), ierr )
        obufsize = reduced(1)



        IF ( obufsize .LT. 0 ) THEN
          IF ( obufsize .EQ. -100 ) THEN         
            CALL ext_ncd_ioexit( Status )
            CALL ext_int_ioexit( Status )
            CALL ext_gr1_ioexit( Status )
            CALL wrf_message ( 'I/O QUILT SERVERS DONE' )
               CALL mpi_finalize(ierr)
            STOP
          ELSE
            WRITE(mess,*)'Possible 32-bit overflow on output server. Try larger nio_tasks_per_group in namelist.'
            CALL wrf_error_fatal3("<stdin>",677,&
mess)
          ENDIF
        ENDIF







        IF ( obufsize .GT. 0 ) THEN
          ALLOCATE( obuf( (obufsize+1)/itypesize ) )


          CALL collect_on_comm_debug("module_io_quilt_old.F",722, mpi_comm_io_groups(1),        &
                                onebyte,                      &
                                dummy, 0,                     &
                                obuf, obufsize )

        ELSE
          ! Necessarily, the compute processes send the ioclose signal,
          
          ! will stall on the ioclose message waiting for the quilt 
          
          
          
          
          
          ! Then a header representing the ioclose message is constructed
          
          
          
          
          
          
          ALLOCATE( obuf( 4096 ) )
          
          CALL int_gen_handle_header( obuf, obufsize, itypesize, &
                                      reduced(2) , int_ioclose )

          if(poll_servers) then 
             
             
             call_server_ready=.true.
          endif
        ENDIF











        CALL init_store_piece_of_field
        CALL mpi_type_size ( MPI_INTEGER , itypesize , ierr )



        vid = 0
        icurs = itypesize
        num_noops = 0 
        num_commit_messages = 0 
        num_field_training_msgs = 0 
        DO WHILE ( icurs .lt. obufsize ) 
          hdr_tag = get_hdr_tag( obuf ( icurs / itypesize ) )
          SELECT CASE ( hdr_tag )
            CASE ( int_field )
              CALL int_get_write_field_header ( obuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                DomainStart , DomainEnd ,                                    &
                                                MemoryStart , MemoryEnd ,                                    &
                                                PatchStart , PatchEnd )
              chunksize = (PatchEnd(1)-PatchStart(1)+1)*(PatchEnd(2)-PatchStart(2)+1)* &
                          (PatchEnd(3)-PatchStart(3)+1)*ftypesize

              IF ( DomainDesc .EQ. 333933 ) THEN  
                 IF ( num_field_training_msgs .EQ. 0 ) THEN
                   call add_to_bufsize_for_field( VarName, hdrbufsize )

                 ENDIF
                 num_field_training_msgs = num_field_training_msgs + 1
              ELSE
                 call add_to_bufsize_for_field( VarName, hdrbufsize )

              ENDIF
              icurs = icurs + hdrbufsize



              
              
              IF ( DomainDesc .NE. 333933 ) THEN   

                call add_to_bufsize_for_field( VarName, chunksize )
                icurs = icurs + chunksize
              ENDIF
            CASE ( int_open_for_write_commit )  
              hdrbufsize = obuf(icurs/itypesize)
              IF (num_commit_messages.EQ.0) THEN
                call add_to_bufsize_for_field( 'COMMIT', hdrbufsize )
              ENDIF
              num_commit_messages = num_commit_messages + 1
              icurs = icurs + hdrbufsize
            CASE DEFAULT
              hdrbufsize = obuf(icurs/itypesize)




































              IF ((hdr_tag.EQ.int_noop.AND.mytask_local.NE.0.AND.num_noops.LE.0)  &
                  .OR.hdr_tag.NE.int_noop) THEN
                write(VarName,'(I5.5)')vid 

                call add_to_bufsize_for_field( VarName, hdrbufsize )
                vid = vid+1
              ENDIF
              IF ( hdr_tag .EQ. int_noop ) num_noops = num_noops + 1
              icurs = icurs + hdrbufsize
          END SELECT
        ENDDO 



        vid = 0
        icurs = itypesize
        num_noops = 0 
        num_commit_messages = 0 
        num_field_training_msgs = 0 
        DO WHILE ( icurs .lt. obufsize ) 

          hdr_tag = get_hdr_tag( obuf ( icurs / itypesize ) )
          SELECT CASE ( hdr_tag )
            CASE ( int_field )
              CALL int_get_write_field_header ( obuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                DomainStart , DomainEnd ,                                    &
                                                MemoryStart , MemoryEnd ,                                    &
                                                PatchStart , PatchEnd )
              chunksize = (PatchEnd(1)-PatchStart(1)+1)*(PatchEnd(2)-PatchStart(2)+1)* &
                          (PatchEnd(3)-PatchStart(3)+1)*ftypesize

              IF ( DomainDesc .EQ. 333933 ) THEN  
                 IF ( num_field_training_msgs .EQ. 0 ) THEN
                   call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )

                 ENDIF
                 num_field_training_msgs = num_field_training_msgs + 1
              ELSE
                 call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )

              ENDIF
              icurs = icurs + hdrbufsize
              
              
              IF ( DomainDesc .NE. 333933 ) THEN   

                call store_piece_of_field( obuf(icurs/itypesize), VarName, chunksize )
                icurs = icurs + chunksize
              ENDIF
            CASE ( int_open_for_write_commit )  
              hdrbufsize = obuf(icurs/itypesize)
              IF (num_commit_messages.EQ.0) THEN
                call store_piece_of_field( obuf(icurs/itypesize), 'COMMIT', hdrbufsize )
              ENDIF
              num_commit_messages = num_commit_messages + 1
              icurs = icurs + hdrbufsize
            CASE DEFAULT
              hdrbufsize = obuf(icurs/itypesize)
              IF ((hdr_tag.EQ.int_noop.AND.mytask_local.NE.0.AND.num_noops.LE.0)  &
                  .OR.hdr_tag.NE.int_noop) THEN
                write(VarName,'(I5.5)')vid 

                call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )
                vid = vid+1
              ENDIF
              IF ( hdr_tag .EQ. int_noop ) num_noops = num_noops + 1
              icurs = icurs + hdrbufsize
          END SELECT
        ENDDO 



        CALL init_retrieve_pieces_of_field


        CALL retrieve_pieces_of_field ( obuf , VarName, obufsize, sz, retval )


        CALL MPI_Reduce( sz, bigbufsize, 1, MPI_INTEGER, MPI_SUM, ntasks_local_group-1, mpi_comm_local, ierr )



        DO WHILE ( retval ) 



          IF ( mytask_local .EQ. ntasks_local_group-1 ) THEN
            ALLOCATE( bigbuf( (bigbufsize+1)/itypesize ) )
         else
            ALLOCATE( bigbuf(1) )
          ENDIF



          CALL collect_on_comm_debug2("module_io_quilt_old.F",952,Trim(VarName),        &
                                get_hdr_tag(obuf),sz,get_hdr_rec_size(obuf),  &
                                mpi_comm_local,                               &
                                onebyte,                                      &
                                obuf, sz,                                     &
                                bigbuf, bigbufsize )


          IF ( mytask_local .EQ. ntasks_local_group-1 ) THEN






            icurs = itypesize  

            stored_write_record = .false.


            DO WHILE ( icurs .lt. bigbufsize ) 
              CALL mpi_type_size ( MPI_INTEGER , itypesize , ierr )





              SELECT CASE ( get_hdr_tag( bigbuf(icurs/itypesize) ) )




                CASE ( int_noop )
                  CALL int_get_noop_header( bigbuf(icurs/itypesize), hdrbufsize, itypesize )
                  icurs = icurs + hdrbufsize


                CASE ( int_dom_td_real )
                  CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                  ALLOCATE( RData( bigbuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_td_header( bigbuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, DateStr, Element, RData, Count, code )
                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                     CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( RData )

                CASE ( int_dom_ti_real )

                  CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                  ALLOCATE( RData( bigbuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_ti_header( bigbuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, Element, RData, Count, code )
                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )

                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( RData )


                CASE ( int_dom_td_integer )

                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  ALLOCATE( IData( bigbuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_td_header( bigbuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, DateStr, Element, IData, Count, code )
                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( IData )


                CASE ( int_dom_ti_integer )


                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  ALLOCATE( IData( bigbuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_ti_header( bigbuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, Element, IData, Count, code )
                  icurs = icurs + hdrbufsize
                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )

                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )

                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( IData)
 

                CASE ( int_set_time )

                  CALL int_get_ti_header_char( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )
                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_INTIO   )
                      CALL ext_int_set_time ( handle(DataHandle), TRIM(CData), Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  icurs = icurs + hdrbufsize


                CASE ( int_dom_ti_char )

                  CALL int_get_ti_header_char( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )


                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  icurs = icurs + hdrbufsize


                CASE ( int_var_ti_char )

                  CALL int_get_ti_header_char( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  icurs = icurs + hdrbufsize

                CASE ( int_ioexit )

                  CALL wrf_error_fatal3("<stdin>",1100,&
                         "quilt: should have handled int_ioexit already")
! The I/O server "root" handles the "ioclose" request.
                CASE ( int_ioclose )
                  CALL int_get_handle_header( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                              DataHandle , code )
                  icurs = icurs + hdrbufsize

                  IF ( DataHandle .GE. 1 ) THEN

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_inquire_filename( handle(DataHandle), fname, fstat, Status )
                      IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                        CALL ext_ncd_ioclose(handle(DataHandle),Status)
                      ENDIF
                    CASE ( IO_INTIO   )
                      CALL ext_int_inquire_filename( handle(DataHandle), fname, fstat, Status )
                      IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                        CALL ext_int_ioclose(handle(DataHandle),Status)
                      ENDIF
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_inquire_filename( handle(DataHandle), fname, fstat, Status )
                    IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                      CALL ext_gr1_ioclose(handle(DataHandle),Status)
                    ENDIF
                    CASE DEFAULT
                      Status = 0
                  END SELECT
                  ENDIF



                  IF (fname(1:6) .EQ. 'wrfout' .AND. config_flags%output_ready_flag ) THEN
                    OPEN (unit=99,file='wrfoutReady' // fname(7:30), status='unknown', access='sequential')
                    CLOSE (99)
                  ENDIF


                CASE ( int_open_for_write_begin )

                  CALL int_get_ofwb_header( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                            FileName,SysDepInfo,io_form_arg,DataHandle )





                  icurs = icurs + hdrbufsize

                
                  io_form(DataHandle) = io_form_arg

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)

                    CASE ( IO_INTIO   )
                      CALL ext_int_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)
                    CASE ( IO_GRIB1 )
                       CALL ext_gr1_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT
                
                  okay_to_write(DataHandle) = .false.





                CASE ( int_open_for_write_commit )

                  CALL int_get_handle_header( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                              DataHandle , code )
                  icurs = icurs + hdrbufsize
                  okay_to_commit(DataHandle) = .true.











                CASE ( int_field )
                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  CALL int_get_write_field_header ( bigbuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                    DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                    DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                    DomainStart , DomainEnd ,                                    &
                                                    MemoryStart , MemoryEnd ,                                    &
                                                    PatchStart , PatchEnd )

                  icurs = icurs + hdrbufsize

                  IF ( okay_to_write(DataHandle) ) THEN




                    IF ( FieldType .EQ. WRF_FLOAT .OR. FieldType .EQ. WRF_DOUBLE)  THEN
                      
                      
                      IF ( FieldType .EQ. WRF_DOUBLE)  THEN

                        CALL mpi_type_size( MPI_DOUBLE_PRECISION, ftypesize, ierr )
                      ELSE
                        CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                      ENDIF
                      stored_write_record = .true.
                      CALL store_patch_in_outbuf ( bigbuf(icurs/itypesize), dummybuf, TRIM(DateStr), TRIM(VarName) , &
                                                   FieldType, TRIM(MemoryOrder), TRIM(Stagger), DimNames, &
                                                   DomainStart , DomainEnd , &
                                                   MemoryStart , MemoryEnd , &
                                                   PatchStart , PatchEnd )

                    ELSE IF ( FieldType .EQ. WRF_INTEGER ) THEN
                      CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                      stored_write_record = .true.
                      CALL store_patch_in_outbuf ( dummybuf, bigbuf(icurs/itypesize), TRIM(DateStr), TRIM(VarName) , &
                                                   FieldType, TRIM(MemoryOrder), TRIM(Stagger), DimNames, &
                                                   DomainStart , DomainEnd , &
                                                   MemoryStart , MemoryEnd , &
                                                   PatchStart , PatchEnd )
                    ELSE IF ( FieldType .EQ. WRF_LOGICAL ) THEN
                      ftypesize = 4
                    ENDIF
                    icurs = icurs + (PatchEnd(1)-PatchStart(1)+1)*(PatchEnd(2)-PatchStart(2)+1)* &
                                    (PatchEnd(3)-PatchStart(3)+1)*ftypesize
                  ELSE
                    SELECT CASE (use_package(io_form(DataHandle)))
                      CASE ( IO_NETCDF   )
                        CALL ext_ncd_write_field ( handle(DataHandle) , TRIM(DateStr) ,         &
                                   TRIM(VarName) , dummy , FieldType , Comm , IOComm,           &
                                   DomainDesc , TRIM(MemoryOrder) , TRIM(Stagger) , DimNames ,  &
                                   DomainStart , DomainEnd ,                                    &
                                   DomainStart , DomainEnd ,                                    &
                                   DomainStart , DomainEnd ,                                    &
                                   Status )
                      CASE DEFAULT
                        Status = 0
                    END SELECT
                  ENDIF
                CASE ( int_iosync )
                  CALL int_get_handle_header( bigbuf(icurs/itypesize), hdrbufsize, itypesize, &
                                            DataHandle , code )
                  icurs = icurs + hdrbufsize
                CASE DEFAULT
                  WRITE(mess,*)'quilt: bad tag: ',get_hdr_tag( bigbuf(icurs/itypesize) ),' icurs ',icurs/itypesize
                  CALL wrf_error_fatal3("<stdin>",1253,&
mess )
              END SELECT

            ENDDO 



            IF (stored_write_record) THEN







              CALL write_outbuf ( handle(DataHandle), use_package(io_form(DataHandle))) 

            ENDIF




            IF (okay_to_commit(DataHandle)) THEN

              SELECT CASE (use_package(io_form(DataHandle)))
                CASE ( IO_NETCDF   )
                  CALL ext_ncd_inquire_filename( handle(DataHandle), fname, fstat, Status )
                  IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                    CALL ext_ncd_open_for_write_commit(handle(DataHandle),Status)
                    okay_to_write(DataHandle) = .true.
                  ENDIF
                CASE ( IO_INTIO   )
                  CALL ext_int_inquire_filename( handle(DataHandle), fname, fstat, Status )
                  IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                    CALL ext_int_open_for_write_commit(handle(DataHandle),Status)
                    okay_to_write(DataHandle) = .true.
                  ENDIF
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_inquire_filename( handle(DataHandle), fname, fstat, Status )
                    IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                       CALL ext_gr1_open_for_write_commit(handle(DataHandle),Status)
                       okay_to_write(DataHandle) = .true.
                    ENDIF

                CASE DEFAULT
                  Status = 0
              END SELECT

            okay_to_commit(DataHandle) = .false.
          ENDIF
          DEALLOCATE( bigbuf )
        ENDIF
        if(allocated(bigbuf)) deallocate(bigbuf)


        CALL retrieve_pieces_of_field ( obuf , VarName, obufsize, sz, retval )


        CALL MPI_Reduce( sz, bigbufsize, 1, MPI_INTEGER,MPI_SUM, ntasks_local_group-1,mpi_comm_local, ierr )



      END DO 

      DEALLOCATE( obuf )

      
      IF (stored_write_record) THEN

        SELECT CASE ( use_package(io_form) )
          CASE ( IO_NETCDF   )
            CALL ext_ncd_iosync( handle(DataHandle), Status )
          CASE ( IO_GRIB1   )
            CALL ext_gr1_iosync( handle(DataHandle), Status )
          CASE ( IO_INTIO   )
            CALL ext_int_iosync( handle(DataHandle), Status )
          CASE DEFAULT
            Status = 0
        END SELECT

      ENDIF

      END DO 

    END SUBROUTINE quilt

    SUBROUTINE quilt_pnc





      USE module_state_description
      USE module_quilt_outbuf_ops
      IMPLICIT NONE
      INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101
      integer, parameter  :: WRF_FILE_OPENED_FOR_WRITE            = 102
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103
      integer, parameter  :: WRF_REAL                             = 105
      integer, parameter  :: WRF_DOUBLE                           = 105
      integer, parameter  :: WRF_FLOAT=WRF_REAL
      integer, parameter  :: WRF_INTEGER                          = 106
      integer, parameter  :: WRF_LOGICAL                          = 107
      integer, parameter  :: WRF_COMPLEX                          = 108
      integer, parameter  :: WRF_DOUBLE_COMPLEX                   = 109
      integer, parameter  :: WRF_FILE_OPENED_FOR_UPDATE           = 110


      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102
      INTEGER itag, ninbuf, ntasks_io_group, ntasks_local_group, mytask_local, ierr
      INTEGER istat
      INTEGER mytask_io_group
      INTEGER   :: nout_set = 0
      INTEGER   :: obufsize, bigbufsize, chunksize, sz
      REAL,                 DIMENSION(1) :: dummy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: obuf, bigbuf
      REAL,    ALLOCATABLE, DIMENSION(:) :: RDATA
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDATA
      CHARACTER (LEN=512) :: CDATA
      CHARACTER (LEN=80) :: fname
      INTEGER icurs, hdrbufsize, itypesize, ftypesize, rtypesize, Status, fstat, io_form_arg
      INTEGER :: DataHandle, FieldType, Comm, IOComm, DomainDesc, code, Count
      INTEGER, DIMENSION(3) :: DomainStart , DomainEnd , MemoryStart , MemoryEnd , PatchStart , PatchEnd
      INTEGER :: dummybuf(1)
      INTEGER :: num_noops, num_commit_messages, num_field_training_msgs, hdr_tag
      CHARACTER (len=256) :: DateStr , Element, VarName, MemoryOrder , Stagger , DimNames(3), FileName, SysDepInfo, mess
      INTEGER, EXTERNAL :: use_package
      LOGICAL           :: stored_write_record, retval, written_record
      INTEGER iii, jjj, vid, CC, DD





      SysDepInfo = " "
      CALL ext_ncd_ioinit( SysDepInfo, ierr)
      CALL ext_int_ioinit( SysDepInfo, ierr )
      CALL ext_gr1_ioinit( SysDepInfo, ierr)

      okay_to_commit = .false.
      stored_write_record = .false.
      ninbuf = 0
      
      
      CALL Mpi_Comm_Size (  mpi_comm_io_groups(1), ntasks_io_group,    ierr  )
      CALL MPI_COMM_RANK( mpi_comm_io_groups(1), mytask_io_group,    ierr )
      CALL Mpi_Comm_Size (  mpi_comm_local,        ntasks_local_group, ierr  )
      CALL MPI_COMM_RANK( mpi_comm_local,        mytask_local,       ierr )

      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
      IF ( itypesize <= 0 ) THEN
        CALL wrf_error_fatal3("<stdin>",1441,&
"external/RSL/module_dm.F: quilt: type size <= 0 invalid")
      ENDIF







       CC = ntasks_io_group - 1

       DD = ncompute_tasks / ntasks_local_group









      DO WHILE (.TRUE.)  













        
        

        
        ! if needed (currently needed only for ioclose).
        reduced_dummy = 0
        CALL MPI_Reduce( reduced_dummy, reduced, 2, MPI_INTEGER, MPI_SUM, mytask_io_group, mpi_comm_io_groups(1), ierr )
        obufsize = reduced(1)



        IF ( obufsize .LT. 0 ) THEN
          IF ( obufsize .EQ. -100 ) THEN         
            CALL ext_ncd_ioexit( Status )
            CALL ext_int_ioexit( Status )
            CALL ext_gr1_ioexit( Status )
            CALL wrf_message ( 'I/O QUILT SERVERS DONE' )
            CALL mpi_finalize(ierr)
            STOP
          ELSE
            WRITE(mess,*)'Possible 32-bit overflow on output server. Try larger nio_tasks_per_group in namelist.'
            CALL wrf_error_fatal3("<stdin>",1498,&
mess)
          ENDIF
        ENDIF







        IF ( obufsize .GT. 0 ) THEN
          ALLOCATE( obuf( (obufsize+1)/itypesize ) )


          CALL collect_on_comm_debug("module_io_quilt_old.F",1724, mpi_comm_io_groups(1),        &
                                onebyte,                      &
                                dummy, 0,                     &
                                obuf, obufsize )

        ELSE
          ! Necessarily, the compute processes send the ioclose signal,
          
          ! will stall on the ioclose message waiting for the quilt 
          
          
          
          
          
          ! Then a header representing the ioclose message is constructed
          
          
          
          
          
          
          ALLOCATE( obuf( 4096 ) )
          
          CALL int_gen_handle_header( obuf, obufsize, itypesize, &
                                      reduced(2) , int_ioclose )
        ENDIF












        CALL init_store_piece_of_field
        CALL mpi_type_size ( MPI_INTEGER , itypesize , ierr )



        vid = 0
        icurs = itypesize
        num_noops = 0 
        num_commit_messages = 0 
        num_field_training_msgs = 0 
        DO WHILE ( icurs .lt. obufsize ) 
          hdr_tag = get_hdr_tag( obuf ( icurs / itypesize ) )
          SELECT CASE ( hdr_tag )
            CASE ( int_field )
              CALL int_get_write_field_header ( obuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                DomainStart , DomainEnd ,                                    &
                                                MemoryStart , MemoryEnd ,                                    &
                                                PatchStart , PatchEnd )
              chunksize = (PatchEnd(1)-PatchStart(1)+1)*(PatchEnd(2)-PatchStart(2)+1)* &
                          (PatchEnd(3)-PatchStart(3)+1)*ftypesize

              IF ( DomainDesc .EQ. 333933 ) THEN  
                 IF ( num_field_training_msgs .EQ. 0 ) THEN
                   call add_to_bufsize_for_field( VarName, hdrbufsize )

                 ENDIF
                 num_field_training_msgs = num_field_training_msgs + 1
              ELSE
                 call add_to_bufsize_for_field( VarName, hdrbufsize )

              ENDIF
              icurs = icurs + hdrbufsize



              
              
              IF ( DomainDesc .NE. 333933 ) THEN   

                call add_to_bufsize_for_field( VarName, chunksize )
                icurs = icurs + chunksize
              ENDIF
            CASE ( int_open_for_write_commit )  
              hdrbufsize = obuf(icurs/itypesize)
              IF (num_commit_messages.EQ.0) THEN
                call add_to_bufsize_for_field( 'COMMIT', hdrbufsize )
              ENDIF
              num_commit_messages = num_commit_messages + 1
              icurs = icurs + hdrbufsize
            CASE DEFAULT
              hdrbufsize = obuf(icurs/itypesize)





























              IF (hdr_tag.NE.int_noop) THEN
                write(VarName,'(I5.5)')vid 

                call add_to_bufsize_for_field( VarName, hdrbufsize )
                vid = vid+1
              ENDIF
              IF ( hdr_tag .EQ. int_noop ) num_noops = num_noops + 1
              icurs = icurs + hdrbufsize

          END SELECT
        ENDDO 



        vid = 0
        icurs = itypesize
        num_noops = 0 
        num_commit_messages = 0 
        num_field_training_msgs = 0 
        DO WHILE ( icurs .lt. obufsize ) 

          hdr_tag = get_hdr_tag( obuf ( icurs / itypesize ) )
          SELECT CASE ( hdr_tag )
            CASE ( int_field )
              CALL int_get_write_field_header ( obuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                DomainStart , DomainEnd ,                                    &
                                                MemoryStart , MemoryEnd ,                                    &
                                                PatchStart , PatchEnd )
              chunksize = (PatchEnd(1)-PatchStart(1)+1)*(PatchEnd(2)-PatchStart(2)+1)* &
                          (PatchEnd(3)-PatchStart(3)+1)*ftypesize

              IF ( DomainDesc .EQ. 333933 ) THEN  
                 IF ( num_field_training_msgs .EQ. 0 ) THEN
                   call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )

                 ENDIF
                 num_field_training_msgs = num_field_training_msgs + 1
              ELSE
                 call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )

              ENDIF
              icurs = icurs + hdrbufsize
              
              
              IF ( DomainDesc .NE. 333933 ) THEN   
                call store_piece_of_field( obuf(icurs/itypesize), VarName, chunksize )
                icurs = icurs + chunksize

              ENDIF
            CASE ( int_open_for_write_commit )  
              hdrbufsize = obuf(icurs/itypesize)
              IF (num_commit_messages.EQ.0) THEN
                call store_piece_of_field( obuf(icurs/itypesize), 'COMMIT', hdrbufsize )
              ENDIF
              num_commit_messages = num_commit_messages + 1
              icurs = icurs + hdrbufsize
            CASE DEFAULT
              hdrbufsize = obuf(icurs/itypesize)
              IF (hdr_tag.NE.int_noop) THEN

                write(VarName,'(I5.5)')vid 

                call store_piece_of_field( obuf(icurs/itypesize), VarName, hdrbufsize )
                vid = vid+1
              ENDIF
              IF ( hdr_tag .EQ. int_noop ) num_noops = num_noops + 1
              icurs = icurs + hdrbufsize
          END SELECT
       ENDDO 



       CALL init_retrieve_pieces_of_field


       CALL retrieve_pieces_of_field ( obuf , VarName, obufsize, sz, retval )
       written_record = .false.


       DO WHILE ( retval ) 




            icurs = itypesize  

            stored_write_record = .false.



            DO WHILE ( icurs .lt. sz)





              SELECT CASE ( get_hdr_tag( obuf(icurs/itypesize) ) )


                CASE ( int_noop )
                  CALL int_get_noop_header( obuf(icurs/itypesize), &
                                            hdrbufsize, itypesize )
                  icurs = icurs + hdrbufsize


                CASE ( int_dom_td_real )
                  CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                  ALLOCATE( RData( obuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_td_header( obuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, DateStr, Element, RData, Count, code )
                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_td_real( handle(DataHandle),TRIM(Element),TRIM(DateStr),RData, Count, Status )
                     CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( RData )

                CASE ( int_dom_ti_real )

                  CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                  ALLOCATE( RData( obuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_ti_header( obuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, Element, RData, Count, code )
                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_ti_real( handle(DataHandle),TRIM(Element), RData, Count, Status )
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( RData )


                CASE ( int_dom_td_integer )

                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  ALLOCATE( IData( obuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_td_header( obuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, DateStr, Element, IData, Count, code )
                  icurs = icurs + hdrbufsize

                  SELECT CASE (use_package(io_form(DataHandle)))
                   CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                   CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                   CASE ( IO_GRIB1 )
                      CALL ext_gr1_put_dom_td_integer( handle(DataHandle),TRIM(Element), Trim(DateStr), IData, Count, Status )
                   CASE DEFAULT
                      Status = 0
                   END SELECT

                   DEALLOCATE( IData )


                CASE ( int_dom_ti_integer )

                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  ALLOCATE( IData( obuf(icurs/itypesize + 4 ) ) )      
                  CALL int_get_ti_header( obuf(icurs/itypesize:), hdrbufsize, itypesize, ftypesize, &
                                          DataHandle, Element, IData, Count, code )
                  icurs = icurs + hdrbufsize
                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_put_dom_ti_integer( handle(DataHandle),TRIM(Element), IData, Count, Status )

                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  DEALLOCATE( IData)
 

                CASE ( int_set_time )

                  CALL int_get_ti_header_char( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )
                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_INTIO   )
                      CALL ext_int_set_time ( handle(DataHandle), TRIM(CData), Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT

                  icurs = icurs + hdrbufsize


                CASE ( int_dom_ti_char )

                  CALL int_get_ti_header_char( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                   CASE ( IO_GRIB1 )
                      CALL ext_gr1_put_dom_ti_char ( handle(DataHandle), TRIM(Element), TRIM(CData), Status)
                   CASE DEFAULT
                      Status = 0
                   END SELECT

                  icurs = icurs + hdrbufsize


                CASE ( int_var_ti_char )

                  CALL int_get_ti_header_char( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                               DataHandle, Element, VarName, CData, code )

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                    CASE ( IO_INTIO   )
                      CALL ext_int_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                   CASE ( IO_GRIB1 )
                      CALL ext_gr1_put_var_ti_char ( handle(DataHandle), TRIM(Element), TRIM(VarName), TRIM(CData), Status)
                   CASE DEFAULT
                      Status = 0
                   END SELECT

                  icurs = icurs + hdrbufsize

                CASE ( int_ioexit )

                  CALL wrf_error_fatal3("<stdin>",1879,&
                         "quilt: should have handled int_ioexit already")
! Every I/O server handles the "ioclose" request.
                CASE ( int_ioclose )
                  CALL int_get_handle_header( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                              DataHandle , code )
                  icurs = icurs + hdrbufsize

                  IF ( DataHandle .GE. 1 ) THEN

                     SELECT CASE (use_package(io_form(DataHandle)))
                     CASE ( IO_NETCDF   )
                        CALL ext_ncd_inquire_filename( handle(DataHandle), fname, fstat, Status )
                        IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                           CALL ext_ncd_ioclose(handle(DataHandle),Status)
                        ENDIF
                     CASE ( IO_INTIO   )
                        CALL ext_int_inquire_filename( handle(DataHandle), fname, fstat, Status )
                        IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                           CALL ext_int_ioclose(handle(DataHandle),Status)
                        ENDIF
                     CASE ( IO_GRIB1 )
                        CALL ext_gr1_inquire_filename( handle(DataHandle), fname, fstat, Status )
                        IF ( fstat .EQ. WRF_FILE_OPENED_FOR_WRITE .OR. fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                           CALL ext_gr1_ioclose(handle(DataHandle),Status)
                        ENDIF
                     CASE DEFAULT
                        Status = 0
                     END SELECT
                  ENDIF


                CASE ( int_open_for_write_begin )

                  CALL int_get_ofwb_header( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                            FileName,SysDepInfo,io_form_arg,DataHandle )





                  icurs = icurs + hdrbufsize

                
                  io_form(DataHandle) = io_form_arg

                  SELECT CASE (use_package(io_form(DataHandle)))
                    CASE ( IO_NETCDF   )
                      CALL ext_ncd_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)

                    CASE ( IO_INTIO   )
                      CALL ext_int_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)
                    CASE ( IO_GRIB1 )
                       CALL ext_gr1_open_for_write_begin(FileName,Comm,IOComm,SysDepInfo,handle(DataHandle),Status)
                    CASE DEFAULT
                      Status = 0
                  END SELECT
                
                  okay_to_write(DataHandle) = .false.





                CASE ( int_open_for_write_commit )

                  CALL int_get_handle_header( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                              DataHandle , code )
                  icurs = icurs + hdrbufsize
                  okay_to_commit(DataHandle) = .true.










                CASE ( int_field )
                  CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                  CALL int_get_write_field_header ( obuf(icurs/itypesize), hdrbufsize, itypesize, ftypesize,  &
                                                    DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                                    DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                                    DomainStart , DomainEnd ,                                    &
                                                    MemoryStart , MemoryEnd ,                                    &
                                                    PatchStart , PatchEnd )

                  icurs = icurs + hdrbufsize

                  IF ( okay_to_write(DataHandle) ) THEN









                    IF ( FieldType .EQ. WRF_FLOAT .OR. FieldType .EQ. WRF_DOUBLE)  THEN
                      
                      
                      IF ( FieldType .EQ. WRF_DOUBLE)  THEN

                        CALL mpi_type_size( MPI_DOUBLE_PRECISION, ftypesize, ierr )
                      ELSE
                        CALL mpi_type_size( MPI_REAL, ftypesize, ierr )
                      ENDIF

                    ELSE IF ( FieldType .EQ. WRF_INTEGER ) THEN
                      CALL mpi_type_size( MPI_INTEGER, ftypesize, ierr )
                    ELSE IF ( FieldType .EQ. WRF_LOGICAL ) THEN
                      ftypesize = 4
                    ENDIF

                    icurs = icurs + (PatchEnd(1)-PatchStart(1)+1)* &
                                    (PatchEnd(2)-PatchStart(2)+1)* &
                                    (PatchEnd(3)-PatchStart(3)+1)*ftypesize

                  ELSE 

                    SELECT CASE (use_package(io_form(DataHandle)))

                      CASE ( IO_NETCDF   )
                        CALL ext_ncd_write_field ( handle(DataHandle) , TRIM(DateStr) ,         &
                                   TRIM(VarName) , dummy , FieldType , Comm , IOComm,           &
                                   DomainDesc , TRIM(MemoryOrder) , TRIM(Stagger) , DimNames ,  &
                                   DomainStart , DomainEnd ,                                    &
                                   DomainStart , DomainEnd ,                                    &
                                   DomainStart , DomainEnd ,                                    &
                                   Status )
                      CASE DEFAULT
                        Status = 0
                    END SELECT
                  ENDIF
                CASE ( int_iosync )
                  CALL int_get_handle_header( obuf(icurs/itypesize), hdrbufsize, itypesize, &
                                            DataHandle , code )
                  icurs = icurs + hdrbufsize
                CASE DEFAULT
                  WRITE(mess,*)'quilt: bad tag: ',                            &
                               get_hdr_tag( obuf(icurs/itypesize) ),' icurs ',&
                               icurs/itypesize
                  CALL wrf_error_fatal3("<stdin>",2024,&
mess )
              END SELECT

            ENDDO 



            IF (stored_write_record) THEN








              stored_write_record = .false.
              written_record = .true.
            ENDIF




            IF (okay_to_commit(DataHandle)) THEN

              SELECT CASE (use_package(io_form(DataHandle)))
                CASE ( IO_NETCDF   )
                  CALL ext_ncd_inquire_filename( handle(DataHandle), fname, fstat, Status )
                  IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                    CALL ext_ncd_open_for_write_commit(handle(DataHandle),Status)
                    okay_to_write(DataHandle) = .true.
                  ENDIF
                CASE ( IO_INTIO   )
                  CALL ext_int_inquire_filename( handle(DataHandle), fname, fstat, Status )
                  IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                    CALL ext_int_open_for_write_commit(handle(DataHandle),Status)
                    okay_to_write(DataHandle) = .true.
                  ENDIF
                 CASE ( IO_GRIB1 )
                    CALL ext_gr1_inquire_filename( handle(DataHandle), fname, fstat, Status )
                    IF ( fstat .EQ. WRF_FILE_OPENED_NOT_COMMITTED ) THEN
                       CALL ext_gr1_open_for_write_commit(handle(DataHandle),Status)
                       okay_to_write(DataHandle) = .true.
                    ENDIF

                CASE DEFAULT
                  Status = 0
              END SELECT

            okay_to_commit(DataHandle) = .false.
          ENDIF




        CALL retrieve_pieces_of_field ( obuf , VarName, obufsize, sz, retval )
      END DO 

      DEALLOCATE( obuf )

      
      IF (written_record) THEN

        SELECT CASE ( use_package(io_form) )
          CASE DEFAULT
            Status = 0
        END SELECT
        written_record = .false.

      ENDIF

      END DO 

    END SUBROUTINE quilt_pnc



    SUBROUTINE init_module_wrf_quilt
      USE module_wrf_error, only: init_module_wrf_error
      USE module_driver_constants
      USE module_dm, only: mpi_comm_allcompute








      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER i
      NAMELIST /namelist_quilt/ nio_tasks_per_group, nio_groups, poll_servers
      INTEGER ntasks, mytask, ierr, io_status
      INTEGER mpi_comm_here, temp_poll
      LOGICAL mpi_inited
      LOGICAL esmf_coupling









integer    :: first_item_in_struct
integer :: run_days
integer :: run_hours
integer :: run_minutes
integer :: run_seconds
integer , DIMENSION(max_domains) :: start_year
integer , DIMENSION(max_domains) :: start_month
integer , DIMENSION(max_domains) :: start_day
integer , DIMENSION(max_domains) :: start_hour
integer , DIMENSION(max_domains) :: start_minute
integer , DIMENSION(max_domains) :: start_second
integer , DIMENSION(max_domains) :: end_year
integer , DIMENSION(max_domains) :: end_month
integer , DIMENSION(max_domains) :: end_day
integer , DIMENSION(max_domains) :: end_hour
integer , DIMENSION(max_domains) :: end_minute
integer , DIMENSION(max_domains) :: end_second
integer :: interval_seconds
logical , DIMENSION(max_domains) :: input_from_file
integer , DIMENSION(max_domains) :: fine_input_stream
logical , DIMENSION(max_domains) :: input_from_hires
character*256 :: rsmas_data_path
logical :: all_ic_times
integer , DIMENSION(max_domains) :: julyr
integer , DIMENSION(max_domains) :: julday
real , DIMENSION(max_domains) :: gmt
character*256 :: input_inname
character*256 :: input_outname
character*256 :: bdy_inname
character*256 :: bdy_outname
character*256 :: rst_inname
character*256 :: rst_outname
logical :: write_input
logical :: write_restart_at_0h
logical :: write_hist_at_0h_rst
logical :: adjust_output_times
logical :: adjust_input_times
integer :: diag_print
logical :: nocolons
logical :: cycling
integer :: output_diagnostics
integer :: nwp_diagnostics
logical :: output_ready_flag
logical :: usepio
integer :: pioprocs
integer :: piostart
integer :: piostride
integer :: pioshift
integer :: dfi_opt
integer :: dfi_savehydmeteors
integer :: dfi_nfilter
logical :: dfi_write_filtered_input
logical :: dfi_write_dfi_history
integer :: dfi_cutoff_seconds
integer :: dfi_time_dim
integer :: dfi_fwdstop_year
integer :: dfi_fwdstop_month
integer :: dfi_fwdstop_day
integer :: dfi_fwdstop_hour
integer :: dfi_fwdstop_minute
integer :: dfi_fwdstop_second
integer :: dfi_bckstop_year
integer :: dfi_bckstop_month
integer :: dfi_bckstop_day
integer :: dfi_bckstop_hour
integer :: dfi_bckstop_minute
integer :: dfi_bckstop_second
integer :: time_step
integer :: time_step_fract_num
integer :: time_step_fract_den
integer :: time_step_dfi
integer , DIMENSION(max_domains) :: min_time_step
integer , DIMENSION(max_domains) :: min_time_step_den
integer , DIMENSION(max_domains) :: max_time_step
integer , DIMENSION(max_domains) :: max_time_step_den
real , DIMENSION(max_domains) :: target_cfl
real , DIMENSION(max_domains) :: target_hcfl
integer , DIMENSION(max_domains) :: max_step_increase_pct
integer , DIMENSION(max_domains) :: starting_time_step
integer , DIMENSION(max_domains) :: starting_time_step_den
logical :: step_to_output_time
integer :: adaptation_domain
logical :: use_adaptive_time_step
logical :: use_adaptive_time_step_dfi
integer :: max_dom
integer :: lats_to_mic
integer , DIMENSION(max_domains) :: s_we
integer , DIMENSION(max_domains) :: e_we
integer , DIMENSION(max_domains) :: s_sn
integer , DIMENSION(max_domains) :: e_sn
integer , DIMENSION(max_domains) :: s_vert
integer , DIMENSION(max_domains) :: e_vert
integer :: num_metgrid_levels
integer :: num_metgrid_soil_levels
real :: p_top_requested
logical :: interp_theta
integer :: interp_type
integer :: rebalance
integer , DIMENSION(max_domains) :: vert_refine_method
integer :: vert_refine_fact
integer :: extrap_type
integer :: t_extrap_type
integer :: hypsometric_opt
logical :: lowest_lev_from_sfc
logical :: use_levels_below_ground
logical :: use_tavg_for_tsk
logical :: use_surface
integer :: lagrange_order
integer :: force_sfc_in_vinterp
real :: zap_close_levels
real :: maxw_horiz_pres_diff
real :: trop_horiz_pres_diff
real :: maxw_above_this_level
integer :: use_maxw_level
integer :: use_trop_level
logical :: sfcp_to_sfcp
logical :: adjust_heights
logical :: smooth_cg_topo
integer :: nest_interp_coord
integer :: interp_method_type
logical :: aggregate_lu
logical :: rh2qv_wrt_liquid
integer :: rh2qv_method
real :: qv_max_p_safe
real :: qv_max_flag
real :: qv_max_value
real :: qv_min_p_safe
real :: qv_min_flag
real :: qv_min_value
integer :: ideal_init_method
real , DIMENSION(max_domains) :: dx
real , DIMENSION(max_domains) :: dy
integer , DIMENSION(max_domains) :: grid_id
logical , DIMENSION(max_domains) :: grid_allowed
integer , DIMENSION(max_domains) :: parent_id
integer , DIMENSION(max_domains) :: i_parent_start
integer , DIMENSION(max_domains) :: j_parent_start
integer , DIMENSION(max_domains) :: parent_grid_ratio
integer , DIMENSION(max_domains) :: parent_time_step_ratio
integer :: feedback
integer :: smooth_option
integer :: blend_width
real , DIMENSION(max_domains) :: ztop
integer , DIMENSION(max_domains) :: moad_grid_ratio
integer , DIMENSION(max_domains) :: moad_time_step_ratio
integer , DIMENSION(max_domains) :: shw
integer :: tile_sz_x
integer :: tile_sz_y
integer :: numtiles
integer :: numtiles_inc
integer :: numtiles_x
integer :: numtiles_y
integer :: tile_strategy
integer :: nproc_x
integer :: nproc_y
integer :: irand
real , DIMENSION(max_domains) :: dt
integer :: fft_used
integer :: cu_used
integer :: shcu_used
integer :: cam_used
integer :: alloc_qndropsource
integer :: num_moves
integer :: ts_buf_size
integer :: max_ts_locs
integer , DIMENSION(max_domains) :: vortex_interval
integer , DIMENSION(max_domains) :: max_vortex_speed
integer , DIMENSION(max_domains) :: corral_dist
integer :: track_level
real , DIMENSION(max_domains) :: time_to_move
integer , DIMENSION(max_moves) :: move_id
integer , DIMENSION(max_moves) :: move_interval
integer , DIMENSION(max_moves) :: move_cd_x
integer , DIMENSION(max_moves) :: move_cd_y
logical , DIMENSION(max_domains) :: swap_x
logical , DIMENSION(max_domains) :: swap_y
logical , DIMENSION(max_domains) :: cycle_x
logical , DIMENSION(max_domains) :: cycle_y
logical :: reorder_mesh
logical :: perturb_input
real , DIMENSION(max_eta) :: eta_levels
real :: max_dz
integer :: ocean_levels
real , DIMENSION(max_ocean) :: ocean_z
real , DIMENSION(max_ocean) :: ocean_t
real , DIMENSION(max_ocean) :: ocean_s
integer :: num_traj
integer :: max_ts_level
integer :: track_loc_in
integer :: num_ext_model_couple_dom
logical :: insert_bogus_storm
logical :: remove_storm
integer :: num_storm
real , DIMENSION(max_bogus) :: latc_loc
real , DIMENSION(max_bogus) :: lonc_loc
real , DIMENSION(max_bogus) :: vmax_meters_per_second
real , DIMENSION(max_bogus) :: rmax
real , DIMENSION(max_bogus) :: vmax_ratio
real :: rankine_lid
logical :: force_read_thompson
logical :: write_thompson_tables
integer , DIMENSION(max_domains) :: mp_physics
real , DIMENSION(max_domains) :: nssl_cccn
real , DIMENSION(max_domains) :: nssl_alphah
real , DIMENSION(max_domains) :: nssl_alphahl
real , DIMENSION(max_domains) :: nssl_cnoh
real , DIMENSION(max_domains) :: nssl_cnohl
real , DIMENSION(max_domains) :: nssl_cnor
real , DIMENSION(max_domains) :: nssl_cnos
real , DIMENSION(max_domains) :: nssl_rho_qh
real , DIMENSION(max_domains) :: nssl_rho_qhl
real , DIMENSION(max_domains) :: nssl_rho_qs
integer , DIMENSION(max_domains) :: nudge_lightning
integer , DIMENSION(max_domains) :: nudge_light_times
integer , DIMENSION(max_domains) :: nudge_light_timee
integer , DIMENSION(max_domains) :: nudge_light_int
character*256 :: path_to_files
integer :: gsfcgce_hail
integer :: gsfcgce_2ice
integer , DIMENSION(max_domains) :: progn
real :: accum_mode
real :: aitken_mode
real :: coarse_mode
integer :: do_radar_ref
integer :: compute_radar_ref
integer , DIMENSION(max_domains) :: ra_lw_physics
integer , DIMENSION(max_domains) :: ra_sw_physics
real , DIMENSION(max_domains) :: radt
real , DIMENSION(max_domains) :: naer
integer , DIMENSION(max_domains) :: sf_sfclay_physics
integer , DIMENSION(max_domains) :: sf_surface_physics
integer , DIMENSION(max_domains) :: bl_pbl_physics
integer , DIMENSION(max_domains) :: bl_mynn_tkebudget
integer :: ysu_topdown_pblmix
integer , DIMENSION(max_domains) :: shinhong_tke_diag
logical , DIMENSION(max_domains) :: bl_mynn_tkeadvect
integer :: bl_mynn_cloudpdf
integer :: bl_mynn_mixlength
integer , DIMENSION(max_domains) :: bl_mynn_edmf
integer , DIMENSION(max_domains) :: bl_mynn_edmf_mom
integer , DIMENSION(max_domains) :: bl_mynn_edmf_tke
integer , DIMENSION(max_domains) :: bl_mynn_edmf_part
integer , DIMENSION(max_domains) :: bl_mynn_cloudmix
integer , DIMENSION(max_domains) :: bl_mynn_mixqt
integer :: icloud_bl
integer , DIMENSION(max_domains) :: mfshconv
integer , DIMENSION(max_domains) :: sf_urban_physics
real , DIMENSION(max_domains) :: bldt
integer , DIMENSION(max_domains) :: cu_physics
integer , DIMENSION(max_domains) :: shcu_physics
integer , DIMENSION(max_domains) :: cu_diag
integer , DIMENSION(max_domains) :: kf_edrates
integer :: kfeta_trigger
integer :: nsas_dx_factor
real , DIMENSION(max_domains) :: cudt
real , DIMENSION(max_domains) :: gsmdt
integer :: isfflx
integer :: ifsnow
integer :: icloud
integer :: ideal_xland
real :: swrad_scat
integer :: surface_input_source
integer :: num_soil_layers
integer :: maxpatch
integer :: num_snow_layers
integer :: num_snso_layers
integer :: num_urban_layers
integer :: num_urban_hi
integer :: num_months
integer :: sf_surface_mosaic
integer :: mosaic_cat
integer :: mosaic_cat_soil
integer :: mosaic_lu
integer :: mosaic_soil
integer :: maxiens
integer :: maxens
integer :: maxens2
integer :: maxens3
integer :: ensdim
integer :: cugd_avedx
integer :: clos_choice
integer :: imomentum
integer :: ishallow
real :: convtrans_avglen_m
integer :: num_land_cat
integer :: num_soil_cat
integer :: mp_zero_out
real :: mp_zero_out_thresh
real :: seaice_threshold
integer :: sst_update
integer :: sst_skin
integer :: tmn_update
logical :: usemonalb
logical :: rdmaxalb
logical :: rdlai2d
logical :: ua_phys
integer :: opt_thcnd
integer :: co2tf
integer :: ra_call_offset
real :: cam_abs_freq_s
integer :: levsiz
integer :: paerlev
integer :: cam_abs_dim1
integer :: cam_abs_dim2
integer :: lagday
integer :: no_src_types
integer :: alevsiz
integer :: o3input
integer :: aer_opt
integer :: swint_opt
integer , DIMENSION(max_domains) :: aer_type
integer , DIMENSION(max_domains) :: aer_aod550_opt
integer , DIMENSION(max_domains) :: aer_angexp_opt
integer , DIMENSION(max_domains) :: aer_ssa_opt
integer , DIMENSION(max_domains) :: aer_asy_opt
real , DIMENSION(max_domains) :: aer_aod550_val
real , DIMENSION(max_domains) :: aer_angexp_val
real , DIMENSION(max_domains) :: aer_ssa_val
real , DIMENSION(max_domains) :: aer_asy_val
logical , DIMENSION(max_domains) :: cu_rad_feedback
logical , DIMENSION(max_domains) :: shallowcu_forced_ra
integer , DIMENSION(max_domains) :: numbins
real , DIMENSION(max_domains) :: thbinsize
real , DIMENSION(max_domains) :: rbinsize
real , DIMENSION(max_domains) :: mindeepfreq
real , DIMENSION(max_domains) :: minshallowfreq
integer , DIMENSION(max_domains) :: shcu_aerosols_opt
integer , DIMENSION(max_domains) :: icloud_cu
integer , DIMENSION(max_domains) :: pxlsm_smois_init
integer :: omlcall
integer :: sf_ocean_physics
integer :: traj_opt
integer :: tracercall
real :: omdt
real :: oml_hml0
real :: oml_gamma
real :: oml_relaxation_time
integer :: isftcflx
integer :: iz0tlnd
real :: shadlen
integer , DIMENSION(max_domains) :: slope_rad
integer , DIMENSION(max_domains) :: topo_shading
integer , DIMENSION(max_domains) :: topo_wind
integer :: no_mp_heating
integer :: fractional_seaice
integer :: seaice_snowdepth_opt
real :: seaice_snowdepth_max
real :: seaice_snowdepth_min
integer :: seaice_albedo_opt
real :: seaice_albedo_default
integer :: seaice_thickness_opt
real :: seaice_thickness_default
logical :: tice2tsk_if2cold
real :: bucket_mm
real :: bucket_j
real :: mp_tend_lim
real , DIMENSION(max_domains) :: prec_acc_dt
integer :: prec_acc_opt
integer :: bucketr_opt
integer :: process_time_series
integer , DIMENSION(max_domains) :: grav_settling
real , DIMENSION(max_domains) :: sas_pgcon
integer , DIMENSION(max_domains) :: scalar_pblmix
integer , DIMENSION(max_domains) :: tracer_pblmix
logical :: use_aero_icbc
logical :: use_rap_aero_icbc
integer :: use_mp_re
real :: ccn_conc
integer :: hail_opt
integer :: dveg
integer :: opt_crs
integer :: opt_btr
integer :: opt_run
integer :: opt_sfc
integer :: opt_frz
integer :: opt_inf
integer :: opt_rad
integer :: opt_alb
integer :: opt_snf
integer :: opt_tbot
integer :: opt_stc
integer :: opt_gla
integer :: opt_rsf
real , DIMENSION(max_domains) :: wtddt
integer :: wrf_hydro
real , DIMENSION(max_domains) :: fgdt
integer , DIMENSION(max_domains) :: fgdtzero
integer , DIMENSION(max_domains) :: grid_fdda
integer , DIMENSION(max_domains) :: grid_sfdda
integer , DIMENSION(max_domains) :: if_no_pbl_nudging_uv
integer , DIMENSION(max_domains) :: if_no_pbl_nudging_t
integer , DIMENSION(max_domains) :: if_no_pbl_nudging_ph
integer , DIMENSION(max_domains) :: if_no_pbl_nudging_q
integer , DIMENSION(max_domains) :: if_zfac_uv
integer , DIMENSION(max_domains) :: k_zfac_uv
integer , DIMENSION(max_domains) :: if_zfac_t
integer , DIMENSION(max_domains) :: k_zfac_t
integer , DIMENSION(max_domains) :: if_zfac_ph
integer , DIMENSION(max_domains) :: k_zfac_ph
integer , DIMENSION(max_domains) :: if_zfac_q
integer , DIMENSION(max_domains) :: k_zfac_q
integer , DIMENSION(max_domains) :: dk_zfac_uv
integer , DIMENSION(max_domains) :: dk_zfac_t
integer , DIMENSION(max_domains) :: dk_zfac_ph
real , DIMENSION(max_domains) :: guv
real , DIMENSION(max_domains) :: guv_sfc
real , DIMENSION(max_domains) :: gt
real , DIMENSION(max_domains) :: gt_sfc
real , DIMENSION(max_domains) :: gq
real , DIMENSION(max_domains) :: gq_sfc
real , DIMENSION(max_domains) :: gph
real :: dtramp_min
integer :: if_ramping
real , DIMENSION(max_domains) :: rinblw
integer , DIMENSION(max_domains) :: xwavenum
integer , DIMENSION(max_domains) :: ywavenum
integer , DIMENSION(max_domains) :: pxlsm_soil_nudge
integer , DIMENSION(max_domains) :: fasdas
integer , DIMENSION(max_domains) :: obs_nudge_opt
integer :: max_obs
real , DIMENSION(max_domains) :: fdda_start
real , DIMENSION(max_domains) :: fdda_end
integer , DIMENSION(max_domains) :: obs_nudge_wind
real , DIMENSION(max_domains) :: obs_coef_wind
integer , DIMENSION(max_domains) :: obs_nudge_temp
real , DIMENSION(max_domains) :: obs_coef_temp
integer , DIMENSION(max_domains) :: obs_nudge_mois
real , DIMENSION(max_domains) :: obs_coef_mois
integer , DIMENSION(max_domains) :: obs_nudge_pstr
real , DIMENSION(max_domains) :: obs_coef_pstr
integer , DIMENSION(max_domains) :: obs_no_pbl_nudge_uv
integer , DIMENSION(max_domains) :: obs_no_pbl_nudge_t
integer , DIMENSION(max_domains) :: obs_no_pbl_nudge_q
integer :: obs_sfc_scheme_horiz
integer :: obs_sfc_scheme_vert
real :: obs_max_sndng_gap
real :: obs_nudgezfullr1_uv
real :: obs_nudgezrampr1_uv
real :: obs_nudgezfullr2_uv
real :: obs_nudgezrampr2_uv
real :: obs_nudgezfullr4_uv
real :: obs_nudgezrampr4_uv
real :: obs_nudgezfullr1_t
real :: obs_nudgezrampr1_t
real :: obs_nudgezfullr2_t
real :: obs_nudgezrampr2_t
real :: obs_nudgezfullr4_t
real :: obs_nudgezrampr4_t
real :: obs_nudgezfullr1_q
real :: obs_nudgezrampr1_q
real :: obs_nudgezfullr2_q
real :: obs_nudgezrampr2_q
real :: obs_nudgezfullr4_q
real :: obs_nudgezrampr4_q
real :: obs_nudgezfullmin
real :: obs_nudgezrampmin
real :: obs_nudgezmax
real :: obs_sfcfact
real :: obs_sfcfacr
real :: obs_dpsmx
real , DIMENSION(max_domains) :: obs_rinxy
real :: obs_rinsig
real , DIMENSION(max_domains) :: obs_twindo
integer :: obs_npfi
integer , DIMENSION(max_domains) :: obs_ionf
integer :: obs_idynin
real :: obs_dtramp
integer :: obs_prt_max
integer , DIMENSION(max_domains) :: obs_prt_freq
logical :: obs_ipf_in4dob
logical :: obs_ipf_errob
logical :: obs_ipf_nudob
logical :: obs_ipf_init
integer :: obs_scl_neg_qv_innov
integer :: scm_force
real :: scm_force_dx
integer :: num_force_layers
integer :: scm_lu_index
integer :: scm_isltyp
real :: scm_vegfra
real :: scm_canwat
real :: scm_lat
real :: scm_lon
logical :: scm_th_t_tend
logical :: scm_qv_t_tend
logical :: scm_th_adv
logical :: scm_wind_adv
logical :: scm_qv_adv
logical :: scm_ql_adv
logical :: scm_vert_adv
integer :: num_force_soil_layers
logical :: scm_soilt_force
logical :: scm_soilq_force
logical :: scm_force_th_largescale
logical :: scm_force_qv_largescale
logical :: scm_force_ql_largescale
logical :: scm_force_wind_largescale
integer :: scm_force_skintemp
integer :: scm_force_flux
integer :: dyn_opt
integer :: rk_ord
integer :: w_damping
integer , DIMENSION(max_domains) :: diff_opt
integer , DIMENSION(max_domains) :: diff_opt_dfi
integer , DIMENSION(max_domains) :: km_opt
integer , DIMENSION(max_domains) :: km_opt_dfi
integer :: damp_opt
integer :: rad_nudge
integer :: gwd_opt
real , DIMENSION(max_domains) :: zdamp
real , DIMENSION(max_domains) :: dampcoef
real , DIMENSION(max_domains) :: khdif
real , DIMENSION(max_domains) :: kvdif
real , DIMENSION(max_domains) :: diff_6th_factor
integer , DIMENSION(max_domains) :: diff_6th_opt
integer :: use_theta_m
integer :: use_q_diabatic
real , DIMENSION(max_domains) :: c_s
real , DIMENSION(max_domains) :: c_k
real , DIMENSION(max_domains) :: smdiv
real , DIMENSION(max_domains) :: emdiv
real , DIMENSION(max_domains) :: epssm
logical , DIMENSION(max_domains) :: non_hydrostatic
logical :: use_input_w
integer , DIMENSION(max_domains) :: time_step_sound
integer , DIMENSION(max_domains) :: h_mom_adv_order
integer , DIMENSION(max_domains) :: v_mom_adv_order
integer , DIMENSION(max_domains) :: h_sca_adv_order
integer , DIMENSION(max_domains) :: v_sca_adv_order
integer , DIMENSION(max_domains) :: momentum_adv_opt
integer , DIMENSION(max_domains) :: moist_adv_opt
integer , DIMENSION(max_domains) :: moist_adv_dfi_opt
integer , DIMENSION(max_domains) :: chem_adv_opt
integer , DIMENSION(max_domains) :: tracer_adv_opt
integer , DIMENSION(max_domains) :: scalar_adv_opt
integer , DIMENSION(max_domains) :: tke_adv_opt
logical , DIMENSION(max_domains) :: top_radiation
integer , DIMENSION(max_domains) :: mix_isotropic
real , DIMENSION(max_domains) :: mix_upper_bound
logical , DIMENSION(max_domains) :: top_lid
real , DIMENSION(max_domains) :: tke_upper_bound
real , DIMENSION(max_domains) :: tke_drag_coefficient
real , DIMENSION(max_domains) :: tke_heat_flux
logical , DIMENSION(max_domains) :: pert_coriolis
logical , DIMENSION(max_domains) :: coriolis2d
logical , DIMENSION(max_domains) :: mix_full_fields
real :: base_pres
real :: base_temp
real :: base_lapse
real :: iso_temp
real :: base_pres_strat
real :: base_lapse_strat
logical :: use_baseparam_fr_nml
real :: fft_filter_lat
logical :: coupled_filtering
logical :: pos_def
logical :: swap_pole_with_next_j
logical :: actual_distance_average
logical :: rotated_pole
logical , DIMENSION(max_domains) :: do_coriolis
logical , DIMENSION(max_domains) :: do_curvature
logical , DIMENSION(max_domains) :: do_gradp
integer , DIMENSION(max_domains) :: tracer_opt
integer , DIMENSION(max_domains) :: tenddiag
integer :: spec_bdy_width
integer :: spec_zone
integer :: relax_zone
logical , DIMENSION(max_domains) :: specified
logical :: constant_bc
logical , DIMENSION(max_domains) :: periodic_x
logical , DIMENSION(max_domains) :: symmetric_xs
logical , DIMENSION(max_domains) :: symmetric_xe
logical , DIMENSION(max_domains) :: open_xs
logical , DIMENSION(max_domains) :: open_xe
logical , DIMENSION(max_domains) :: periodic_y
logical , DIMENSION(max_domains) :: symmetric_ys
logical , DIMENSION(max_domains) :: symmetric_ye
logical , DIMENSION(max_domains) :: open_ys
logical , DIMENSION(max_domains) :: open_ye
logical , DIMENSION(max_domains) :: polar
logical , DIMENSION(max_domains) :: nested
real :: spec_exp
integer :: spec_bdy_final_mu
integer :: real_data_init_type
logical , DIMENSION(max_domains) :: have_bcs_moist
logical , DIMENSION(max_domains) :: have_bcs_scalar
integer :: background_proc_id
integer :: forecast_proc_id
integer :: production_status
integer :: compression
integer :: nobs_ndg_vars
integer :: nobs_err_flds
real , DIMENSION(max_domains) :: cen_lat
real , DIMENSION(max_domains) :: cen_lon
real , DIMENSION(max_domains) :: truelat1
real , DIMENSION(max_domains) :: truelat2
real , DIMENSION(max_domains) :: moad_cen_lat
real , DIMENSION(max_domains) :: stand_lon
real , DIMENSION(max_domains) :: pole_lat
real , DIMENSION(max_domains) :: pole_lon
integer :: flag_metgrid
integer :: flag_snow
integer :: flag_psfc
integer :: flag_sm000010
integer :: flag_sm010040
integer :: flag_sm040100
integer :: flag_sm100200
integer :: flag_st000010
integer :: flag_st010040
integer :: flag_st040100
integer :: flag_st100200
integer :: flag_soil_layers
integer :: flag_slp
integer :: flag_soilhgt
integer :: flag_mf_xy
integer :: flag_um_soil
real , DIMENSION(max_domains) :: bdyfrq
character*256 , DIMENSION(max_domains) :: mminlu
integer , DIMENSION(max_domains) :: iswater
integer , DIMENSION(max_domains) :: islake
integer , DIMENSION(max_domains) :: isice
integer , DIMENSION(max_domains) :: isurban
integer , DIMENSION(max_domains) :: isoilwater
integer , DIMENSION(max_domains) :: map_proj
integer :: use_wps_input
integer , DIMENSION(max_domains) :: dfi_stage
integer , DIMENSION(max_domains) :: mp_physics_dfi
integer , DIMENSION(max_domains) :: bl_pbl_physics_dfi
integer , DIMENSION(max_domains) :: windfarm_opt
integer :: windfarm_ij
integer , DIMENSION(max_domains) :: lightning_option
real , DIMENSION(max_domains) :: lightning_dt
real , DIMENSION(max_domains) :: lightning_start_seconds
real , DIMENSION(max_domains) :: flashrate_factor
integer , DIMENSION(max_domains) :: iccg_method
real , DIMENSION(max_domains) :: iccg_prescribed_num
real , DIMENSION(max_domains) :: iccg_prescribed_den
integer , DIMENSION(max_domains) :: cellcount_method
real , DIMENSION(max_domains) :: cldtop_adjustment
integer , DIMENSION(max_domains) :: sf_lake_physics
character*256 :: auxinput1_inname
integer :: io_form_auxinput1
logical :: override_restart_timers
character*256 :: auxhist1_inname
character*256 :: auxhist1_outname
integer , DIMENSION(max_domains) :: auxhist1_interval_y
integer , DIMENSION(max_domains) :: auxhist1_interval_d
integer , DIMENSION(max_domains) :: auxhist1_interval_h
integer , DIMENSION(max_domains) :: auxhist1_interval_m
integer , DIMENSION(max_domains) :: auxhist1_interval_s
integer , DIMENSION(max_domains) :: auxhist1_interval
integer , DIMENSION(max_domains) :: auxhist1_begin_y
integer , DIMENSION(max_domains) :: auxhist1_begin_d
integer , DIMENSION(max_domains) :: auxhist1_begin_h
integer , DIMENSION(max_domains) :: auxhist1_begin_m
integer , DIMENSION(max_domains) :: auxhist1_begin_s
integer , DIMENSION(max_domains) :: auxhist1_begin
integer , DIMENSION(max_domains) :: auxhist1_end_y
integer , DIMENSION(max_domains) :: auxhist1_end_d
integer , DIMENSION(max_domains) :: auxhist1_end_h
integer , DIMENSION(max_domains) :: auxhist1_end_m
integer , DIMENSION(max_domains) :: auxhist1_end_s
integer , DIMENSION(max_domains) :: auxhist1_end
integer :: io_form_auxhist1
integer , DIMENSION(max_domains) :: frames_per_auxhist1
character*256 :: auxhist2_inname
character*256 :: auxhist2_outname
integer , DIMENSION(max_domains) :: auxhist2_interval_y
integer , DIMENSION(max_domains) :: auxhist2_interval_d
integer , DIMENSION(max_domains) :: auxhist2_interval_h
integer , DIMENSION(max_domains) :: auxhist2_interval_m
integer , DIMENSION(max_domains) :: auxhist2_interval_s
integer , DIMENSION(max_domains) :: auxhist2_interval
integer , DIMENSION(max_domains) :: auxhist2_begin_y
integer , DIMENSION(max_domains) :: auxhist2_begin_d
integer , DIMENSION(max_domains) :: auxhist2_begin_h
integer , DIMENSION(max_domains) :: auxhist2_begin_m
integer , DIMENSION(max_domains) :: auxhist2_begin_s
integer , DIMENSION(max_domains) :: auxhist2_begin
integer , DIMENSION(max_domains) :: auxhist2_end_y
integer , DIMENSION(max_domains) :: auxhist2_end_d
integer , DIMENSION(max_domains) :: auxhist2_end_h
integer , DIMENSION(max_domains) :: auxhist2_end_m
integer , DIMENSION(max_domains) :: auxhist2_end_s
integer , DIMENSION(max_domains) :: auxhist2_end
integer :: io_form_auxhist2
integer , DIMENSION(max_domains) :: frames_per_auxhist2
character*256 :: auxhist3_inname
character*256 :: auxhist3_outname
integer , DIMENSION(max_domains) :: auxhist3_interval_y
integer , DIMENSION(max_domains) :: auxhist3_interval_d
integer , DIMENSION(max_domains) :: auxhist3_interval_h
integer , DIMENSION(max_domains) :: auxhist3_interval_m
integer , DIMENSION(max_domains) :: auxhist3_interval_s
integer , DIMENSION(max_domains) :: auxhist3_interval
integer , DIMENSION(max_domains) :: auxhist3_begin_y
integer , DIMENSION(max_domains) :: auxhist3_begin_d
integer , DIMENSION(max_domains) :: auxhist3_begin_h
integer , DIMENSION(max_domains) :: auxhist3_begin_m
integer , DIMENSION(max_domains) :: auxhist3_begin_s
integer , DIMENSION(max_domains) :: auxhist3_begin
integer , DIMENSION(max_domains) :: auxhist3_end_y
integer , DIMENSION(max_domains) :: auxhist3_end_d
integer , DIMENSION(max_domains) :: auxhist3_end_h
integer , DIMENSION(max_domains) :: auxhist3_end_m
integer , DIMENSION(max_domains) :: auxhist3_end_s
integer , DIMENSION(max_domains) :: auxhist3_end
integer :: io_form_auxhist3
integer , DIMENSION(max_domains) :: frames_per_auxhist3
character*256 :: auxhist4_inname
character*256 :: auxhist4_outname
integer , DIMENSION(max_domains) :: auxhist4_interval_y
integer , DIMENSION(max_domains) :: auxhist4_interval_d
integer , DIMENSION(max_domains) :: auxhist4_interval_h
integer , DIMENSION(max_domains) :: auxhist4_interval_m
integer , DIMENSION(max_domains) :: auxhist4_interval_s
integer , DIMENSION(max_domains) :: auxhist4_interval
integer , DIMENSION(max_domains) :: auxhist4_begin_y
integer , DIMENSION(max_domains) :: auxhist4_begin_d
integer , DIMENSION(max_domains) :: auxhist4_begin_h
integer , DIMENSION(max_domains) :: auxhist4_begin_m
integer , DIMENSION(max_domains) :: auxhist4_begin_s
integer , DIMENSION(max_domains) :: auxhist4_begin
integer , DIMENSION(max_domains) :: auxhist4_end_y
integer , DIMENSION(max_domains) :: auxhist4_end_d
integer , DIMENSION(max_domains) :: auxhist4_end_h
integer , DIMENSION(max_domains) :: auxhist4_end_m
integer , DIMENSION(max_domains) :: auxhist4_end_s
integer , DIMENSION(max_domains) :: auxhist4_end
integer :: io_form_auxhist4
integer , DIMENSION(max_domains) :: frames_per_auxhist4
character*256 :: auxhist5_inname
character*256 :: auxhist5_outname
integer , DIMENSION(max_domains) :: auxhist5_interval_y
integer , DIMENSION(max_domains) :: auxhist5_interval_d
integer , DIMENSION(max_domains) :: auxhist5_interval_h
integer , DIMENSION(max_domains) :: auxhist5_interval_m
integer , DIMENSION(max_domains) :: auxhist5_interval_s
integer , DIMENSION(max_domains) :: auxhist5_interval
integer , DIMENSION(max_domains) :: auxhist5_begin_y
integer , DIMENSION(max_domains) :: auxhist5_begin_d
integer , DIMENSION(max_domains) :: auxhist5_begin_h
integer , DIMENSION(max_domains) :: auxhist5_begin_m
integer , DIMENSION(max_domains) :: auxhist5_begin_s
integer , DIMENSION(max_domains) :: auxhist5_begin
integer , DIMENSION(max_domains) :: auxhist5_end_y
integer , DIMENSION(max_domains) :: auxhist5_end_d
integer , DIMENSION(max_domains) :: auxhist5_end_h
integer , DIMENSION(max_domains) :: auxhist5_end_m
integer , DIMENSION(max_domains) :: auxhist5_end_s
integer , DIMENSION(max_domains) :: auxhist5_end
integer :: io_form_auxhist5
integer , DIMENSION(max_domains) :: frames_per_auxhist5
character*256 :: auxhist6_inname
character*256 :: auxhist6_outname
integer , DIMENSION(max_domains) :: auxhist6_interval_y
integer , DIMENSION(max_domains) :: auxhist6_interval_d
integer , DIMENSION(max_domains) :: auxhist6_interval_h
integer , DIMENSION(max_domains) :: auxhist6_interval_m
integer , DIMENSION(max_domains) :: auxhist6_interval_s
integer , DIMENSION(max_domains) :: auxhist6_interval
integer , DIMENSION(max_domains) :: auxhist6_begin_y
integer , DIMENSION(max_domains) :: auxhist6_begin_d
integer , DIMENSION(max_domains) :: auxhist6_begin_h
integer , DIMENSION(max_domains) :: auxhist6_begin_m
integer , DIMENSION(max_domains) :: auxhist6_begin_s
integer , DIMENSION(max_domains) :: auxhist6_begin
integer , DIMENSION(max_domains) :: auxhist6_end_y
integer , DIMENSION(max_domains) :: auxhist6_end_d
integer , DIMENSION(max_domains) :: auxhist6_end_h
integer , DIMENSION(max_domains) :: auxhist6_end_m
integer , DIMENSION(max_domains) :: auxhist6_end_s
integer , DIMENSION(max_domains) :: auxhist6_end
integer :: io_form_auxhist6
integer , DIMENSION(max_domains) :: frames_per_auxhist6
character*256 :: auxhist7_inname
character*256 :: auxhist7_outname
integer , DIMENSION(max_domains) :: auxhist7_interval_y
integer , DIMENSION(max_domains) :: auxhist7_interval_d
integer , DIMENSION(max_domains) :: auxhist7_interval_h
integer , DIMENSION(max_domains) :: auxhist7_interval_m
integer , DIMENSION(max_domains) :: auxhist7_interval_s
integer , DIMENSION(max_domains) :: auxhist7_interval
integer , DIMENSION(max_domains) :: auxhist7_begin_y
integer , DIMENSION(max_domains) :: auxhist7_begin_d
integer , DIMENSION(max_domains) :: auxhist7_begin_h
integer , DIMENSION(max_domains) :: auxhist7_begin_m
integer , DIMENSION(max_domains) :: auxhist7_begin_s
integer , DIMENSION(max_domains) :: auxhist7_begin
integer , DIMENSION(max_domains) :: auxhist7_end_y
integer , DIMENSION(max_domains) :: auxhist7_end_d
integer , DIMENSION(max_domains) :: auxhist7_end_h
integer , DIMENSION(max_domains) :: auxhist7_end_m
integer , DIMENSION(max_domains) :: auxhist7_end_s
integer , DIMENSION(max_domains) :: auxhist7_end
integer :: io_form_auxhist7
integer , DIMENSION(max_domains) :: frames_per_auxhist7
character*256 :: auxhist8_inname
character*256 :: auxhist8_outname
integer , DIMENSION(max_domains) :: auxhist8_interval_y
integer , DIMENSION(max_domains) :: auxhist8_interval_d
integer , DIMENSION(max_domains) :: auxhist8_interval_h
integer , DIMENSION(max_domains) :: auxhist8_interval_m
integer , DIMENSION(max_domains) :: auxhist8_interval_s
integer , DIMENSION(max_domains) :: auxhist8_interval
integer , DIMENSION(max_domains) :: auxhist8_begin_y
integer , DIMENSION(max_domains) :: auxhist8_begin_d
integer , DIMENSION(max_domains) :: auxhist8_begin_h
integer , DIMENSION(max_domains) :: auxhist8_begin_m
integer , DIMENSION(max_domains) :: auxhist8_begin_s
integer , DIMENSION(max_domains) :: auxhist8_begin
integer , DIMENSION(max_domains) :: auxhist8_end_y
integer , DIMENSION(max_domains) :: auxhist8_end_d
integer , DIMENSION(max_domains) :: auxhist8_end_h
integer , DIMENSION(max_domains) :: auxhist8_end_m
integer , DIMENSION(max_domains) :: auxhist8_end_s
integer , DIMENSION(max_domains) :: auxhist8_end
integer :: io_form_auxhist8
integer , DIMENSION(max_domains) :: frames_per_auxhist8
character*256 :: auxhist9_inname
character*256 :: auxhist9_outname
integer , DIMENSION(max_domains) :: auxhist9_interval_y
integer , DIMENSION(max_domains) :: auxhist9_interval_d
integer , DIMENSION(max_domains) :: auxhist9_interval_h
integer , DIMENSION(max_domains) :: auxhist9_interval_m
integer , DIMENSION(max_domains) :: auxhist9_interval_s
integer , DIMENSION(max_domains) :: auxhist9_interval
integer , DIMENSION(max_domains) :: auxhist9_begin_y
integer , DIMENSION(max_domains) :: auxhist9_begin_d
integer , DIMENSION(max_domains) :: auxhist9_begin_h
integer , DIMENSION(max_domains) :: auxhist9_begin_m
integer , DIMENSION(max_domains) :: auxhist9_begin_s
integer , DIMENSION(max_domains) :: auxhist9_begin
integer , DIMENSION(max_domains) :: auxhist9_end_y
integer , DIMENSION(max_domains) :: auxhist9_end_d
integer , DIMENSION(max_domains) :: auxhist9_end_h
integer , DIMENSION(max_domains) :: auxhist9_end_m
integer , DIMENSION(max_domains) :: auxhist9_end_s
integer , DIMENSION(max_domains) :: auxhist9_end
integer :: io_form_auxhist9
integer , DIMENSION(max_domains) :: frames_per_auxhist9
character*256 :: auxhist10_inname
character*256 :: auxhist10_outname
integer , DIMENSION(max_domains) :: auxhist10_interval_y
integer , DIMENSION(max_domains) :: auxhist10_interval_d
integer , DIMENSION(max_domains) :: auxhist10_interval_h
integer , DIMENSION(max_domains) :: auxhist10_interval_m
integer , DIMENSION(max_domains) :: auxhist10_interval_s
integer , DIMENSION(max_domains) :: auxhist10_interval
integer , DIMENSION(max_domains) :: auxhist10_begin_y
integer , DIMENSION(max_domains) :: auxhist10_begin_d
integer , DIMENSION(max_domains) :: auxhist10_begin_h
integer , DIMENSION(max_domains) :: auxhist10_begin_m
integer , DIMENSION(max_domains) :: auxhist10_begin_s
integer , DIMENSION(max_domains) :: auxhist10_begin
integer , DIMENSION(max_domains) :: auxhist10_end_y
integer , DIMENSION(max_domains) :: auxhist10_end_d
integer , DIMENSION(max_domains) :: auxhist10_end_h
integer , DIMENSION(max_domains) :: auxhist10_end_m
integer , DIMENSION(max_domains) :: auxhist10_end_s
integer , DIMENSION(max_domains) :: auxhist10_end
integer :: io_form_auxhist10
integer , DIMENSION(max_domains) :: frames_per_auxhist10
character*256 :: auxhist11_inname
character*256 :: auxhist11_outname
integer , DIMENSION(max_domains) :: auxhist11_interval_y
integer , DIMENSION(max_domains) :: auxhist11_interval_d
integer , DIMENSION(max_domains) :: auxhist11_interval_h
integer , DIMENSION(max_domains) :: auxhist11_interval_m
integer , DIMENSION(max_domains) :: auxhist11_interval_s
integer , DIMENSION(max_domains) :: auxhist11_interval
integer , DIMENSION(max_domains) :: auxhist11_begin_y
integer , DIMENSION(max_domains) :: auxhist11_begin_d
integer , DIMENSION(max_domains) :: auxhist11_begin_h
integer , DIMENSION(max_domains) :: auxhist11_begin_m
integer , DIMENSION(max_domains) :: auxhist11_begin_s
integer , DIMENSION(max_domains) :: auxhist11_begin
integer , DIMENSION(max_domains) :: auxhist11_end_y
integer , DIMENSION(max_domains) :: auxhist11_end_d
integer , DIMENSION(max_domains) :: auxhist11_end_h
integer , DIMENSION(max_domains) :: auxhist11_end_m
integer , DIMENSION(max_domains) :: auxhist11_end_s
integer , DIMENSION(max_domains) :: auxhist11_end
integer :: io_form_auxhist11
integer , DIMENSION(max_domains) :: frames_per_auxhist11
character*256 :: auxhist12_inname
character*256 :: auxhist12_outname
integer , DIMENSION(max_domains) :: auxhist12_interval_y
integer , DIMENSION(max_domains) :: auxhist12_interval_d
integer , DIMENSION(max_domains) :: auxhist12_interval_h
integer , DIMENSION(max_domains) :: auxhist12_interval_m
integer , DIMENSION(max_domains) :: auxhist12_interval_s
integer , DIMENSION(max_domains) :: auxhist12_interval
integer , DIMENSION(max_domains) :: auxhist12_begin_y
integer , DIMENSION(max_domains) :: auxhist12_begin_d
integer , DIMENSION(max_domains) :: auxhist12_begin_h
integer , DIMENSION(max_domains) :: auxhist12_begin_m
integer , DIMENSION(max_domains) :: auxhist12_begin_s
integer , DIMENSION(max_domains) :: auxhist12_begin
integer , DIMENSION(max_domains) :: auxhist12_end_y
integer , DIMENSION(max_domains) :: auxhist12_end_d
integer , DIMENSION(max_domains) :: auxhist12_end_h
integer , DIMENSION(max_domains) :: auxhist12_end_m
integer , DIMENSION(max_domains) :: auxhist12_end_s
integer , DIMENSION(max_domains) :: auxhist12_end
integer :: io_form_auxhist12
integer , DIMENSION(max_domains) :: frames_per_auxhist12
character*256 :: auxhist13_inname
character*256 :: auxhist13_outname
integer , DIMENSION(max_domains) :: auxhist13_interval_y
integer , DIMENSION(max_domains) :: auxhist13_interval_d
integer , DIMENSION(max_domains) :: auxhist13_interval_h
integer , DIMENSION(max_domains) :: auxhist13_interval_m
integer , DIMENSION(max_domains) :: auxhist13_interval_s
integer , DIMENSION(max_domains) :: auxhist13_interval
integer , DIMENSION(max_domains) :: auxhist13_begin_y
integer , DIMENSION(max_domains) :: auxhist13_begin_d
integer , DIMENSION(max_domains) :: auxhist13_begin_h
integer , DIMENSION(max_domains) :: auxhist13_begin_m
integer , DIMENSION(max_domains) :: auxhist13_begin_s
integer , DIMENSION(max_domains) :: auxhist13_begin
integer , DIMENSION(max_domains) :: auxhist13_end_y
integer , DIMENSION(max_domains) :: auxhist13_end_d
integer , DIMENSION(max_domains) :: auxhist13_end_h
integer , DIMENSION(max_domains) :: auxhist13_end_m
integer , DIMENSION(max_domains) :: auxhist13_end_s
integer , DIMENSION(max_domains) :: auxhist13_end
integer :: io_form_auxhist13
integer , DIMENSION(max_domains) :: frames_per_auxhist13
character*256 :: auxhist14_inname
character*256 :: auxhist14_outname
integer , DIMENSION(max_domains) :: auxhist14_interval_y
integer , DIMENSION(max_domains) :: auxhist14_interval_d
integer , DIMENSION(max_domains) :: auxhist14_interval_h
integer , DIMENSION(max_domains) :: auxhist14_interval_m
integer , DIMENSION(max_domains) :: auxhist14_interval_s
integer , DIMENSION(max_domains) :: auxhist14_interval
integer , DIMENSION(max_domains) :: auxhist14_begin_y
integer , DIMENSION(max_domains) :: auxhist14_begin_d
integer , DIMENSION(max_domains) :: auxhist14_begin_h
integer , DIMENSION(max_domains) :: auxhist14_begin_m
integer , DIMENSION(max_domains) :: auxhist14_begin_s
integer , DIMENSION(max_domains) :: auxhist14_begin
integer , DIMENSION(max_domains) :: auxhist14_end_y
integer , DIMENSION(max_domains) :: auxhist14_end_d
integer , DIMENSION(max_domains) :: auxhist14_end_h
integer , DIMENSION(max_domains) :: auxhist14_end_m
integer , DIMENSION(max_domains) :: auxhist14_end_s
integer , DIMENSION(max_domains) :: auxhist14_end
integer :: io_form_auxhist14
integer , DIMENSION(max_domains) :: frames_per_auxhist14
character*256 :: auxhist15_inname
character*256 :: auxhist15_outname
integer , DIMENSION(max_domains) :: auxhist15_interval_y
integer , DIMENSION(max_domains) :: auxhist15_interval_d
integer , DIMENSION(max_domains) :: auxhist15_interval_h
integer , DIMENSION(max_domains) :: auxhist15_interval_m
integer , DIMENSION(max_domains) :: auxhist15_interval_s
integer , DIMENSION(max_domains) :: auxhist15_interval
integer , DIMENSION(max_domains) :: auxhist15_begin_y
integer , DIMENSION(max_domains) :: auxhist15_begin_d
integer , DIMENSION(max_domains) :: auxhist15_begin_h
integer , DIMENSION(max_domains) :: auxhist15_begin_m
integer , DIMENSION(max_domains) :: auxhist15_begin_s
integer , DIMENSION(max_domains) :: auxhist15_begin
integer , DIMENSION(max_domains) :: auxhist15_end_y
integer , DIMENSION(max_domains) :: auxhist15_end_d
integer , DIMENSION(max_domains) :: auxhist15_end_h
integer , DIMENSION(max_domains) :: auxhist15_end_m
integer , DIMENSION(max_domains) :: auxhist15_end_s
integer , DIMENSION(max_domains) :: auxhist15_end
integer :: io_form_auxhist15
integer , DIMENSION(max_domains) :: frames_per_auxhist15
character*256 :: auxhist16_inname
character*256 :: auxhist16_outname
integer , DIMENSION(max_domains) :: auxhist16_interval_y
integer , DIMENSION(max_domains) :: auxhist16_interval_d
integer , DIMENSION(max_domains) :: auxhist16_interval_h
integer , DIMENSION(max_domains) :: auxhist16_interval_m
integer , DIMENSION(max_domains) :: auxhist16_interval_s
integer , DIMENSION(max_domains) :: auxhist16_interval
integer , DIMENSION(max_domains) :: auxhist16_begin_y
integer , DIMENSION(max_domains) :: auxhist16_begin_d
integer , DIMENSION(max_domains) :: auxhist16_begin_h
integer , DIMENSION(max_domains) :: auxhist16_begin_m
integer , DIMENSION(max_domains) :: auxhist16_begin_s
integer , DIMENSION(max_domains) :: auxhist16_begin
integer , DIMENSION(max_domains) :: auxhist16_end_y
integer , DIMENSION(max_domains) :: auxhist16_end_d
integer , DIMENSION(max_domains) :: auxhist16_end_h
integer , DIMENSION(max_domains) :: auxhist16_end_m
integer , DIMENSION(max_domains) :: auxhist16_end_s
integer , DIMENSION(max_domains) :: auxhist16_end
integer :: io_form_auxhist16
integer , DIMENSION(max_domains) :: frames_per_auxhist16
character*256 :: auxhist17_inname
character*256 :: auxhist17_outname
integer , DIMENSION(max_domains) :: auxhist17_interval_y
integer , DIMENSION(max_domains) :: auxhist17_interval_d
integer , DIMENSION(max_domains) :: auxhist17_interval_h
integer , DIMENSION(max_domains) :: auxhist17_interval_m
integer , DIMENSION(max_domains) :: auxhist17_interval_s
integer , DIMENSION(max_domains) :: auxhist17_interval
integer , DIMENSION(max_domains) :: auxhist17_begin_y
integer , DIMENSION(max_domains) :: auxhist17_begin_d
integer , DIMENSION(max_domains) :: auxhist17_begin_h
integer , DIMENSION(max_domains) :: auxhist17_begin_m
integer , DIMENSION(max_domains) :: auxhist17_begin_s
integer , DIMENSION(max_domains) :: auxhist17_begin
integer , DIMENSION(max_domains) :: auxhist17_end_y
integer , DIMENSION(max_domains) :: auxhist17_end_d
integer , DIMENSION(max_domains) :: auxhist17_end_h
integer , DIMENSION(max_domains) :: auxhist17_end_m
integer , DIMENSION(max_domains) :: auxhist17_end_s
integer , DIMENSION(max_domains) :: auxhist17_end
integer :: io_form_auxhist17
integer , DIMENSION(max_domains) :: frames_per_auxhist17
character*256 :: auxhist18_inname
character*256 :: auxhist18_outname
integer , DIMENSION(max_domains) :: auxhist18_interval_y
integer , DIMENSION(max_domains) :: auxhist18_interval_d
integer , DIMENSION(max_domains) :: auxhist18_interval_h
integer , DIMENSION(max_domains) :: auxhist18_interval_m
integer , DIMENSION(max_domains) :: auxhist18_interval_s
integer , DIMENSION(max_domains) :: auxhist18_interval
integer , DIMENSION(max_domains) :: auxhist18_begin_y
integer , DIMENSION(max_domains) :: auxhist18_begin_d
integer , DIMENSION(max_domains) :: auxhist18_begin_h
integer , DIMENSION(max_domains) :: auxhist18_begin_m
integer , DIMENSION(max_domains) :: auxhist18_begin_s
integer , DIMENSION(max_domains) :: auxhist18_begin
integer , DIMENSION(max_domains) :: auxhist18_end_y
integer , DIMENSION(max_domains) :: auxhist18_end_d
integer , DIMENSION(max_domains) :: auxhist18_end_h
integer , DIMENSION(max_domains) :: auxhist18_end_m
integer , DIMENSION(max_domains) :: auxhist18_end_s
integer , DIMENSION(max_domains) :: auxhist18_end
integer :: io_form_auxhist18
integer , DIMENSION(max_domains) :: frames_per_auxhist18
character*256 :: auxhist19_inname
character*256 :: auxhist19_outname
integer , DIMENSION(max_domains) :: auxhist19_interval_y
integer , DIMENSION(max_domains) :: auxhist19_interval_d
integer , DIMENSION(max_domains) :: auxhist19_interval_h
integer , DIMENSION(max_domains) :: auxhist19_interval_m
integer , DIMENSION(max_domains) :: auxhist19_interval_s
integer , DIMENSION(max_domains) :: auxhist19_interval
integer , DIMENSION(max_domains) :: auxhist19_begin_y
integer , DIMENSION(max_domains) :: auxhist19_begin_d
integer , DIMENSION(max_domains) :: auxhist19_begin_h
integer , DIMENSION(max_domains) :: auxhist19_begin_m
integer , DIMENSION(max_domains) :: auxhist19_begin_s
integer , DIMENSION(max_domains) :: auxhist19_begin
integer , DIMENSION(max_domains) :: auxhist19_end_y
integer , DIMENSION(max_domains) :: auxhist19_end_d
integer , DIMENSION(max_domains) :: auxhist19_end_h
integer , DIMENSION(max_domains) :: auxhist19_end_m
integer , DIMENSION(max_domains) :: auxhist19_end_s
integer , DIMENSION(max_domains) :: auxhist19_end
integer :: io_form_auxhist19
integer , DIMENSION(max_domains) :: frames_per_auxhist19
character*256 :: auxhist20_inname
character*256 :: auxhist20_outname
integer , DIMENSION(max_domains) :: auxhist20_interval_y
integer , DIMENSION(max_domains) :: auxhist20_interval_d
integer , DIMENSION(max_domains) :: auxhist20_interval_h
integer , DIMENSION(max_domains) :: auxhist20_interval_m
integer , DIMENSION(max_domains) :: auxhist20_interval_s
integer , DIMENSION(max_domains) :: auxhist20_interval
integer , DIMENSION(max_domains) :: auxhist20_begin_y
integer , DIMENSION(max_domains) :: auxhist20_begin_d
integer , DIMENSION(max_domains) :: auxhist20_begin_h
integer , DIMENSION(max_domains) :: auxhist20_begin_m
integer , DIMENSION(max_domains) :: auxhist20_begin_s
integer , DIMENSION(max_domains) :: auxhist20_begin
integer , DIMENSION(max_domains) :: auxhist20_end_y
integer , DIMENSION(max_domains) :: auxhist20_end_d
integer , DIMENSION(max_domains) :: auxhist20_end_h
integer , DIMENSION(max_domains) :: auxhist20_end_m
integer , DIMENSION(max_domains) :: auxhist20_end_s
integer , DIMENSION(max_domains) :: auxhist20_end
integer :: io_form_auxhist20
integer , DIMENSION(max_domains) :: frames_per_auxhist20
character*256 :: auxhist21_inname
character*256 :: auxhist21_outname
integer , DIMENSION(max_domains) :: auxhist21_interval_y
integer , DIMENSION(max_domains) :: auxhist21_interval_d
integer , DIMENSION(max_domains) :: auxhist21_interval_h
integer , DIMENSION(max_domains) :: auxhist21_interval_m
integer , DIMENSION(max_domains) :: auxhist21_interval_s
integer , DIMENSION(max_domains) :: auxhist21_interval
integer , DIMENSION(max_domains) :: auxhist21_begin_y
integer , DIMENSION(max_domains) :: auxhist21_begin_d
integer , DIMENSION(max_domains) :: auxhist21_begin_h
integer , DIMENSION(max_domains) :: auxhist21_begin_m
integer , DIMENSION(max_domains) :: auxhist21_begin_s
integer , DIMENSION(max_domains) :: auxhist21_begin
integer , DIMENSION(max_domains) :: auxhist21_end_y
integer , DIMENSION(max_domains) :: auxhist21_end_d
integer , DIMENSION(max_domains) :: auxhist21_end_h
integer , DIMENSION(max_domains) :: auxhist21_end_m
integer , DIMENSION(max_domains) :: auxhist21_end_s
integer , DIMENSION(max_domains) :: auxhist21_end
integer :: io_form_auxhist21
integer , DIMENSION(max_domains) :: frames_per_auxhist21
character*256 :: auxhist22_inname
character*256 :: auxhist22_outname
integer , DIMENSION(max_domains) :: auxhist22_interval_y
integer , DIMENSION(max_domains) :: auxhist22_interval_d
integer , DIMENSION(max_domains) :: auxhist22_interval_h
integer , DIMENSION(max_domains) :: auxhist22_interval_m
integer , DIMENSION(max_domains) :: auxhist22_interval_s
integer , DIMENSION(max_domains) :: auxhist22_interval
integer , DIMENSION(max_domains) :: auxhist22_begin_y
integer , DIMENSION(max_domains) :: auxhist22_begin_d
integer , DIMENSION(max_domains) :: auxhist22_begin_h
integer , DIMENSION(max_domains) :: auxhist22_begin_m
integer , DIMENSION(max_domains) :: auxhist22_begin_s
integer , DIMENSION(max_domains) :: auxhist22_begin
integer , DIMENSION(max_domains) :: auxhist22_end_y
integer , DIMENSION(max_domains) :: auxhist22_end_d
integer , DIMENSION(max_domains) :: auxhist22_end_h
integer , DIMENSION(max_domains) :: auxhist22_end_m
integer , DIMENSION(max_domains) :: auxhist22_end_s
integer , DIMENSION(max_domains) :: auxhist22_end
integer :: io_form_auxhist22
integer , DIMENSION(max_domains) :: frames_per_auxhist22
character*256 :: auxhist23_inname
character*256 :: auxhist23_outname
integer , DIMENSION(max_domains) :: auxhist23_interval_y
integer , DIMENSION(max_domains) :: auxhist23_interval_d
integer , DIMENSION(max_domains) :: auxhist23_interval_h
integer , DIMENSION(max_domains) :: auxhist23_interval_m
integer , DIMENSION(max_domains) :: auxhist23_interval_s
integer , DIMENSION(max_domains) :: auxhist23_interval
integer , DIMENSION(max_domains) :: auxhist23_begin_y
integer , DIMENSION(max_domains) :: auxhist23_begin_d
integer , DIMENSION(max_domains) :: auxhist23_begin_h
integer , DIMENSION(max_domains) :: auxhist23_begin_m
integer , DIMENSION(max_domains) :: auxhist23_begin_s
integer , DIMENSION(max_domains) :: auxhist23_begin
integer , DIMENSION(max_domains) :: auxhist23_end_y
integer , DIMENSION(max_domains) :: auxhist23_end_d
integer , DIMENSION(max_domains) :: auxhist23_end_h
integer , DIMENSION(max_domains) :: auxhist23_end_m
integer , DIMENSION(max_domains) :: auxhist23_end_s
integer , DIMENSION(max_domains) :: auxhist23_end
integer :: io_form_auxhist23
integer , DIMENSION(max_domains) :: frames_per_auxhist23
character*256 :: auxhist24_inname
character*256 :: auxhist24_outname
integer , DIMENSION(max_domains) :: auxhist24_interval_y
integer , DIMENSION(max_domains) :: auxhist24_interval_d
integer , DIMENSION(max_domains) :: auxhist24_interval_h
integer , DIMENSION(max_domains) :: auxhist24_interval_m
integer , DIMENSION(max_domains) :: auxhist24_interval_s
integer , DIMENSION(max_domains) :: auxhist24_interval
integer , DIMENSION(max_domains) :: auxhist24_begin_y
integer , DIMENSION(max_domains) :: auxhist24_begin_d
integer , DIMENSION(max_domains) :: auxhist24_begin_h
integer , DIMENSION(max_domains) :: auxhist24_begin_m
integer , DIMENSION(max_domains) :: auxhist24_begin_s
integer , DIMENSION(max_domains) :: auxhist24_begin
integer , DIMENSION(max_domains) :: auxhist24_end_y
integer , DIMENSION(max_domains) :: auxhist24_end_d
integer , DIMENSION(max_domains) :: auxhist24_end_h
integer , DIMENSION(max_domains) :: auxhist24_end_m
integer , DIMENSION(max_domains) :: auxhist24_end_s
integer , DIMENSION(max_domains) :: auxhist24_end
integer :: io_form_auxhist24
integer , DIMENSION(max_domains) :: frames_per_auxhist24
character*256 :: auxinput1_outname
integer , DIMENSION(max_domains) :: auxinput1_interval_y
integer , DIMENSION(max_domains) :: auxinput1_interval_d
integer , DIMENSION(max_domains) :: auxinput1_interval_h
integer , DIMENSION(max_domains) :: auxinput1_interval_m
integer , DIMENSION(max_domains) :: auxinput1_interval_s
integer , DIMENSION(max_domains) :: auxinput1_interval
integer , DIMENSION(max_domains) :: auxinput1_begin_y
integer , DIMENSION(max_domains) :: auxinput1_begin_d
integer , DIMENSION(max_domains) :: auxinput1_begin_h
integer , DIMENSION(max_domains) :: auxinput1_begin_m
integer , DIMENSION(max_domains) :: auxinput1_begin_s
integer , DIMENSION(max_domains) :: auxinput1_begin
integer , DIMENSION(max_domains) :: auxinput1_end_y
integer , DIMENSION(max_domains) :: auxinput1_end_d
integer , DIMENSION(max_domains) :: auxinput1_end_h
integer , DIMENSION(max_domains) :: auxinput1_end_m
integer , DIMENSION(max_domains) :: auxinput1_end_s
integer , DIMENSION(max_domains) :: auxinput1_end
integer , DIMENSION(max_domains) :: frames_per_auxinput1
character*256 :: auxinput2_inname
character*256 :: auxinput2_outname
integer , DIMENSION(max_domains) :: auxinput2_interval_y
integer , DIMENSION(max_domains) :: auxinput2_interval_d
integer , DIMENSION(max_domains) :: auxinput2_interval_h
integer , DIMENSION(max_domains) :: auxinput2_interval_m
integer , DIMENSION(max_domains) :: auxinput2_interval_s
integer , DIMENSION(max_domains) :: auxinput2_interval
integer , DIMENSION(max_domains) :: auxinput2_begin_y
integer , DIMENSION(max_domains) :: auxinput2_begin_d
integer , DIMENSION(max_domains) :: auxinput2_begin_h
integer , DIMENSION(max_domains) :: auxinput2_begin_m
integer , DIMENSION(max_domains) :: auxinput2_begin_s
integer , DIMENSION(max_domains) :: auxinput2_begin
integer , DIMENSION(max_domains) :: auxinput2_end_y
integer , DIMENSION(max_domains) :: auxinput2_end_d
integer , DIMENSION(max_domains) :: auxinput2_end_h
integer , DIMENSION(max_domains) :: auxinput2_end_m
integer , DIMENSION(max_domains) :: auxinput2_end_s
integer , DIMENSION(max_domains) :: auxinput2_end
integer :: io_form_auxinput2
integer , DIMENSION(max_domains) :: frames_per_auxinput2
character*256 :: auxinput3_inname
character*256 :: auxinput3_outname
integer , DIMENSION(max_domains) :: auxinput3_interval_y
integer , DIMENSION(max_domains) :: auxinput3_interval_d
integer , DIMENSION(max_domains) :: auxinput3_interval_h
integer , DIMENSION(max_domains) :: auxinput3_interval_m
integer , DIMENSION(max_domains) :: auxinput3_interval_s
integer , DIMENSION(max_domains) :: auxinput3_interval
integer , DIMENSION(max_domains) :: auxinput3_begin_y
integer , DIMENSION(max_domains) :: auxinput3_begin_d
integer , DIMENSION(max_domains) :: auxinput3_begin_h
integer , DIMENSION(max_domains) :: auxinput3_begin_m
integer , DIMENSION(max_domains) :: auxinput3_begin_s
integer , DIMENSION(max_domains) :: auxinput3_begin
integer , DIMENSION(max_domains) :: auxinput3_end_y
integer , DIMENSION(max_domains) :: auxinput3_end_d
integer , DIMENSION(max_domains) :: auxinput3_end_h
integer , DIMENSION(max_domains) :: auxinput3_end_m
integer , DIMENSION(max_domains) :: auxinput3_end_s
integer , DIMENSION(max_domains) :: auxinput3_end
integer :: io_form_auxinput3
integer , DIMENSION(max_domains) :: frames_per_auxinput3
character*256 :: auxinput4_inname
character*256 :: auxinput4_outname
integer , DIMENSION(max_domains) :: auxinput4_interval_y
integer , DIMENSION(max_domains) :: auxinput4_interval_d
integer , DIMENSION(max_domains) :: auxinput4_interval_h
integer , DIMENSION(max_domains) :: auxinput4_interval_m
integer , DIMENSION(max_domains) :: auxinput4_interval_s
integer , DIMENSION(max_domains) :: auxinput4_interval
integer , DIMENSION(max_domains) :: auxinput4_begin_y
integer , DIMENSION(max_domains) :: auxinput4_begin_d
integer , DIMENSION(max_domains) :: auxinput4_begin_h
integer , DIMENSION(max_domains) :: auxinput4_begin_m
integer , DIMENSION(max_domains) :: auxinput4_begin_s
integer , DIMENSION(max_domains) :: auxinput4_begin
integer , DIMENSION(max_domains) :: auxinput4_end_y
integer , DIMENSION(max_domains) :: auxinput4_end_d
integer , DIMENSION(max_domains) :: auxinput4_end_h
integer , DIMENSION(max_domains) :: auxinput4_end_m
integer , DIMENSION(max_domains) :: auxinput4_end_s
integer , DIMENSION(max_domains) :: auxinput4_end
integer :: io_form_auxinput4
integer , DIMENSION(max_domains) :: frames_per_auxinput4
character*256 :: auxinput5_inname
character*256 :: auxinput5_outname
integer , DIMENSION(max_domains) :: auxinput5_interval_y
integer , DIMENSION(max_domains) :: auxinput5_interval_d
integer , DIMENSION(max_domains) :: auxinput5_interval_h
integer , DIMENSION(max_domains) :: auxinput5_interval_m
integer , DIMENSION(max_domains) :: auxinput5_interval_s
integer , DIMENSION(max_domains) :: auxinput5_interval
integer , DIMENSION(max_domains) :: auxinput5_begin_y
integer , DIMENSION(max_domains) :: auxinput5_begin_d
integer , DIMENSION(max_domains) :: auxinput5_begin_h
integer , DIMENSION(max_domains) :: auxinput5_begin_m
integer , DIMENSION(max_domains) :: auxinput5_begin_s
integer , DIMENSION(max_domains) :: auxinput5_begin
integer , DIMENSION(max_domains) :: auxinput5_end_y
integer , DIMENSION(max_domains) :: auxinput5_end_d
integer , DIMENSION(max_domains) :: auxinput5_end_h
integer , DIMENSION(max_domains) :: auxinput5_end_m
integer , DIMENSION(max_domains) :: auxinput5_end_s
integer , DIMENSION(max_domains) :: auxinput5_end
integer :: io_form_auxinput5
integer , DIMENSION(max_domains) :: frames_per_auxinput5
character*256 :: auxinput6_inname
character*256 :: auxinput6_outname
integer , DIMENSION(max_domains) :: auxinput6_interval_y
integer , DIMENSION(max_domains) :: auxinput6_interval_d
integer , DIMENSION(max_domains) :: auxinput6_interval_h
integer , DIMENSION(max_domains) :: auxinput6_interval_m
integer , DIMENSION(max_domains) :: auxinput6_interval_s
integer , DIMENSION(max_domains) :: auxinput6_interval
integer , DIMENSION(max_domains) :: auxinput6_begin_y
integer , DIMENSION(max_domains) :: auxinput6_begin_d
integer , DIMENSION(max_domains) :: auxinput6_begin_h
integer , DIMENSION(max_domains) :: auxinput6_begin_m
integer , DIMENSION(max_domains) :: auxinput6_begin_s
integer , DIMENSION(max_domains) :: auxinput6_begin
integer , DIMENSION(max_domains) :: auxinput6_end_y
integer , DIMENSION(max_domains) :: auxinput6_end_d
integer , DIMENSION(max_domains) :: auxinput6_end_h
integer , DIMENSION(max_domains) :: auxinput6_end_m
integer , DIMENSION(max_domains) :: auxinput6_end_s
integer , DIMENSION(max_domains) :: auxinput6_end
integer :: io_form_auxinput6
integer , DIMENSION(max_domains) :: frames_per_auxinput6
character*256 :: auxinput7_inname
character*256 :: auxinput7_outname
integer , DIMENSION(max_domains) :: auxinput7_interval_y
integer , DIMENSION(max_domains) :: auxinput7_interval_d
integer , DIMENSION(max_domains) :: auxinput7_interval_h
integer , DIMENSION(max_domains) :: auxinput7_interval_m
integer , DIMENSION(max_domains) :: auxinput7_interval_s
integer , DIMENSION(max_domains) :: auxinput7_interval
integer , DIMENSION(max_domains) :: auxinput7_begin_y
integer , DIMENSION(max_domains) :: auxinput7_begin_d
integer , DIMENSION(max_domains) :: auxinput7_begin_h
integer , DIMENSION(max_domains) :: auxinput7_begin_m
integer , DIMENSION(max_domains) :: auxinput7_begin_s
integer , DIMENSION(max_domains) :: auxinput7_begin
integer , DIMENSION(max_domains) :: auxinput7_end_y
integer , DIMENSION(max_domains) :: auxinput7_end_d
integer , DIMENSION(max_domains) :: auxinput7_end_h
integer , DIMENSION(max_domains) :: auxinput7_end_m
integer , DIMENSION(max_domains) :: auxinput7_end_s
integer , DIMENSION(max_domains) :: auxinput7_end
integer :: io_form_auxinput7
integer , DIMENSION(max_domains) :: frames_per_auxinput7
character*256 :: auxinput8_inname
character*256 :: auxinput8_outname
integer , DIMENSION(max_domains) :: auxinput8_interval_y
integer , DIMENSION(max_domains) :: auxinput8_interval_d
integer , DIMENSION(max_domains) :: auxinput8_interval_h
integer , DIMENSION(max_domains) :: auxinput8_interval_m
integer , DIMENSION(max_domains) :: auxinput8_interval_s
integer , DIMENSION(max_domains) :: auxinput8_interval
integer , DIMENSION(max_domains) :: auxinput8_begin_y
integer , DIMENSION(max_domains) :: auxinput8_begin_d
integer , DIMENSION(max_domains) :: auxinput8_begin_h
integer , DIMENSION(max_domains) :: auxinput8_begin_m
integer , DIMENSION(max_domains) :: auxinput8_begin_s
integer , DIMENSION(max_domains) :: auxinput8_begin
integer , DIMENSION(max_domains) :: auxinput8_end_y
integer , DIMENSION(max_domains) :: auxinput8_end_d
integer , DIMENSION(max_domains) :: auxinput8_end_h
integer , DIMENSION(max_domains) :: auxinput8_end_m
integer , DIMENSION(max_domains) :: auxinput8_end_s
integer , DIMENSION(max_domains) :: auxinput8_end
integer :: io_form_auxinput8
integer , DIMENSION(max_domains) :: frames_per_auxinput8
character*256 :: auxinput9_inname
character*256 :: auxinput9_outname
integer , DIMENSION(max_domains) :: auxinput9_interval_y
integer , DIMENSION(max_domains) :: auxinput9_interval_d
integer , DIMENSION(max_domains) :: auxinput9_interval_h
integer , DIMENSION(max_domains) :: auxinput9_interval_m
integer , DIMENSION(max_domains) :: auxinput9_interval_s
integer , DIMENSION(max_domains) :: auxinput9_interval
integer , DIMENSION(max_domains) :: auxinput9_begin_y
integer , DIMENSION(max_domains) :: auxinput9_begin_d
integer , DIMENSION(max_domains) :: auxinput9_begin_h
integer , DIMENSION(max_domains) :: auxinput9_begin_m
integer , DIMENSION(max_domains) :: auxinput9_begin_s
integer , DIMENSION(max_domains) :: auxinput9_begin
integer , DIMENSION(max_domains) :: auxinput9_end_y
integer , DIMENSION(max_domains) :: auxinput9_end_d
integer , DIMENSION(max_domains) :: auxinput9_end_h
integer , DIMENSION(max_domains) :: auxinput9_end_m
integer , DIMENSION(max_domains) :: auxinput9_end_s
integer , DIMENSION(max_domains) :: auxinput9_end
integer :: io_form_auxinput9
integer , DIMENSION(max_domains) :: frames_per_auxinput9
character*256 :: auxinput10_inname
character*256 :: auxinput10_outname
integer , DIMENSION(max_domains) :: auxinput10_interval_y
integer , DIMENSION(max_domains) :: auxinput10_interval_d
integer , DIMENSION(max_domains) :: auxinput10_interval_h
integer , DIMENSION(max_domains) :: auxinput10_interval_m
integer , DIMENSION(max_domains) :: auxinput10_interval_s
integer , DIMENSION(max_domains) :: auxinput10_interval
integer , DIMENSION(max_domains) :: auxinput10_begin_y
integer , DIMENSION(max_domains) :: auxinput10_begin_d
integer , DIMENSION(max_domains) :: auxinput10_begin_h
integer , DIMENSION(max_domains) :: auxinput10_begin_m
integer , DIMENSION(max_domains) :: auxinput10_begin_s
integer , DIMENSION(max_domains) :: auxinput10_begin
integer , DIMENSION(max_domains) :: auxinput10_end_y
integer , DIMENSION(max_domains) :: auxinput10_end_d
integer , DIMENSION(max_domains) :: auxinput10_end_h
integer , DIMENSION(max_domains) :: auxinput10_end_m
integer , DIMENSION(max_domains) :: auxinput10_end_s
integer , DIMENSION(max_domains) :: auxinput10_end
integer :: io_form_auxinput10
integer , DIMENSION(max_domains) :: frames_per_auxinput10
character*256 :: auxinput11_inname
character*256 :: auxinput11_outname
integer , DIMENSION(max_domains) :: auxinput11_interval_y
integer , DIMENSION(max_domains) :: auxinput11_interval_d
integer , DIMENSION(max_domains) :: auxinput11_interval_h
integer , DIMENSION(max_domains) :: auxinput11_interval_m
integer , DIMENSION(max_domains) :: auxinput11_interval_s
integer , DIMENSION(max_domains) :: auxinput11_interval
integer , DIMENSION(max_domains) :: auxinput11_begin_y
integer , DIMENSION(max_domains) :: auxinput11_begin_d
integer , DIMENSION(max_domains) :: auxinput11_begin_h
integer , DIMENSION(max_domains) :: auxinput11_begin_m
integer , DIMENSION(max_domains) :: auxinput11_begin_s
integer , DIMENSION(max_domains) :: auxinput11_begin
integer , DIMENSION(max_domains) :: auxinput11_end_y
integer , DIMENSION(max_domains) :: auxinput11_end_d
integer , DIMENSION(max_domains) :: auxinput11_end_h
integer , DIMENSION(max_domains) :: auxinput11_end_m
integer , DIMENSION(max_domains) :: auxinput11_end_s
integer , DIMENSION(max_domains) :: auxinput11_end
integer :: io_form_auxinput11
integer , DIMENSION(max_domains) :: frames_per_auxinput11
character*256 :: auxinput12_inname
character*256 :: auxinput12_outname
integer , DIMENSION(max_domains) :: auxinput12_interval_y
integer , DIMENSION(max_domains) :: auxinput12_interval_d
integer , DIMENSION(max_domains) :: auxinput12_interval_h
integer , DIMENSION(max_domains) :: auxinput12_interval_m
integer , DIMENSION(max_domains) :: auxinput12_interval_s
integer , DIMENSION(max_domains) :: auxinput12_interval
integer , DIMENSION(max_domains) :: auxinput12_begin_y
integer , DIMENSION(max_domains) :: auxinput12_begin_d
integer , DIMENSION(max_domains) :: auxinput12_begin_h
integer , DIMENSION(max_domains) :: auxinput12_begin_m
integer , DIMENSION(max_domains) :: auxinput12_begin_s
integer , DIMENSION(max_domains) :: auxinput12_begin
integer , DIMENSION(max_domains) :: auxinput12_end_y
integer , DIMENSION(max_domains) :: auxinput12_end_d
integer , DIMENSION(max_domains) :: auxinput12_end_h
integer , DIMENSION(max_domains) :: auxinput12_end_m
integer , DIMENSION(max_domains) :: auxinput12_end_s
integer , DIMENSION(max_domains) :: auxinput12_end
integer :: io_form_auxinput12
integer , DIMENSION(max_domains) :: frames_per_auxinput12
character*256 :: auxinput13_inname
character*256 :: auxinput13_outname
integer , DIMENSION(max_domains) :: auxinput13_interval_y
integer , DIMENSION(max_domains) :: auxinput13_interval_d
integer , DIMENSION(max_domains) :: auxinput13_interval_h
integer , DIMENSION(max_domains) :: auxinput13_interval_m
integer , DIMENSION(max_domains) :: auxinput13_interval_s
integer , DIMENSION(max_domains) :: auxinput13_interval
integer , DIMENSION(max_domains) :: auxinput13_begin_y
integer , DIMENSION(max_domains) :: auxinput13_begin_d
integer , DIMENSION(max_domains) :: auxinput13_begin_h
integer , DIMENSION(max_domains) :: auxinput13_begin_m
integer , DIMENSION(max_domains) :: auxinput13_begin_s
integer , DIMENSION(max_domains) :: auxinput13_begin
integer , DIMENSION(max_domains) :: auxinput13_end_y
integer , DIMENSION(max_domains) :: auxinput13_end_d
integer , DIMENSION(max_domains) :: auxinput13_end_h
integer , DIMENSION(max_domains) :: auxinput13_end_m
integer , DIMENSION(max_domains) :: auxinput13_end_s
integer , DIMENSION(max_domains) :: auxinput13_end
integer :: io_form_auxinput13
integer , DIMENSION(max_domains) :: frames_per_auxinput13
character*256 :: auxinput14_inname
character*256 :: auxinput14_outname
integer , DIMENSION(max_domains) :: auxinput14_interval_y
integer , DIMENSION(max_domains) :: auxinput14_interval_d
integer , DIMENSION(max_domains) :: auxinput14_interval_h
integer , DIMENSION(max_domains) :: auxinput14_interval_m
integer , DIMENSION(max_domains) :: auxinput14_interval_s
integer , DIMENSION(max_domains) :: auxinput14_interval
integer , DIMENSION(max_domains) :: auxinput14_begin_y
integer , DIMENSION(max_domains) :: auxinput14_begin_d
integer , DIMENSION(max_domains) :: auxinput14_begin_h
integer , DIMENSION(max_domains) :: auxinput14_begin_m
integer , DIMENSION(max_domains) :: auxinput14_begin_s
integer , DIMENSION(max_domains) :: auxinput14_begin
integer , DIMENSION(max_domains) :: auxinput14_end_y
integer , DIMENSION(max_domains) :: auxinput14_end_d
integer , DIMENSION(max_domains) :: auxinput14_end_h
integer , DIMENSION(max_domains) :: auxinput14_end_m
integer , DIMENSION(max_domains) :: auxinput14_end_s
integer , DIMENSION(max_domains) :: auxinput14_end
integer :: io_form_auxinput14
integer , DIMENSION(max_domains) :: frames_per_auxinput14
character*256 :: auxinput15_inname
character*256 :: auxinput15_outname
integer , DIMENSION(max_domains) :: auxinput15_interval_y
integer , DIMENSION(max_domains) :: auxinput15_interval_d
integer , DIMENSION(max_domains) :: auxinput15_interval_h
integer , DIMENSION(max_domains) :: auxinput15_interval_m
integer , DIMENSION(max_domains) :: auxinput15_interval_s
integer , DIMENSION(max_domains) :: auxinput15_interval
integer , DIMENSION(max_domains) :: auxinput15_begin_y
integer , DIMENSION(max_domains) :: auxinput15_begin_d
integer , DIMENSION(max_domains) :: auxinput15_begin_h
integer , DIMENSION(max_domains) :: auxinput15_begin_m
integer , DIMENSION(max_domains) :: auxinput15_begin_s
integer , DIMENSION(max_domains) :: auxinput15_begin
integer , DIMENSION(max_domains) :: auxinput15_end_y
integer , DIMENSION(max_domains) :: auxinput15_end_d
integer , DIMENSION(max_domains) :: auxinput15_end_h
integer , DIMENSION(max_domains) :: auxinput15_end_m
integer , DIMENSION(max_domains) :: auxinput15_end_s
integer , DIMENSION(max_domains) :: auxinput15_end
integer :: io_form_auxinput15
integer , DIMENSION(max_domains) :: frames_per_auxinput15
character*256 :: auxinput16_inname
character*256 :: auxinput16_outname
integer , DIMENSION(max_domains) :: auxinput16_interval_y
integer , DIMENSION(max_domains) :: auxinput16_interval_d
integer , DIMENSION(max_domains) :: auxinput16_interval_h
integer , DIMENSION(max_domains) :: auxinput16_interval_m
integer , DIMENSION(max_domains) :: auxinput16_interval_s
integer , DIMENSION(max_domains) :: auxinput16_interval
integer , DIMENSION(max_domains) :: auxinput16_begin_y
integer , DIMENSION(max_domains) :: auxinput16_begin_d
integer , DIMENSION(max_domains) :: auxinput16_begin_h
integer , DIMENSION(max_domains) :: auxinput16_begin_m
integer , DIMENSION(max_domains) :: auxinput16_begin_s
integer , DIMENSION(max_domains) :: auxinput16_begin
integer , DIMENSION(max_domains) :: auxinput16_end_y
integer , DIMENSION(max_domains) :: auxinput16_end_d
integer , DIMENSION(max_domains) :: auxinput16_end_h
integer , DIMENSION(max_domains) :: auxinput16_end_m
integer , DIMENSION(max_domains) :: auxinput16_end_s
integer , DIMENSION(max_domains) :: auxinput16_end
integer :: io_form_auxinput16
integer , DIMENSION(max_domains) :: frames_per_auxinput16
character*256 :: auxinput17_inname
character*256 :: auxinput17_outname
integer , DIMENSION(max_domains) :: auxinput17_interval_y
integer , DIMENSION(max_domains) :: auxinput17_interval_d
integer , DIMENSION(max_domains) :: auxinput17_interval_h
integer , DIMENSION(max_domains) :: auxinput17_interval_m
integer , DIMENSION(max_domains) :: auxinput17_interval_s
integer , DIMENSION(max_domains) :: auxinput17_interval
integer , DIMENSION(max_domains) :: auxinput17_begin_y
integer , DIMENSION(max_domains) :: auxinput17_begin_d
integer , DIMENSION(max_domains) :: auxinput17_begin_h
integer , DIMENSION(max_domains) :: auxinput17_begin_m
integer , DIMENSION(max_domains) :: auxinput17_begin_s
integer , DIMENSION(max_domains) :: auxinput17_begin
integer , DIMENSION(max_domains) :: auxinput17_end_y
integer , DIMENSION(max_domains) :: auxinput17_end_d
integer , DIMENSION(max_domains) :: auxinput17_end_h
integer , DIMENSION(max_domains) :: auxinput17_end_m
integer , DIMENSION(max_domains) :: auxinput17_end_s
integer , DIMENSION(max_domains) :: auxinput17_end
integer :: io_form_auxinput17
integer , DIMENSION(max_domains) :: frames_per_auxinput17
character*256 :: auxinput18_inname
character*256 :: auxinput18_outname
integer , DIMENSION(max_domains) :: auxinput18_interval_y
integer , DIMENSION(max_domains) :: auxinput18_interval_d
integer , DIMENSION(max_domains) :: auxinput18_interval_h
integer , DIMENSION(max_domains) :: auxinput18_interval_m
integer , DIMENSION(max_domains) :: auxinput18_interval_s
integer , DIMENSION(max_domains) :: auxinput18_interval
integer , DIMENSION(max_domains) :: auxinput18_begin_y
integer , DIMENSION(max_domains) :: auxinput18_begin_d
integer , DIMENSION(max_domains) :: auxinput18_begin_h
integer , DIMENSION(max_domains) :: auxinput18_begin_m
integer , DIMENSION(max_domains) :: auxinput18_begin_s
integer , DIMENSION(max_domains) :: auxinput18_begin
integer , DIMENSION(max_domains) :: auxinput18_end_y
integer , DIMENSION(max_domains) :: auxinput18_end_d
integer , DIMENSION(max_domains) :: auxinput18_end_h
integer , DIMENSION(max_domains) :: auxinput18_end_m
integer , DIMENSION(max_domains) :: auxinput18_end_s
integer , DIMENSION(max_domains) :: auxinput18_end
integer :: io_form_auxinput18
integer , DIMENSION(max_domains) :: frames_per_auxinput18
character*256 :: auxinput19_inname
character*256 :: auxinput19_outname
integer , DIMENSION(max_domains) :: auxinput19_interval_y
integer , DIMENSION(max_domains) :: auxinput19_interval_d
integer , DIMENSION(max_domains) :: auxinput19_interval_h
integer , DIMENSION(max_domains) :: auxinput19_interval_m
integer , DIMENSION(max_domains) :: auxinput19_interval_s
integer , DIMENSION(max_domains) :: auxinput19_interval
integer , DIMENSION(max_domains) :: auxinput19_begin_y
integer , DIMENSION(max_domains) :: auxinput19_begin_d
integer , DIMENSION(max_domains) :: auxinput19_begin_h
integer , DIMENSION(max_domains) :: auxinput19_begin_m
integer , DIMENSION(max_domains) :: auxinput19_begin_s
integer , DIMENSION(max_domains) :: auxinput19_begin
integer , DIMENSION(max_domains) :: auxinput19_end_y
integer , DIMENSION(max_domains) :: auxinput19_end_d
integer , DIMENSION(max_domains) :: auxinput19_end_h
integer , DIMENSION(max_domains) :: auxinput19_end_m
integer , DIMENSION(max_domains) :: auxinput19_end_s
integer , DIMENSION(max_domains) :: auxinput19_end
integer :: io_form_auxinput19
integer , DIMENSION(max_domains) :: frames_per_auxinput19
character*256 :: auxinput20_inname
character*256 :: auxinput20_outname
integer , DIMENSION(max_domains) :: auxinput20_interval_y
integer , DIMENSION(max_domains) :: auxinput20_interval_d
integer , DIMENSION(max_domains) :: auxinput20_interval_h
integer , DIMENSION(max_domains) :: auxinput20_interval_m
integer , DIMENSION(max_domains) :: auxinput20_interval_s
integer , DIMENSION(max_domains) :: auxinput20_interval
integer , DIMENSION(max_domains) :: auxinput20_begin_y
integer , DIMENSION(max_domains) :: auxinput20_begin_d
integer , DIMENSION(max_domains) :: auxinput20_begin_h
integer , DIMENSION(max_domains) :: auxinput20_begin_m
integer , DIMENSION(max_domains) :: auxinput20_begin_s
integer , DIMENSION(max_domains) :: auxinput20_begin
integer , DIMENSION(max_domains) :: auxinput20_end_y
integer , DIMENSION(max_domains) :: auxinput20_end_d
integer , DIMENSION(max_domains) :: auxinput20_end_h
integer , DIMENSION(max_domains) :: auxinput20_end_m
integer , DIMENSION(max_domains) :: auxinput20_end_s
integer , DIMENSION(max_domains) :: auxinput20_end
integer :: io_form_auxinput20
integer , DIMENSION(max_domains) :: frames_per_auxinput20
character*256 :: auxinput21_inname
character*256 :: auxinput21_outname
integer , DIMENSION(max_domains) :: auxinput21_interval_y
integer , DIMENSION(max_domains) :: auxinput21_interval_d
integer , DIMENSION(max_domains) :: auxinput21_interval_h
integer , DIMENSION(max_domains) :: auxinput21_interval_m
integer , DIMENSION(max_domains) :: auxinput21_interval_s
integer , DIMENSION(max_domains) :: auxinput21_interval
integer , DIMENSION(max_domains) :: auxinput21_begin_y
integer , DIMENSION(max_domains) :: auxinput21_begin_d
integer , DIMENSION(max_domains) :: auxinput21_begin_h
integer , DIMENSION(max_domains) :: auxinput21_begin_m
integer , DIMENSION(max_domains) :: auxinput21_begin_s
integer , DIMENSION(max_domains) :: auxinput21_begin
integer , DIMENSION(max_domains) :: auxinput21_end_y
integer , DIMENSION(max_domains) :: auxinput21_end_d
integer , DIMENSION(max_domains) :: auxinput21_end_h
integer , DIMENSION(max_domains) :: auxinput21_end_m
integer , DIMENSION(max_domains) :: auxinput21_end_s
integer , DIMENSION(max_domains) :: auxinput21_end
integer :: io_form_auxinput21
integer , DIMENSION(max_domains) :: frames_per_auxinput21
character*256 :: auxinput22_inname
character*256 :: auxinput22_outname
integer , DIMENSION(max_domains) :: auxinput22_interval_y
integer , DIMENSION(max_domains) :: auxinput22_interval_d
integer , DIMENSION(max_domains) :: auxinput22_interval_h
integer , DIMENSION(max_domains) :: auxinput22_interval_m
integer , DIMENSION(max_domains) :: auxinput22_interval_s
integer , DIMENSION(max_domains) :: auxinput22_interval
integer , DIMENSION(max_domains) :: auxinput22_begin_y
integer , DIMENSION(max_domains) :: auxinput22_begin_d
integer , DIMENSION(max_domains) :: auxinput22_begin_h
integer , DIMENSION(max_domains) :: auxinput22_begin_m
integer , DIMENSION(max_domains) :: auxinput22_begin_s
integer , DIMENSION(max_domains) :: auxinput22_begin
integer , DIMENSION(max_domains) :: auxinput22_end_y
integer , DIMENSION(max_domains) :: auxinput22_end_d
integer , DIMENSION(max_domains) :: auxinput22_end_h
integer , DIMENSION(max_domains) :: auxinput22_end_m
integer , DIMENSION(max_domains) :: auxinput22_end_s
integer , DIMENSION(max_domains) :: auxinput22_end
integer :: io_form_auxinput22
integer , DIMENSION(max_domains) :: frames_per_auxinput22
character*256 :: auxinput23_inname
character*256 :: auxinput23_outname
integer , DIMENSION(max_domains) :: auxinput23_interval_y
integer , DIMENSION(max_domains) :: auxinput23_interval_d
integer , DIMENSION(max_domains) :: auxinput23_interval_h
integer , DIMENSION(max_domains) :: auxinput23_interval_m
integer , DIMENSION(max_domains) :: auxinput23_interval_s
integer , DIMENSION(max_domains) :: auxinput23_interval
integer , DIMENSION(max_domains) :: auxinput23_begin_y
integer , DIMENSION(max_domains) :: auxinput23_begin_d
integer , DIMENSION(max_domains) :: auxinput23_begin_h
integer , DIMENSION(max_domains) :: auxinput23_begin_m
integer , DIMENSION(max_domains) :: auxinput23_begin_s
integer , DIMENSION(max_domains) :: auxinput23_begin
integer , DIMENSION(max_domains) :: auxinput23_end_y
integer , DIMENSION(max_domains) :: auxinput23_end_d
integer , DIMENSION(max_domains) :: auxinput23_end_h
integer , DIMENSION(max_domains) :: auxinput23_end_m
integer , DIMENSION(max_domains) :: auxinput23_end_s
integer , DIMENSION(max_domains) :: auxinput23_end
integer :: io_form_auxinput23
integer , DIMENSION(max_domains) :: frames_per_auxinput23
character*256 :: auxinput24_inname
character*256 :: auxinput24_outname
integer , DIMENSION(max_domains) :: auxinput24_interval_y
integer , DIMENSION(max_domains) :: auxinput24_interval_d
integer , DIMENSION(max_domains) :: auxinput24_interval_h
integer , DIMENSION(max_domains) :: auxinput24_interval_m
integer , DIMENSION(max_domains) :: auxinput24_interval_s
integer , DIMENSION(max_domains) :: auxinput24_interval
integer , DIMENSION(max_domains) :: auxinput24_begin_y
integer , DIMENSION(max_domains) :: auxinput24_begin_d
integer , DIMENSION(max_domains) :: auxinput24_begin_h
integer , DIMENSION(max_domains) :: auxinput24_begin_m
integer , DIMENSION(max_domains) :: auxinput24_begin_s
integer , DIMENSION(max_domains) :: auxinput24_begin
integer , DIMENSION(max_domains) :: auxinput24_end_y
integer , DIMENSION(max_domains) :: auxinput24_end_d
integer , DIMENSION(max_domains) :: auxinput24_end_h
integer , DIMENSION(max_domains) :: auxinput24_end_m
integer , DIMENSION(max_domains) :: auxinput24_end_s
integer , DIMENSION(max_domains) :: auxinput24_end
integer :: io_form_auxinput24
integer , DIMENSION(max_domains) :: frames_per_auxinput24
integer , DIMENSION(max_domains) :: history_interval
integer , DIMENSION(max_domains) :: frames_per_outfile
logical :: restart
integer :: restart_interval
integer :: io_form_input
integer :: io_form_history
integer :: io_form_restart
integer :: io_form_boundary
integer :: debug_level
logical :: self_test_domain
character*256 :: history_outname
character*256 :: history_inname
logical :: use_netcdf_classic
integer , DIMENSION(max_domains) :: history_interval_d
integer , DIMENSION(max_domains) :: history_interval_h
integer , DIMENSION(max_domains) :: history_interval_m
integer , DIMENSION(max_domains) :: history_interval_s
integer , DIMENSION(max_domains) :: inputout_interval_d
integer , DIMENSION(max_domains) :: inputout_interval_h
integer , DIMENSION(max_domains) :: inputout_interval_m
integer , DIMENSION(max_domains) :: inputout_interval_s
integer , DIMENSION(max_domains) :: inputout_interval
integer :: restart_interval_d
integer :: restart_interval_h
integer :: restart_interval_m
integer :: restart_interval_s
integer , DIMENSION(max_domains) :: history_begin_y
integer , DIMENSION(max_domains) :: history_begin_d
integer , DIMENSION(max_domains) :: history_begin_h
integer , DIMENSION(max_domains) :: history_begin_m
integer , DIMENSION(max_domains) :: history_begin_s
integer , DIMENSION(max_domains) :: history_begin
integer , DIMENSION(max_domains) :: inputout_begin_y
integer , DIMENSION(max_domains) :: inputout_begin_d
integer , DIMENSION(max_domains) :: inputout_begin_h
integer , DIMENSION(max_domains) :: inputout_begin_m
integer , DIMENSION(max_domains) :: inputout_begin_s
integer :: restart_begin_y
integer :: restart_begin_d
integer :: restart_begin_h
integer :: restart_begin_m
integer :: restart_begin_s
integer :: restart_begin
integer , DIMENSION(max_domains) :: history_end_y
integer , DIMENSION(max_domains) :: history_end_d
integer , DIMENSION(max_domains) :: history_end_h
integer , DIMENSION(max_domains) :: history_end_m
integer , DIMENSION(max_domains) :: history_end_s
integer , DIMENSION(max_domains) :: history_end
integer , DIMENSION(max_domains) :: inputout_end_y
integer , DIMENSION(max_domains) :: inputout_end_d
integer , DIMENSION(max_domains) :: inputout_end_h
integer , DIMENSION(max_domains) :: inputout_end_m
integer , DIMENSION(max_domains) :: inputout_end_s
integer :: simulation_start_year
integer :: simulation_start_month
integer :: simulation_start_day
integer :: simulation_start_hour
integer :: simulation_start_minute
integer :: simulation_start_second
logical :: reset_simulation_start
integer , DIMENSION(max_domains) :: sr_x
integer , DIMENSION(max_domains) :: sr_y
character*256 :: sgfdda_inname
character*256 :: gfdda_inname
integer , DIMENSION(max_domains) :: sgfdda_interval_d
integer , DIMENSION(max_domains) :: sgfdda_interval_h
integer , DIMENSION(max_domains) :: sgfdda_interval_m
integer , DIMENSION(max_domains) :: sgfdda_interval_s
integer , DIMENSION(max_domains) :: sgfdda_interval_y
integer , DIMENSION(max_domains) :: sgfdda_interval
integer , DIMENSION(max_domains) :: gfdda_interval_d
integer , DIMENSION(max_domains) :: gfdda_interval_h
integer , DIMENSION(max_domains) :: gfdda_interval_m
integer , DIMENSION(max_domains) :: gfdda_interval_s
integer , DIMENSION(max_domains) :: gfdda_interval_y
integer , DIMENSION(max_domains) :: gfdda_interval
integer , DIMENSION(max_domains) :: sgfdda_begin_y
integer , DIMENSION(max_domains) :: sgfdda_begin_d
integer , DIMENSION(max_domains) :: sgfdda_begin_h
integer , DIMENSION(max_domains) :: sgfdda_begin_m
integer , DIMENSION(max_domains) :: sgfdda_begin_s
integer , DIMENSION(max_domains) :: gfdda_begin_y
integer , DIMENSION(max_domains) :: gfdda_begin_d
integer , DIMENSION(max_domains) :: gfdda_begin_h
integer , DIMENSION(max_domains) :: gfdda_begin_m
integer , DIMENSION(max_domains) :: gfdda_begin_s
integer , DIMENSION(max_domains) :: sgfdda_end_y
integer , DIMENSION(max_domains) :: sgfdda_end_d
integer , DIMENSION(max_domains) :: sgfdda_end_h
integer , DIMENSION(max_domains) :: sgfdda_end_m
integer , DIMENSION(max_domains) :: sgfdda_end_s
integer , DIMENSION(max_domains) :: gfdda_end_y
integer , DIMENSION(max_domains) :: gfdda_end_d
integer , DIMENSION(max_domains) :: gfdda_end_h
integer , DIMENSION(max_domains) :: gfdda_end_m
integer , DIMENSION(max_domains) :: gfdda_end_s
integer :: io_form_sgfdda
integer :: io_form_gfdda
character*256 , DIMENSION(max_domains) :: iofields_filename
logical :: ignore_iofields_warning
logical :: ncd_nofill
logical :: update_sfcdiags
logical :: use_wrf_sfcinfo
logical :: use_background_errors
logical :: write_increments
logical :: var4d
integer :: var4d_bin
integer :: var4d_bin_rain
logical :: var4d_lbc
integer :: multi_inc
logical :: print_detail_radar
logical :: print_detail_rain
logical :: print_detail_rad
logical :: print_detail_xa
logical :: print_detail_xb
logical :: print_detail_obs
logical :: print_detail_f_obs
logical :: print_detail_map
logical :: print_detail_grad
logical :: print_detail_regression
logical :: print_detail_spectral
logical :: print_detail_testing
logical :: print_detail_parallel
logical :: print_detail_be
logical :: print_detail_outerloop
logical :: check_max_iv_print
logical :: check_buddy_print
integer :: analysis_accu
logical :: calc_w_increment
logical :: dt_cloud_model
logical :: write_mod_filtered_obs
logical :: wind_sd
logical :: wind_sd_buoy
logical :: wind_sd_synop
logical :: wind_sd_ships
logical :: wind_sd_metar
logical :: wind_sd_sound
logical :: wind_sd_pilot
logical :: wind_sd_airep
logical :: wind_sd_qscat
logical :: wind_sd_tamdar
logical :: wind_sd_geoamv
logical :: wind_sd_mtgirs
logical :: wind_sd_polaramv
logical :: wind_sd_profiler
logical :: wind_stats_sd
logical :: qc_rej_both
integer :: fg_format
integer :: ob_format
integer :: ob_format_gpsro
integer :: num_fgat_time
logical :: thin_conv
logical :: thin_conv_ascii
real , DIMENSION(num_ob_indexes) :: thin_mesh_conv
logical :: thin_rainobs
logical :: use_synopobs
logical :: use_shipsobs
logical :: use_metarobs
logical :: use_soundobs
logical :: use_mtgirsobs
logical :: use_tamdarobs
logical :: use_pilotobs
logical :: use_airepobs
logical :: use_geoamvobs
logical :: use_polaramvobs
logical :: use_bogusobs
logical :: use_buoyobs
logical :: use_profilerobs
logical :: use_satemobs
logical :: use_gpsztdobs
logical :: use_gpspwobs
logical :: use_gpsrefobs
real :: top_km_gpsro
real :: bot_km_gpsro
logical :: use_ssmiretrievalobs
logical :: use_ssmitbobs
logical :: use_ssmt1obs
logical :: use_ssmt2obs
logical :: use_qscatobs
logical :: use_radarobs
logical :: use_radar_rv
logical :: use_radar_rf
logical :: use_radar_rqv
logical :: use_radar_rhv
logical :: use_3dvar_phy
logical :: use_rainobs
logical :: use_hirs2obs
logical :: use_hirs3obs
logical :: use_hirs4obs
logical :: use_mhsobs
logical :: use_msuobs
logical :: use_amsuaobs
logical :: use_amsubobs
logical :: use_airsobs
logical :: use_airsretobs
logical :: use_eos_amsuaobs
logical :: use_hsbobs
logical :: use_ssmisobs
logical :: use_iasiobs
logical :: use_seviriobs
logical :: use_amsr2obs
logical :: use_kma1dvar
logical :: use_filtered_rad
logical :: use_obs_errfac
logical :: use_atmsobs
logical :: use_mwtsobs
logical :: use_mwhsobs
logical :: check_max_iv
real :: max_error_t
real :: max_error_uv
real :: max_error_spd
real :: max_error_dir
real :: max_omb_spd
real :: max_omb_dir
real :: max_error_pw
real :: max_error_ref
real :: max_error_rh
real :: max_error_q
real :: max_error_p
real :: max_error_tb
real :: max_error_thickness
real :: max_error_rv
real :: max_error_rf
real :: max_error_rain
real :: max_error_buv
real :: max_error_bt
real :: max_error_bq
real :: max_error_slp
logical :: check_buddy
logical :: put_rand_seed
logical :: omb_set_rand
logical :: omb_add_noise
logical :: position_lev_dependant
integer :: obs_qc_pointer
integer :: qmarker_retain
integer :: max_sound_input
integer :: max_mtgirs_input
integer :: max_tamdar_input
integer :: max_synop_input
integer :: max_geoamv_input
integer :: max_polaramv_input
integer :: max_airep_input
integer :: max_satem_input
integer :: max_pilot_input
integer :: max_radar_input
integer :: max_rain_input
integer :: max_metar_input
integer :: max_gpspw_input
integer :: max_ships_input
integer :: max_profiler_input
integer :: max_bogus_input
integer :: max_buoy_input
integer :: max_ssmi_rv_input
integer :: max_ssmi_tb_input
integer :: max_ssmt1_input
integer :: max_ssmt2_input
integer :: max_qscat_input
integer :: max_gpsref_input
integer :: max_airsr_input
integer :: max_tovs_input
integer :: max_ssmis_input
integer :: report_start
integer :: report_end
integer :: tovs_start
integer :: tovs_end
logical :: gpsref_thinning
logical :: outer_loop_restart
integer :: max_ext_its
integer , DIMENSION(max_outer_iterations) :: ntmax
integer :: nsave
integer :: write_interval
real , DIMENSION(max_outer_iterations) :: eps
logical :: precondition_cg
real :: precondition_factor
logical :: use_lanczos
logical :: read_lanczos
logical :: write_lanczos
logical :: orthonorm_gradient
integer :: cv_options
integer :: cloud_cv_options
real , DIMENSION(3*max_outer_iterations) :: as1
real , DIMENSION(3*max_outer_iterations) :: as2
real , DIMENSION(3*max_outer_iterations) :: as3
real , DIMENSION(3*max_outer_iterations) :: as4
real , DIMENSION(3*max_outer_iterations) :: as5
logical :: do_normalize
logical :: use_rf
integer :: rf_passes
real , DIMENSION(max_outer_iterations) :: var_scaling1
real , DIMENSION(max_outer_iterations) :: var_scaling2
real , DIMENSION(max_outer_iterations) :: var_scaling3
real , DIMENSION(max_outer_iterations) :: var_scaling4
real , DIMENSION(max_outer_iterations) :: var_scaling5
real , DIMENSION(max_outer_iterations) :: var_scaling6
real , DIMENSION(max_outer_iterations) :: var_scaling7
real , DIMENSION(max_outer_iterations) :: var_scaling8
real , DIMENSION(max_outer_iterations) :: var_scaling9
real , DIMENSION(max_outer_iterations) :: var_scaling10
real , DIMENSION(max_outer_iterations) :: var_scaling11
real , DIMENSION(max_outer_iterations) :: len_scaling1
real , DIMENSION(max_outer_iterations) :: len_scaling2
real , DIMENSION(max_outer_iterations) :: len_scaling3
real , DIMENSION(max_outer_iterations) :: len_scaling4
real , DIMENSION(max_outer_iterations) :: len_scaling5
real , DIMENSION(max_outer_iterations) :: len_scaling6
real , DIMENSION(max_outer_iterations) :: len_scaling7
real , DIMENSION(max_outer_iterations) :: len_scaling8
real , DIMENSION(max_outer_iterations) :: len_scaling9
real , DIMENSION(max_outer_iterations) :: len_scaling10
real , DIMENSION(max_outer_iterations) :: len_scaling11
real :: je_factor
real :: power_truncation
logical :: def_sub_domain
real :: x_start_sub_domain
real :: y_start_sub_domain
real :: x_end_sub_domain
real :: y_end_sub_domain
integer :: stdout
integer :: stderr
integer :: trace_unit
integer :: trace_pe
integer :: trace_repeat_head
integer :: trace_repeat_body
integer :: trace_max_depth
logical :: trace_use
logical :: trace_use_frequent
logical :: trace_use_dull
logical :: trace_memory
logical :: trace_all_pes
logical :: trace_csv
logical :: use_html
logical :: warnings_are_fatal
logical :: test_transforms
logical :: test_gradient
logical :: test_statistics
logical :: interpolate_stats
real , DIMENSION(99) :: be_eta
logical :: test_dm_exact
integer :: cv_options_hum
integer :: check_rh
real :: set_omb_rand_fac
integer :: seed_array1
integer :: seed_array2
integer :: sfc_assi_options
logical :: psfc_from_slp
logical :: calculate_cg_cost_fn
logical :: lat_stats_option
integer :: interp_option
integer :: balance_type
logical :: use_wpec
real :: wpec_factor
integer :: vert_corr
integer :: vertical_ip
integer :: vert_evalue
real :: max_vert_var1
real :: max_vert_var2
real :: max_vert_var3
real :: max_vert_var4
real :: max_vert_var5
real :: max_vert_var6
real :: max_vert_var7
real :: max_vert_var8
real :: max_vert_var9
real :: max_vert_var10
real :: max_vert_var11
real :: max_vert_var_alpha
real :: psi_chi_factor
real :: psi_t_factor
real :: psi_ps_factor
real :: psi_rh_factor
real :: chi_u_t_factor
real :: chi_u_ps_factor
real :: chi_u_rh_factor
real :: t_u_rh_factor
real :: ps_u_rh_factor
integer :: rttov_emis_atlas_ir
integer :: rttov_emis_atlas_mw
integer :: rtminit_print
integer :: rtminit_nsensor
integer , DIMENSION(max_instruments) :: rtminit_platform
integer , DIMENSION(max_instruments) :: rtminit_satid
integer , DIMENSION(max_instruments) :: rtminit_sensor
integer , DIMENSION(max_instruments) :: rad_monitoring
real , DIMENSION(max_instruments) :: thinning_mesh
logical :: thinning
logical :: read_biascoef
logical :: biascorr
logical :: biasprep
logical :: rttov_scatt
logical :: write_profile
logical :: write_jacobian
logical :: qc_rad
logical :: write_iv_rad_ascii
logical :: write_oa_rad_ascii
logical :: write_filtered_rad
logical :: use_error_factor_rad
logical :: use_landem
logical , DIMENSION(max_instruments) :: use_antcorr
logical , DIMENSION(max_instruments) :: use_mspps_emis
logical , DIMENSION(max_instruments) :: use_mspps_ts
integer :: mw_emis_sea
integer :: tovs_min_transfer
logical :: tovs_batch
integer :: rtm_option
logical :: use_crtm_kmatrix
logical :: use_rttov_kmatrix
logical :: crtm_cloud
logical :: only_sea_rad
logical :: use_pseudo_rad
integer :: pseudo_rad_platid
integer :: pseudo_rad_satid
integer :: pseudo_rad_senid
integer :: pseudo_rad_ichan
real :: pseudo_rad_lat
real :: pseudo_rad_lon
real :: pseudo_rad_inv
real :: pseudo_rad_err
logical :: use_simulated_rad
logical :: simulated_rad_io
integer :: simulated_rad_ngrid
logical :: use_varbc
logical :: freeze_varbc
real :: varbc_factor
integer :: varbc_nbgerr
integer :: varbc_nobsmin
logical :: use_clddet_mmr
logical :: use_clddet_ecmwf
logical :: airs_warmest_fov
logical , DIMENSION(2) :: use_satcv
logical :: use_blacklist_rad
logical :: calc_weightfunc
character*256 :: crtm_coef_path
character*256 :: crtm_irwater_coef
character*256 :: crtm_mwwater_coef
character*256 :: crtm_irland_coef
character*256 :: crtm_visland_coef
integer :: num_pseudo
real :: pseudo_x
real :: pseudo_y
real :: pseudo_z
real :: pseudo_val
real :: pseudo_err
integer :: alphacv_method
integer :: ensdim_alpha
integer :: alpha_truncation
integer :: alpha_corr_type
real :: alpha_corr_scale
real :: alpha_std_dev
logical :: alpha_vertloc
logical :: alpha_hydrometeors
logical :: hybrid_dual_res
integer :: dual_res_upscale_opt
character*256 :: analysis_type
integer :: sensitivity_option
logical :: adj_sens
character*256 :: analysis_date
character*256 :: pseudo_var
character*256 :: documentation_url
character*256 :: time_window_min
character*256 :: time_window_max
logical :: jcdfi_use
integer :: jcdfi_diag
real :: jcdfi_penalty
logical :: enable_identity
logical :: trajectory_io
logical :: var4d_detail_out
logical :: var4d_run
integer , DIMENSION(max_domains) :: mp_physics_ad
integer , DIMENSION(max_domains) :: mp_physics_4dvar
integer , DIMENSION(max_domains) :: chem_opt
integer    :: last_item_in_struct









NAMELIST /time_control/ run_days
NAMELIST /time_control/ run_hours
NAMELIST /time_control/ run_minutes
NAMELIST /time_control/ run_seconds
NAMELIST /time_control/ start_year
NAMELIST /time_control/ start_month
NAMELIST /time_control/ start_day
NAMELIST /time_control/ start_hour
NAMELIST /time_control/ start_minute
NAMELIST /time_control/ start_second
NAMELIST /time_control/ end_year
NAMELIST /time_control/ end_month
NAMELIST /time_control/ end_day
NAMELIST /time_control/ end_hour
NAMELIST /time_control/ end_minute
NAMELIST /time_control/ end_second
NAMELIST /time_control/ interval_seconds
NAMELIST /time_control/ input_from_file
NAMELIST /time_control/ fine_input_stream
NAMELIST /time_control/ input_from_hires
NAMELIST /time_control/ rsmas_data_path
NAMELIST /time_control/ all_ic_times
NAMELIST /time_control/ julyr
NAMELIST /time_control/ julday
NAMELIST /time_control/ gmt
NAMELIST /time_control/ input_inname
NAMELIST /time_control/ input_outname
NAMELIST /time_control/ bdy_inname
NAMELIST /time_control/ bdy_outname
NAMELIST /time_control/ rst_inname
NAMELIST /time_control/ rst_outname
NAMELIST /time_control/ write_input
NAMELIST /time_control/ write_restart_at_0h
NAMELIST /time_control/ write_hist_at_0h_rst
NAMELIST /time_control/ adjust_output_times
NAMELIST /time_control/ adjust_input_times
NAMELIST /time_control/ diag_print
NAMELIST /time_control/ nocolons
NAMELIST /time_control/ cycling
NAMELIST /time_control/ output_diagnostics
NAMELIST /time_control/ nwp_diagnostics
NAMELIST /time_control/ output_ready_flag
NAMELIST /pio_control/ usepio
NAMELIST /pio_control/ pioprocs
NAMELIST /pio_control/ piostart
NAMELIST /pio_control/ piostride
NAMELIST /pio_control/ pioshift
NAMELIST /dfi_control/ dfi_opt
NAMELIST /dfi_control/ dfi_savehydmeteors
NAMELIST /dfi_control/ dfi_nfilter
NAMELIST /dfi_control/ dfi_write_filtered_input
NAMELIST /dfi_control/ dfi_write_dfi_history
NAMELIST /dfi_control/ dfi_cutoff_seconds
NAMELIST /dfi_control/ dfi_time_dim
NAMELIST /dfi_control/ dfi_fwdstop_year
NAMELIST /dfi_control/ dfi_fwdstop_month
NAMELIST /dfi_control/ dfi_fwdstop_day
NAMELIST /dfi_control/ dfi_fwdstop_hour
NAMELIST /dfi_control/ dfi_fwdstop_minute
NAMELIST /dfi_control/ dfi_fwdstop_second
NAMELIST /dfi_control/ dfi_bckstop_year
NAMELIST /dfi_control/ dfi_bckstop_month
NAMELIST /dfi_control/ dfi_bckstop_day
NAMELIST /dfi_control/ dfi_bckstop_hour
NAMELIST /dfi_control/ dfi_bckstop_minute
NAMELIST /dfi_control/ dfi_bckstop_second
NAMELIST /domains/ time_step
NAMELIST /domains/ time_step_fract_num
NAMELIST /domains/ time_step_fract_den
NAMELIST /domains/ time_step_dfi
NAMELIST /domains/ min_time_step
NAMELIST /domains/ min_time_step_den
NAMELIST /domains/ max_time_step
NAMELIST /domains/ max_time_step_den
NAMELIST /domains/ target_cfl
NAMELIST /domains/ target_hcfl
NAMELIST /domains/ max_step_increase_pct
NAMELIST /domains/ starting_time_step
NAMELIST /domains/ starting_time_step_den
NAMELIST /domains/ step_to_output_time
NAMELIST /domains/ adaptation_domain
NAMELIST /domains/ use_adaptive_time_step
NAMELIST /domains/ use_adaptive_time_step_dfi
NAMELIST /domains/ max_dom
NAMELIST /domains/ lats_to_mic
NAMELIST /domains/ s_we
NAMELIST /domains/ e_we
NAMELIST /domains/ s_sn
NAMELIST /domains/ e_sn
NAMELIST /domains/ s_vert
NAMELIST /domains/ e_vert
NAMELIST /domains/ num_metgrid_levels
NAMELIST /domains/ num_metgrid_soil_levels
NAMELIST /domains/ p_top_requested
NAMELIST /domains/ interp_theta
NAMELIST /domains/ interp_type
NAMELIST /domains/ rebalance
NAMELIST /domains/ vert_refine_method
NAMELIST /domains/ vert_refine_fact
NAMELIST /domains/ extrap_type
NAMELIST /domains/ t_extrap_type
NAMELIST /domains/ hypsometric_opt
NAMELIST /domains/ lowest_lev_from_sfc
NAMELIST /domains/ use_levels_below_ground
NAMELIST /domains/ use_tavg_for_tsk
NAMELIST /domains/ use_surface
NAMELIST /domains/ lagrange_order
NAMELIST /domains/ force_sfc_in_vinterp
NAMELIST /domains/ zap_close_levels
NAMELIST /domains/ maxw_horiz_pres_diff
NAMELIST /domains/ trop_horiz_pres_diff
NAMELIST /domains/ maxw_above_this_level
NAMELIST /domains/ use_maxw_level
NAMELIST /domains/ use_trop_level
NAMELIST /domains/ sfcp_to_sfcp
NAMELIST /domains/ adjust_heights
NAMELIST /domains/ smooth_cg_topo
NAMELIST /domains/ nest_interp_coord
NAMELIST /domains/ interp_method_type
NAMELIST /domains/ aggregate_lu
NAMELIST /domains/ rh2qv_wrt_liquid
NAMELIST /domains/ rh2qv_method
NAMELIST /domains/ qv_max_p_safe
NAMELIST /domains/ qv_max_flag
NAMELIST /domains/ qv_max_value
NAMELIST /domains/ qv_min_p_safe
NAMELIST /domains/ qv_min_flag
NAMELIST /domains/ qv_min_value
NAMELIST /domains/ ideal_init_method
NAMELIST /domains/ dx
NAMELIST /domains/ dy
NAMELIST /domains/ grid_id
NAMELIST /domains/ grid_allowed
NAMELIST /domains/ parent_id
NAMELIST /domains/ i_parent_start
NAMELIST /domains/ j_parent_start
NAMELIST /domains/ parent_grid_ratio
NAMELIST /domains/ parent_time_step_ratio
NAMELIST /domains/ feedback
NAMELIST /domains/ smooth_option
NAMELIST /domains/ blend_width
NAMELIST /domains/ ztop
NAMELIST /domains/ moad_grid_ratio
NAMELIST /domains/ moad_time_step_ratio
NAMELIST /domains/ shw
NAMELIST /domains/ tile_sz_x
NAMELIST /domains/ tile_sz_y
NAMELIST /domains/ numtiles
NAMELIST /domains/ numtiles_inc
NAMELIST /domains/ numtiles_x
NAMELIST /domains/ numtiles_y
NAMELIST /domains/ tile_strategy
NAMELIST /domains/ nproc_x
NAMELIST /domains/ nproc_y
NAMELIST /domains/ irand
NAMELIST /domains/ num_moves
NAMELIST /domains/ ts_buf_size
NAMELIST /domains/ max_ts_locs
NAMELIST /domains/ vortex_interval
NAMELIST /domains/ max_vortex_speed
NAMELIST /domains/ corral_dist
NAMELIST /domains/ track_level
NAMELIST /domains/ time_to_move
NAMELIST /domains/ move_id
NAMELIST /domains/ move_interval
NAMELIST /domains/ move_cd_x
NAMELIST /domains/ move_cd_y
NAMELIST /domains/ swap_x
NAMELIST /domains/ swap_y
NAMELIST /domains/ cycle_x
NAMELIST /domains/ cycle_y
NAMELIST /domains/ reorder_mesh
NAMELIST /domains/ perturb_input
NAMELIST /domains/ eta_levels
NAMELIST /domains/ max_dz
NAMELIST /domains/ ocean_levels
NAMELIST /domains/ ocean_z
NAMELIST /domains/ ocean_t
NAMELIST /domains/ ocean_s
NAMELIST /domains/ num_traj
NAMELIST /domains/ max_ts_level
NAMELIST /domains/ track_loc_in
NAMELIST /domains/ num_ext_model_couple_dom
NAMELIST /tc/ insert_bogus_storm
NAMELIST /tc/ remove_storm
NAMELIST /tc/ num_storm
NAMELIST /tc/ latc_loc
NAMELIST /tc/ lonc_loc
NAMELIST /tc/ vmax_meters_per_second
NAMELIST /tc/ rmax
NAMELIST /tc/ vmax_ratio
NAMELIST /tc/ rankine_lid
NAMELIST /physics/ force_read_thompson
NAMELIST /physics/ write_thompson_tables
NAMELIST /physics/ mp_physics
NAMELIST /physics/ nssl_cccn
NAMELIST /physics/ nssl_alphah
NAMELIST /physics/ nssl_alphahl
NAMELIST /physics/ nssl_cnoh
NAMELIST /physics/ nssl_cnohl
NAMELIST /physics/ nssl_cnor
NAMELIST /physics/ nssl_cnos
NAMELIST /physics/ nssl_rho_qh
NAMELIST /physics/ nssl_rho_qhl
NAMELIST /physics/ nssl_rho_qs
NAMELIST /physics/ nudge_lightning
NAMELIST /physics/ nudge_light_times
NAMELIST /physics/ nudge_light_timee
NAMELIST /physics/ nudge_light_int
NAMELIST /physics/ path_to_files
NAMELIST /physics/ gsfcgce_hail
NAMELIST /physics/ gsfcgce_2ice
NAMELIST /physics/ progn
NAMELIST /physics/ accum_mode
NAMELIST /physics/ aitken_mode
NAMELIST /physics/ coarse_mode
NAMELIST /physics/ do_radar_ref
NAMELIST /physics/ ra_lw_physics
NAMELIST /physics/ ra_sw_physics
NAMELIST /physics/ radt
NAMELIST /physics/ naer
NAMELIST /physics/ sf_sfclay_physics
NAMELIST /physics/ sf_surface_physics
NAMELIST /physics/ bl_pbl_physics
NAMELIST /physics/ bl_mynn_tkebudget
NAMELIST /physics/ ysu_topdown_pblmix
NAMELIST /physics/ shinhong_tke_diag
NAMELIST /physics/ bl_mynn_tkeadvect
NAMELIST /physics/ bl_mynn_cloudpdf
NAMELIST /physics/ bl_mynn_mixlength
NAMELIST /physics/ bl_mynn_edmf
NAMELIST /physics/ bl_mynn_edmf_mom
NAMELIST /physics/ bl_mynn_edmf_tke
NAMELIST /physics/ bl_mynn_edmf_part
NAMELIST /physics/ bl_mynn_cloudmix
NAMELIST /physics/ bl_mynn_mixqt
NAMELIST /physics/ icloud_bl
NAMELIST /physics/ mfshconv
NAMELIST /physics/ sf_urban_physics
NAMELIST /physics/ bldt
NAMELIST /physics/ cu_physics
NAMELIST /physics/ shcu_physics
NAMELIST /physics/ cu_diag
NAMELIST /physics/ kf_edrates
NAMELIST /physics/ kfeta_trigger
NAMELIST /physics/ nsas_dx_factor
NAMELIST /physics/ cudt
NAMELIST /physics/ gsmdt
NAMELIST /physics/ isfflx
NAMELIST /physics/ ifsnow
NAMELIST /physics/ icloud
NAMELIST /physics/ ideal_xland
NAMELIST /physics/ swrad_scat
NAMELIST /physics/ surface_input_source
NAMELIST /physics/ num_soil_layers
NAMELIST /physics/ maxpatch
NAMELIST /physics/ num_snow_layers
NAMELIST /physics/ num_snso_layers
NAMELIST /physics/ num_urban_layers
NAMELIST /physics/ num_urban_hi
NAMELIST /physics/ num_months
NAMELIST /physics/ sf_surface_mosaic
NAMELIST /physics/ mosaic_cat
NAMELIST /physics/ mosaic_lu
NAMELIST /physics/ mosaic_soil
NAMELIST /physics/ maxiens
NAMELIST /physics/ maxens
NAMELIST /physics/ maxens2
NAMELIST /physics/ maxens3
NAMELIST /physics/ ensdim
NAMELIST /physics/ cugd_avedx
NAMELIST /physics/ clos_choice
NAMELIST /physics/ imomentum
NAMELIST /physics/ ishallow
NAMELIST /physics/ convtrans_avglen_m
NAMELIST /physics/ num_land_cat
NAMELIST /physics/ num_soil_cat
NAMELIST /physics/ mp_zero_out
NAMELIST /physics/ mp_zero_out_thresh
NAMELIST /physics/ seaice_threshold
NAMELIST /physics/ sst_update
NAMELIST /physics/ sst_skin
NAMELIST /physics/ tmn_update
NAMELIST /physics/ usemonalb
NAMELIST /physics/ rdmaxalb
NAMELIST /physics/ rdlai2d
NAMELIST /physics/ ua_phys
NAMELIST /physics/ opt_thcnd
NAMELIST /physics/ co2tf
NAMELIST /physics/ ra_call_offset
NAMELIST /physics/ cam_abs_freq_s
NAMELIST /physics/ levsiz
NAMELIST /physics/ paerlev
NAMELIST /physics/ cam_abs_dim1
NAMELIST /physics/ cam_abs_dim2
NAMELIST /physics/ lagday
NAMELIST /physics/ no_src_types
NAMELIST /physics/ alevsiz
NAMELIST /physics/ o3input
NAMELIST /physics/ aer_opt
NAMELIST /physics/ swint_opt
NAMELIST /physics/ aer_type
NAMELIST /physics/ aer_aod550_opt
NAMELIST /physics/ aer_angexp_opt
NAMELIST /physics/ aer_ssa_opt
NAMELIST /physics/ aer_asy_opt
NAMELIST /physics/ aer_aod550_val
NAMELIST /physics/ aer_angexp_val
NAMELIST /physics/ aer_ssa_val
NAMELIST /physics/ aer_asy_val
NAMELIST /physics/ cu_rad_feedback
NAMELIST /physics/ shallowcu_forced_ra
NAMELIST /physics/ numbins
NAMELIST /physics/ thbinsize
NAMELIST /physics/ rbinsize
NAMELIST /physics/ mindeepfreq
NAMELIST /physics/ minshallowfreq
NAMELIST /physics/ shcu_aerosols_opt
NAMELIST /physics/ pxlsm_smois_init
NAMELIST /physics/ omlcall
NAMELIST /physics/ sf_ocean_physics
NAMELIST /physics/ traj_opt
NAMELIST /physics/ tracercall
NAMELIST /physics/ omdt
NAMELIST /physics/ oml_hml0
NAMELIST /physics/ oml_gamma
NAMELIST /physics/ oml_relaxation_time
NAMELIST /physics/ isftcflx
NAMELIST /physics/ iz0tlnd
NAMELIST /physics/ shadlen
NAMELIST /physics/ slope_rad
NAMELIST /physics/ topo_shading
NAMELIST /physics/ topo_wind
NAMELIST /physics/ no_mp_heating
NAMELIST /physics/ fractional_seaice
NAMELIST /physics/ seaice_snowdepth_opt
NAMELIST /physics/ seaice_snowdepth_max
NAMELIST /physics/ seaice_snowdepth_min
NAMELIST /physics/ seaice_albedo_opt
NAMELIST /physics/ seaice_albedo_default
NAMELIST /physics/ seaice_thickness_opt
NAMELIST /physics/ seaice_thickness_default
NAMELIST /physics/ tice2tsk_if2cold
NAMELIST /physics/ bucket_mm
NAMELIST /physics/ bucket_j
NAMELIST /physics/ mp_tend_lim
NAMELIST /physics/ prec_acc_dt
NAMELIST /physics/ grav_settling
NAMELIST /physics/ sas_pgcon
NAMELIST /physics/ scalar_pblmix
NAMELIST /physics/ tracer_pblmix
NAMELIST /physics/ use_aero_icbc
NAMELIST /physics/ use_rap_aero_icbc
NAMELIST /physics/ use_mp_re
NAMELIST /physics/ ccn_conc
NAMELIST /physics/ hail_opt
NAMELIST /noah_mp/ dveg
NAMELIST /noah_mp/ opt_crs
NAMELIST /noah_mp/ opt_btr
NAMELIST /noah_mp/ opt_run
NAMELIST /noah_mp/ opt_sfc
NAMELIST /noah_mp/ opt_frz
NAMELIST /noah_mp/ opt_inf
NAMELIST /noah_mp/ opt_rad
NAMELIST /noah_mp/ opt_alb
NAMELIST /noah_mp/ opt_snf
NAMELIST /noah_mp/ opt_tbot
NAMELIST /noah_mp/ opt_stc
NAMELIST /noah_mp/ opt_gla
NAMELIST /noah_mp/ opt_rsf
NAMELIST /physics/ wtddt
NAMELIST /fdda/ fgdt
NAMELIST /fdda/ fgdtzero
NAMELIST /fdda/ grid_fdda
NAMELIST /fdda/ grid_sfdda
NAMELIST /fdda/ if_no_pbl_nudging_uv
NAMELIST /fdda/ if_no_pbl_nudging_t
NAMELIST /fdda/ if_no_pbl_nudging_ph
NAMELIST /fdda/ if_no_pbl_nudging_q
NAMELIST /fdda/ if_zfac_uv
NAMELIST /fdda/ k_zfac_uv
NAMELIST /fdda/ if_zfac_t
NAMELIST /fdda/ k_zfac_t
NAMELIST /fdda/ if_zfac_ph
NAMELIST /fdda/ k_zfac_ph
NAMELIST /fdda/ if_zfac_q
NAMELIST /fdda/ k_zfac_q
NAMELIST /fdda/ dk_zfac_uv
NAMELIST /fdda/ dk_zfac_t
NAMELIST /fdda/ dk_zfac_ph
NAMELIST /fdda/ guv
NAMELIST /fdda/ guv_sfc
NAMELIST /fdda/ gt
NAMELIST /fdda/ gt_sfc
NAMELIST /fdda/ gq
NAMELIST /fdda/ gq_sfc
NAMELIST /fdda/ gph
NAMELIST /fdda/ dtramp_min
NAMELIST /fdda/ if_ramping
NAMELIST /fdda/ rinblw
NAMELIST /fdda/ xwavenum
NAMELIST /fdda/ ywavenum
NAMELIST /fdda/ pxlsm_soil_nudge
NAMELIST /fdda/ obs_nudge_opt
NAMELIST /fdda/ max_obs
NAMELIST /fdda/ fdda_start
NAMELIST /fdda/ fdda_end
NAMELIST /fdda/ obs_nudge_wind
NAMELIST /fdda/ obs_coef_wind
NAMELIST /fdda/ obs_nudge_temp
NAMELIST /fdda/ obs_coef_temp
NAMELIST /fdda/ obs_nudge_mois
NAMELIST /fdda/ obs_coef_mois
NAMELIST /fdda/ obs_nudge_pstr
NAMELIST /fdda/ obs_coef_pstr
NAMELIST /fdda/ obs_no_pbl_nudge_uv
NAMELIST /fdda/ obs_no_pbl_nudge_t
NAMELIST /fdda/ obs_no_pbl_nudge_q
NAMELIST /fdda/ obs_sfc_scheme_horiz
NAMELIST /fdda/ obs_sfc_scheme_vert
NAMELIST /fdda/ obs_max_sndng_gap
NAMELIST /fdda/ obs_nudgezfullr1_uv
NAMELIST /fdda/ obs_nudgezrampr1_uv
NAMELIST /fdda/ obs_nudgezfullr2_uv
NAMELIST /fdda/ obs_nudgezrampr2_uv
NAMELIST /fdda/ obs_nudgezfullr4_uv
NAMELIST /fdda/ obs_nudgezrampr4_uv
NAMELIST /fdda/ obs_nudgezfullr1_t
NAMELIST /fdda/ obs_nudgezrampr1_t
NAMELIST /fdda/ obs_nudgezfullr2_t
NAMELIST /fdda/ obs_nudgezrampr2_t
NAMELIST /fdda/ obs_nudgezfullr4_t
NAMELIST /fdda/ obs_nudgezrampr4_t
NAMELIST /fdda/ obs_nudgezfullr1_q
NAMELIST /fdda/ obs_nudgezrampr1_q
NAMELIST /fdda/ obs_nudgezfullr2_q
NAMELIST /fdda/ obs_nudgezrampr2_q
NAMELIST /fdda/ obs_nudgezfullr4_q
NAMELIST /fdda/ obs_nudgezrampr4_q
NAMELIST /fdda/ obs_nudgezfullmin
NAMELIST /fdda/ obs_nudgezrampmin
NAMELIST /fdda/ obs_nudgezmax
NAMELIST /fdda/ obs_sfcfact
NAMELIST /fdda/ obs_sfcfacr
NAMELIST /fdda/ obs_dpsmx
NAMELIST /fdda/ obs_rinxy
NAMELIST /fdda/ obs_rinsig
NAMELIST /fdda/ obs_twindo
NAMELIST /fdda/ obs_npfi
NAMELIST /fdda/ obs_ionf
NAMELIST /fdda/ obs_idynin
NAMELIST /fdda/ obs_dtramp
NAMELIST /fdda/ obs_prt_max
NAMELIST /fdda/ obs_prt_freq
NAMELIST /fdda/ obs_ipf_in4dob
NAMELIST /fdda/ obs_ipf_errob
NAMELIST /fdda/ obs_ipf_nudob
NAMELIST /fdda/ obs_ipf_init
NAMELIST /fdda/ obs_scl_neg_qv_innov
NAMELIST /scm/ scm_force
NAMELIST /scm/ scm_force_dx
NAMELIST /scm/ num_force_layers
NAMELIST /scm/ scm_lu_index
NAMELIST /scm/ scm_isltyp
NAMELIST /scm/ scm_vegfra
NAMELIST /scm/ scm_canwat
NAMELIST /scm/ scm_lat
NAMELIST /scm/ scm_lon
NAMELIST /scm/ scm_th_t_tend
NAMELIST /scm/ scm_qv_t_tend
NAMELIST /scm/ scm_th_adv
NAMELIST /scm/ scm_wind_adv
NAMELIST /scm/ scm_qv_adv
NAMELIST /scm/ scm_ql_adv
NAMELIST /scm/ scm_vert_adv
NAMELIST /scm/ num_force_soil_layers
NAMELIST /scm/ scm_soilt_force
NAMELIST /scm/ scm_soilq_force
NAMELIST /scm/ scm_force_th_largescale
NAMELIST /scm/ scm_force_qv_largescale
NAMELIST /scm/ scm_force_ql_largescale
NAMELIST /scm/ scm_force_wind_largescale
NAMELIST /scm/ scm_force_skintemp
NAMELIST /scm/ scm_force_flux
NAMELIST /dynamics/ dyn_opt
NAMELIST /dynamics/ rk_ord
NAMELIST /dynamics/ w_damping
NAMELIST /dynamics/ diff_opt
NAMELIST /dynamics/ diff_opt_dfi
NAMELIST /dynamics/ km_opt
NAMELIST /dynamics/ km_opt_dfi
NAMELIST /dynamics/ damp_opt
NAMELIST /dynamics/ rad_nudge
NAMELIST /dynamics/ gwd_opt
NAMELIST /dynamics/ zdamp
NAMELIST /dynamics/ dampcoef
NAMELIST /dynamics/ khdif
NAMELIST /dynamics/ kvdif
NAMELIST /dynamics/ diff_6th_factor
NAMELIST /dynamics/ diff_6th_opt
NAMELIST /dynamics/ use_theta_m
NAMELIST /dynamics/ use_q_diabatic
NAMELIST /dynamics/ c_s
NAMELIST /dynamics/ c_k
NAMELIST /dynamics/ smdiv
NAMELIST /dynamics/ emdiv
NAMELIST /dynamics/ epssm
NAMELIST /dynamics/ non_hydrostatic
NAMELIST /dynamics/ use_input_w
NAMELIST /dynamics/ time_step_sound
NAMELIST /dynamics/ h_mom_adv_order
NAMELIST /dynamics/ v_mom_adv_order
NAMELIST /dynamics/ h_sca_adv_order
NAMELIST /dynamics/ v_sca_adv_order
NAMELIST /dynamics/ momentum_adv_opt
NAMELIST /dynamics/ moist_adv_opt
NAMELIST /dynamics/ moist_adv_dfi_opt
NAMELIST /dynamics/ chem_adv_opt
NAMELIST /dynamics/ tracer_adv_opt
NAMELIST /dynamics/ scalar_adv_opt
NAMELIST /dynamics/ tke_adv_opt
NAMELIST /dynamics/ top_radiation
NAMELIST /dynamics/ mix_isotropic
NAMELIST /dynamics/ mix_upper_bound
NAMELIST /dynamics/ top_lid
NAMELIST /dynamics/ tke_upper_bound
NAMELIST /dynamics/ tke_drag_coefficient
NAMELIST /dynamics/ tke_heat_flux
NAMELIST /dynamics/ pert_coriolis
NAMELIST /dynamics/ coriolis2d
NAMELIST /dynamics/ mix_full_fields
NAMELIST /dynamics/ base_pres
NAMELIST /dynamics/ base_temp
NAMELIST /dynamics/ base_lapse
NAMELIST /dynamics/ iso_temp
NAMELIST /dynamics/ base_pres_strat
NAMELIST /dynamics/ base_lapse_strat
NAMELIST /dynamics/ use_baseparam_fr_nml
NAMELIST /dynamics/ fft_filter_lat
NAMELIST /dynamics/ coupled_filtering
NAMELIST /dynamics/ pos_def
NAMELIST /dynamics/ swap_pole_with_next_j
NAMELIST /dynamics/ actual_distance_average
NAMELIST /dynamics/ rotated_pole
NAMELIST /dynamics/ do_coriolis
NAMELIST /dynamics/ do_curvature
NAMELIST /dynamics/ do_gradp
NAMELIST /dynamics/ tracer_opt
NAMELIST /dynamics/ tenddiag
NAMELIST /bdy_control/ spec_bdy_width
NAMELIST /bdy_control/ spec_zone
NAMELIST /bdy_control/ relax_zone
NAMELIST /bdy_control/ specified
NAMELIST /bdy_control/ constant_bc
NAMELIST /bdy_control/ periodic_x
NAMELIST /bdy_control/ symmetric_xs
NAMELIST /bdy_control/ symmetric_xe
NAMELIST /bdy_control/ open_xs
NAMELIST /bdy_control/ open_xe
NAMELIST /bdy_control/ periodic_y
NAMELIST /bdy_control/ symmetric_ys
NAMELIST /bdy_control/ symmetric_ye
NAMELIST /bdy_control/ open_ys
NAMELIST /bdy_control/ open_ye
NAMELIST /bdy_control/ polar
NAMELIST /bdy_control/ nested
NAMELIST /bdy_control/ spec_exp
NAMELIST /bdy_control/ spec_bdy_final_mu
NAMELIST /bdy_control/ real_data_init_type
NAMELIST /bdy_control/ have_bcs_moist
NAMELIST /bdy_control/ have_bcs_scalar
NAMELIST /grib2/ background_proc_id
NAMELIST /grib2/ forecast_proc_id
NAMELIST /grib2/ production_status
NAMELIST /grib2/ compression
NAMELIST /physics/ windfarm_opt
NAMELIST /physics/ windfarm_ij
NAMELIST /physics/ lightning_option
NAMELIST /physics/ lightning_dt
NAMELIST /physics/ lightning_start_seconds
NAMELIST /physics/ flashrate_factor
NAMELIST /physics/ iccg_method
NAMELIST /physics/ iccg_prescribed_num
NAMELIST /physics/ iccg_prescribed_den
NAMELIST /physics/ cellcount_method
NAMELIST /physics/ cldtop_adjustment
NAMELIST /physics/ sf_lake_physics
NAMELIST /time_control/ auxinput1_inname
NAMELIST /time_control/ io_form_auxinput1
NAMELIST /time_control/ override_restart_timers
NAMELIST /time_control/ auxhist1_inname
NAMELIST /time_control/ auxhist1_outname
NAMELIST /time_control/ auxhist1_interval_y
NAMELIST /time_control/ auxhist1_interval_d
NAMELIST /time_control/ auxhist1_interval_h
NAMELIST /time_control/ auxhist1_interval_m
NAMELIST /time_control/ auxhist1_interval_s
NAMELIST /time_control/ auxhist1_interval
NAMELIST /time_control/ auxhist1_begin_y
NAMELIST /time_control/ auxhist1_begin_d
NAMELIST /time_control/ auxhist1_begin_h
NAMELIST /time_control/ auxhist1_begin_m
NAMELIST /time_control/ auxhist1_begin_s
NAMELIST /time_control/ auxhist1_begin
NAMELIST /time_control/ auxhist1_end_y
NAMELIST /time_control/ auxhist1_end_d
NAMELIST /time_control/ auxhist1_end_h
NAMELIST /time_control/ auxhist1_end_m
NAMELIST /time_control/ auxhist1_end_s
NAMELIST /time_control/ auxhist1_end
NAMELIST /time_control/ io_form_auxhist1
NAMELIST /time_control/ frames_per_auxhist1
NAMELIST /time_control/ auxhist2_inname
NAMELIST /time_control/ auxhist2_outname
NAMELIST /time_control/ auxhist2_interval_y
NAMELIST /time_control/ auxhist2_interval_d
NAMELIST /time_control/ auxhist2_interval_h
NAMELIST /time_control/ auxhist2_interval_m
NAMELIST /time_control/ auxhist2_interval_s
NAMELIST /time_control/ auxhist2_interval
NAMELIST /time_control/ auxhist2_begin_y
NAMELIST /time_control/ auxhist2_begin_d
NAMELIST /time_control/ auxhist2_begin_h
NAMELIST /time_control/ auxhist2_begin_m
NAMELIST /time_control/ auxhist2_begin_s
NAMELIST /time_control/ auxhist2_begin
NAMELIST /time_control/ auxhist2_end_y
NAMELIST /time_control/ auxhist2_end_d
NAMELIST /time_control/ auxhist2_end_h
NAMELIST /time_control/ auxhist2_end_m
NAMELIST /time_control/ auxhist2_end_s
NAMELIST /time_control/ auxhist2_end
NAMELIST /time_control/ io_form_auxhist2
NAMELIST /time_control/ frames_per_auxhist2
NAMELIST /time_control/ auxhist3_inname
NAMELIST /time_control/ auxhist3_outname
NAMELIST /time_control/ auxhist3_interval_y
NAMELIST /time_control/ auxhist3_interval_d
NAMELIST /time_control/ auxhist3_interval_h
NAMELIST /time_control/ auxhist3_interval_m
NAMELIST /time_control/ auxhist3_interval_s
NAMELIST /time_control/ auxhist3_interval
NAMELIST /time_control/ auxhist3_begin_y
NAMELIST /time_control/ auxhist3_begin_d
NAMELIST /time_control/ auxhist3_begin_h
NAMELIST /time_control/ auxhist3_begin_m
NAMELIST /time_control/ auxhist3_begin_s
NAMELIST /time_control/ auxhist3_begin
NAMELIST /time_control/ auxhist3_end_y
NAMELIST /time_control/ auxhist3_end_d
NAMELIST /time_control/ auxhist3_end_h
NAMELIST /time_control/ auxhist3_end_m
NAMELIST /time_control/ auxhist3_end_s
NAMELIST /time_control/ auxhist3_end
NAMELIST /time_control/ io_form_auxhist3
NAMELIST /time_control/ frames_per_auxhist3
NAMELIST /time_control/ auxhist4_inname
NAMELIST /time_control/ auxhist4_outname
NAMELIST /time_control/ auxhist4_interval_y
NAMELIST /time_control/ auxhist4_interval_d
NAMELIST /time_control/ auxhist4_interval_h
NAMELIST /time_control/ auxhist4_interval_m
NAMELIST /time_control/ auxhist4_interval_s
NAMELIST /time_control/ auxhist4_interval
NAMELIST /time_control/ auxhist4_begin_y
NAMELIST /time_control/ auxhist4_begin_d
NAMELIST /time_control/ auxhist4_begin_h
NAMELIST /time_control/ auxhist4_begin_m
NAMELIST /time_control/ auxhist4_begin_s
NAMELIST /time_control/ auxhist4_begin
NAMELIST /time_control/ auxhist4_end_y
NAMELIST /time_control/ auxhist4_end_d
NAMELIST /time_control/ auxhist4_end_h
NAMELIST /time_control/ auxhist4_end_m
NAMELIST /time_control/ auxhist4_end_s
NAMELIST /time_control/ auxhist4_end
NAMELIST /time_control/ io_form_auxhist4
NAMELIST /time_control/ frames_per_auxhist4
NAMELIST /time_control/ auxhist5_inname
NAMELIST /time_control/ auxhist5_outname
NAMELIST /time_control/ auxhist5_interval_y
NAMELIST /time_control/ auxhist5_interval_d
NAMELIST /time_control/ auxhist5_interval_h
NAMELIST /time_control/ auxhist5_interval_m
NAMELIST /time_control/ auxhist5_interval_s
NAMELIST /time_control/ auxhist5_interval
NAMELIST /time_control/ auxhist5_begin_y
NAMELIST /time_control/ auxhist5_begin_d
NAMELIST /time_control/ auxhist5_begin_h
NAMELIST /time_control/ auxhist5_begin_m
NAMELIST /time_control/ auxhist5_begin_s
NAMELIST /time_control/ auxhist5_begin
NAMELIST /time_control/ auxhist5_end_y
NAMELIST /time_control/ auxhist5_end_d
NAMELIST /time_control/ auxhist5_end_h
NAMELIST /time_control/ auxhist5_end_m
NAMELIST /time_control/ auxhist5_end_s
NAMELIST /time_control/ auxhist5_end
NAMELIST /time_control/ io_form_auxhist5
NAMELIST /time_control/ frames_per_auxhist5
NAMELIST /time_control/ auxhist6_inname
NAMELIST /time_control/ auxhist6_outname
NAMELIST /time_control/ auxhist6_interval_y
NAMELIST /time_control/ auxhist6_interval_d
NAMELIST /time_control/ auxhist6_interval_h
NAMELIST /time_control/ auxhist6_interval_m
NAMELIST /time_control/ auxhist6_interval_s
NAMELIST /time_control/ auxhist6_interval
NAMELIST /time_control/ auxhist6_begin_y
NAMELIST /time_control/ auxhist6_begin_d
NAMELIST /time_control/ auxhist6_begin_h
NAMELIST /time_control/ auxhist6_begin_m
NAMELIST /time_control/ auxhist6_begin_s
NAMELIST /time_control/ auxhist6_begin
NAMELIST /time_control/ auxhist6_end_y
NAMELIST /time_control/ auxhist6_end_d
NAMELIST /time_control/ auxhist6_end_h
NAMELIST /time_control/ auxhist6_end_m
NAMELIST /time_control/ auxhist6_end_s
NAMELIST /time_control/ auxhist6_end
NAMELIST /time_control/ io_form_auxhist6
NAMELIST /time_control/ frames_per_auxhist6
NAMELIST /time_control/ auxhist7_inname
NAMELIST /time_control/ auxhist7_outname
NAMELIST /time_control/ auxhist7_interval_y
NAMELIST /time_control/ auxhist7_interval_d
NAMELIST /time_control/ auxhist7_interval_h
NAMELIST /time_control/ auxhist7_interval_m
NAMELIST /time_control/ auxhist7_interval_s
NAMELIST /time_control/ auxhist7_interval
NAMELIST /time_control/ auxhist7_begin_y
NAMELIST /time_control/ auxhist7_begin_d
NAMELIST /time_control/ auxhist7_begin_h
NAMELIST /time_control/ auxhist7_begin_m
NAMELIST /time_control/ auxhist7_begin_s
NAMELIST /time_control/ auxhist7_begin
NAMELIST /time_control/ auxhist7_end_y
NAMELIST /time_control/ auxhist7_end_d
NAMELIST /time_control/ auxhist7_end_h
NAMELIST /time_control/ auxhist7_end_m
NAMELIST /time_control/ auxhist7_end_s
NAMELIST /time_control/ auxhist7_end
NAMELIST /time_control/ io_form_auxhist7
NAMELIST /time_control/ frames_per_auxhist7
NAMELIST /time_control/ auxhist8_inname
NAMELIST /time_control/ auxhist8_outname
NAMELIST /time_control/ auxhist8_interval_y
NAMELIST /time_control/ auxhist8_interval_d
NAMELIST /time_control/ auxhist8_interval_h
NAMELIST /time_control/ auxhist8_interval_m
NAMELIST /time_control/ auxhist8_interval_s
NAMELIST /time_control/ auxhist8_interval
NAMELIST /time_control/ auxhist8_begin_y
NAMELIST /time_control/ auxhist8_begin_d
NAMELIST /time_control/ auxhist8_begin_h
NAMELIST /time_control/ auxhist8_begin_m
NAMELIST /time_control/ auxhist8_begin_s
NAMELIST /time_control/ auxhist8_begin
NAMELIST /time_control/ auxhist8_end_y
NAMELIST /time_control/ auxhist8_end_d
NAMELIST /time_control/ auxhist8_end_h
NAMELIST /time_control/ auxhist8_end_m
NAMELIST /time_control/ auxhist8_end_s
NAMELIST /time_control/ auxhist8_end
NAMELIST /time_control/ io_form_auxhist8
NAMELIST /time_control/ frames_per_auxhist8
NAMELIST /time_control/ auxhist9_inname
NAMELIST /time_control/ auxhist9_outname
NAMELIST /time_control/ auxhist9_interval_y
NAMELIST /time_control/ auxhist9_interval_d
NAMELIST /time_control/ auxhist9_interval_h
NAMELIST /time_control/ auxhist9_interval_m
NAMELIST /time_control/ auxhist9_interval_s
NAMELIST /time_control/ auxhist9_interval
NAMELIST /time_control/ auxhist9_begin_y
NAMELIST /time_control/ auxhist9_begin_d
NAMELIST /time_control/ auxhist9_begin_h
NAMELIST /time_control/ auxhist9_begin_m
NAMELIST /time_control/ auxhist9_begin_s
NAMELIST /time_control/ auxhist9_begin
NAMELIST /time_control/ auxhist9_end_y
NAMELIST /time_control/ auxhist9_end_d
NAMELIST /time_control/ auxhist9_end_h
NAMELIST /time_control/ auxhist9_end_m
NAMELIST /time_control/ auxhist9_end_s
NAMELIST /time_control/ auxhist9_end
NAMELIST /time_control/ io_form_auxhist9
NAMELIST /time_control/ frames_per_auxhist9
NAMELIST /time_control/ auxhist10_inname
NAMELIST /time_control/ auxhist10_outname
NAMELIST /time_control/ auxhist10_interval_y
NAMELIST /time_control/ auxhist10_interval_d
NAMELIST /time_control/ auxhist10_interval_h
NAMELIST /time_control/ auxhist10_interval_m
NAMELIST /time_control/ auxhist10_interval_s
NAMELIST /time_control/ auxhist10_interval
NAMELIST /time_control/ auxhist10_begin_y
NAMELIST /time_control/ auxhist10_begin_d
NAMELIST /time_control/ auxhist10_begin_h
NAMELIST /time_control/ auxhist10_begin_m
NAMELIST /time_control/ auxhist10_begin_s
NAMELIST /time_control/ auxhist10_begin
NAMELIST /time_control/ auxhist10_end_y
NAMELIST /time_control/ auxhist10_end_d
NAMELIST /time_control/ auxhist10_end_h
NAMELIST /time_control/ auxhist10_end_m
NAMELIST /time_control/ auxhist10_end_s
NAMELIST /time_control/ auxhist10_end
NAMELIST /time_control/ io_form_auxhist10
NAMELIST /time_control/ frames_per_auxhist10
NAMELIST /time_control/ auxhist11_inname
NAMELIST /time_control/ auxhist11_outname
NAMELIST /time_control/ auxhist11_interval_y
NAMELIST /time_control/ auxhist11_interval_d
NAMELIST /time_control/ auxhist11_interval_h
NAMELIST /time_control/ auxhist11_interval_m
NAMELIST /time_control/ auxhist11_interval_s
NAMELIST /time_control/ auxhist11_interval
NAMELIST /time_control/ auxhist11_begin_y
NAMELIST /time_control/ auxhist11_begin_d
NAMELIST /time_control/ auxhist11_begin_h
NAMELIST /time_control/ auxhist11_begin_m
NAMELIST /time_control/ auxhist11_begin_s
NAMELIST /time_control/ auxhist11_begin
NAMELIST /time_control/ auxhist11_end_y
NAMELIST /time_control/ auxhist11_end_d
NAMELIST /time_control/ auxhist11_end_h
NAMELIST /time_control/ auxhist11_end_m
NAMELIST /time_control/ auxhist11_end_s
NAMELIST /time_control/ auxhist11_end
NAMELIST /time_control/ io_form_auxhist11
NAMELIST /time_control/ frames_per_auxhist11
NAMELIST /time_control/ auxhist12_inname
NAMELIST /time_control/ auxhist12_outname
NAMELIST /time_control/ auxhist12_interval_y
NAMELIST /time_control/ auxhist12_interval_d
NAMELIST /time_control/ auxhist12_interval_h
NAMELIST /time_control/ auxhist12_interval_m
NAMELIST /time_control/ auxhist12_interval_s
NAMELIST /time_control/ auxhist12_interval
NAMELIST /time_control/ auxhist12_begin_y
NAMELIST /time_control/ auxhist12_begin_d
NAMELIST /time_control/ auxhist12_begin_h
NAMELIST /time_control/ auxhist12_begin_m
NAMELIST /time_control/ auxhist12_begin_s
NAMELIST /time_control/ auxhist12_begin
NAMELIST /time_control/ auxhist12_end_y
NAMELIST /time_control/ auxhist12_end_d
NAMELIST /time_control/ auxhist12_end_h
NAMELIST /time_control/ auxhist12_end_m
NAMELIST /time_control/ auxhist12_end_s
NAMELIST /time_control/ auxhist12_end
NAMELIST /time_control/ io_form_auxhist12
NAMELIST /time_control/ frames_per_auxhist12
NAMELIST /time_control/ auxhist13_inname
NAMELIST /time_control/ auxhist13_outname
NAMELIST /time_control/ auxhist13_interval_y
NAMELIST /time_control/ auxhist13_interval_d
NAMELIST /time_control/ auxhist13_interval_h
NAMELIST /time_control/ auxhist13_interval_m
NAMELIST /time_control/ auxhist13_interval_s
NAMELIST /time_control/ auxhist13_interval
NAMELIST /time_control/ auxhist13_begin_y
NAMELIST /time_control/ auxhist13_begin_d
NAMELIST /time_control/ auxhist13_begin_h
NAMELIST /time_control/ auxhist13_begin_m
NAMELIST /time_control/ auxhist13_begin_s
NAMELIST /time_control/ auxhist13_begin
NAMELIST /time_control/ auxhist13_end_y
NAMELIST /time_control/ auxhist13_end_d
NAMELIST /time_control/ auxhist13_end_h
NAMELIST /time_control/ auxhist13_end_m
NAMELIST /time_control/ auxhist13_end_s
NAMELIST /time_control/ auxhist13_end
NAMELIST /time_control/ io_form_auxhist13
NAMELIST /time_control/ frames_per_auxhist13
NAMELIST /time_control/ auxhist14_inname
NAMELIST /time_control/ auxhist14_outname
NAMELIST /time_control/ auxhist14_interval_y
NAMELIST /time_control/ auxhist14_interval_d
NAMELIST /time_control/ auxhist14_interval_h
NAMELIST /time_control/ auxhist14_interval_m
NAMELIST /time_control/ auxhist14_interval_s
NAMELIST /time_control/ auxhist14_interval
NAMELIST /time_control/ auxhist14_begin_y
NAMELIST /time_control/ auxhist14_begin_d
NAMELIST /time_control/ auxhist14_begin_h
NAMELIST /time_control/ auxhist14_begin_m
NAMELIST /time_control/ auxhist14_begin_s
NAMELIST /time_control/ auxhist14_begin
NAMELIST /time_control/ auxhist14_end_y
NAMELIST /time_control/ auxhist14_end_d
NAMELIST /time_control/ auxhist14_end_h
NAMELIST /time_control/ auxhist14_end_m
NAMELIST /time_control/ auxhist14_end_s
NAMELIST /time_control/ auxhist14_end
NAMELIST /time_control/ io_form_auxhist14
NAMELIST /time_control/ frames_per_auxhist14
NAMELIST /time_control/ auxhist15_inname
NAMELIST /time_control/ auxhist15_outname
NAMELIST /time_control/ auxhist15_interval_y
NAMELIST /time_control/ auxhist15_interval_d
NAMELIST /time_control/ auxhist15_interval_h
NAMELIST /time_control/ auxhist15_interval_m
NAMELIST /time_control/ auxhist15_interval_s
NAMELIST /time_control/ auxhist15_interval
NAMELIST /time_control/ auxhist15_begin_y
NAMELIST /time_control/ auxhist15_begin_d
NAMELIST /time_control/ auxhist15_begin_h
NAMELIST /time_control/ auxhist15_begin_m
NAMELIST /time_control/ auxhist15_begin_s
NAMELIST /time_control/ auxhist15_begin
NAMELIST /time_control/ auxhist15_end_y
NAMELIST /time_control/ auxhist15_end_d
NAMELIST /time_control/ auxhist15_end_h
NAMELIST /time_control/ auxhist15_end_m
NAMELIST /time_control/ auxhist15_end_s
NAMELIST /time_control/ auxhist15_end
NAMELIST /time_control/ io_form_auxhist15
NAMELIST /time_control/ frames_per_auxhist15
NAMELIST /time_control/ auxhist16_inname
NAMELIST /time_control/ auxhist16_outname
NAMELIST /time_control/ auxhist16_interval_y
NAMELIST /time_control/ auxhist16_interval_d
NAMELIST /time_control/ auxhist16_interval_h
NAMELIST /time_control/ auxhist16_interval_m
NAMELIST /time_control/ auxhist16_interval_s
NAMELIST /time_control/ auxhist16_interval
NAMELIST /time_control/ auxhist16_begin_y
NAMELIST /time_control/ auxhist16_begin_d
NAMELIST /time_control/ auxhist16_begin_h
NAMELIST /time_control/ auxhist16_begin_m
NAMELIST /time_control/ auxhist16_begin_s
NAMELIST /time_control/ auxhist16_begin
NAMELIST /time_control/ auxhist16_end_y
NAMELIST /time_control/ auxhist16_end_d
NAMELIST /time_control/ auxhist16_end_h
NAMELIST /time_control/ auxhist16_end_m
NAMELIST /time_control/ auxhist16_end_s
NAMELIST /time_control/ auxhist16_end
NAMELIST /time_control/ io_form_auxhist16
NAMELIST /time_control/ frames_per_auxhist16
NAMELIST /time_control/ auxhist17_inname
NAMELIST /time_control/ auxhist17_outname
NAMELIST /time_control/ auxhist17_interval_y
NAMELIST /time_control/ auxhist17_interval_d
NAMELIST /time_control/ auxhist17_interval_h
NAMELIST /time_control/ auxhist17_interval_m
NAMELIST /time_control/ auxhist17_interval_s
NAMELIST /time_control/ auxhist17_interval
NAMELIST /time_control/ auxhist17_begin_y
NAMELIST /time_control/ auxhist17_begin_d
NAMELIST /time_control/ auxhist17_begin_h
NAMELIST /time_control/ auxhist17_begin_m
NAMELIST /time_control/ auxhist17_begin_s
NAMELIST /time_control/ auxhist17_begin
NAMELIST /time_control/ auxhist17_end_y
NAMELIST /time_control/ auxhist17_end_d
NAMELIST /time_control/ auxhist17_end_h
NAMELIST /time_control/ auxhist17_end_m
NAMELIST /time_control/ auxhist17_end_s
NAMELIST /time_control/ auxhist17_end
NAMELIST /time_control/ io_form_auxhist17
NAMELIST /time_control/ frames_per_auxhist17
NAMELIST /time_control/ auxhist18_inname
NAMELIST /time_control/ auxhist18_outname
NAMELIST /time_control/ auxhist18_interval_y
NAMELIST /time_control/ auxhist18_interval_d
NAMELIST /time_control/ auxhist18_interval_h
NAMELIST /time_control/ auxhist18_interval_m
NAMELIST /time_control/ auxhist18_interval_s
NAMELIST /time_control/ auxhist18_interval
NAMELIST /time_control/ auxhist18_begin_y
NAMELIST /time_control/ auxhist18_begin_d
NAMELIST /time_control/ auxhist18_begin_h
NAMELIST /time_control/ auxhist18_begin_m
NAMELIST /time_control/ auxhist18_begin_s
NAMELIST /time_control/ auxhist18_begin
NAMELIST /time_control/ auxhist18_end_y
NAMELIST /time_control/ auxhist18_end_d
NAMELIST /time_control/ auxhist18_end_h
NAMELIST /time_control/ auxhist18_end_m
NAMELIST /time_control/ auxhist18_end_s
NAMELIST /time_control/ auxhist18_end
NAMELIST /time_control/ io_form_auxhist18
NAMELIST /time_control/ frames_per_auxhist18
NAMELIST /time_control/ auxhist19_inname
NAMELIST /time_control/ auxhist19_outname
NAMELIST /time_control/ auxhist19_interval_y
NAMELIST /time_control/ auxhist19_interval_d
NAMELIST /time_control/ auxhist19_interval_h
NAMELIST /time_control/ auxhist19_interval_m
NAMELIST /time_control/ auxhist19_interval_s
NAMELIST /time_control/ auxhist19_interval
NAMELIST /time_control/ auxhist19_begin_y
NAMELIST /time_control/ auxhist19_begin_d
NAMELIST /time_control/ auxhist19_begin_h
NAMELIST /time_control/ auxhist19_begin_m
NAMELIST /time_control/ auxhist19_begin_s
NAMELIST /time_control/ auxhist19_begin
NAMELIST /time_control/ auxhist19_end_y
NAMELIST /time_control/ auxhist19_end_d
NAMELIST /time_control/ auxhist19_end_h
NAMELIST /time_control/ auxhist19_end_m
NAMELIST /time_control/ auxhist19_end_s
NAMELIST /time_control/ auxhist19_end
NAMELIST /time_control/ io_form_auxhist19
NAMELIST /time_control/ frames_per_auxhist19
NAMELIST /time_control/ auxhist20_inname
NAMELIST /time_control/ auxhist20_outname
NAMELIST /time_control/ auxhist20_interval_y
NAMELIST /time_control/ auxhist20_interval_d
NAMELIST /time_control/ auxhist20_interval_h
NAMELIST /time_control/ auxhist20_interval_m
NAMELIST /time_control/ auxhist20_interval_s
NAMELIST /time_control/ auxhist20_interval
NAMELIST /time_control/ auxhist20_begin_y
NAMELIST /time_control/ auxhist20_begin_d
NAMELIST /time_control/ auxhist20_begin_h
NAMELIST /time_control/ auxhist20_begin_m
NAMELIST /time_control/ auxhist20_begin_s
NAMELIST /time_control/ auxhist20_begin
NAMELIST /time_control/ auxhist20_end_y
NAMELIST /time_control/ auxhist20_end_d
NAMELIST /time_control/ auxhist20_end_h
NAMELIST /time_control/ auxhist20_end_m
NAMELIST /time_control/ auxhist20_end_s
NAMELIST /time_control/ auxhist20_end
NAMELIST /time_control/ io_form_auxhist20
NAMELIST /time_control/ frames_per_auxhist20
NAMELIST /time_control/ auxhist21_inname
NAMELIST /time_control/ auxhist21_outname
NAMELIST /time_control/ auxhist21_interval_y
NAMELIST /time_control/ auxhist21_interval_d
NAMELIST /time_control/ auxhist21_interval_h
NAMELIST /time_control/ auxhist21_interval_m
NAMELIST /time_control/ auxhist21_interval_s
NAMELIST /time_control/ auxhist21_interval
NAMELIST /time_control/ auxhist21_begin_y
NAMELIST /time_control/ auxhist21_begin_d
NAMELIST /time_control/ auxhist21_begin_h
NAMELIST /time_control/ auxhist21_begin_m
NAMELIST /time_control/ auxhist21_begin_s
NAMELIST /time_control/ auxhist21_begin
NAMELIST /time_control/ auxhist21_end_y
NAMELIST /time_control/ auxhist21_end_d
NAMELIST /time_control/ auxhist21_end_h
NAMELIST /time_control/ auxhist21_end_m
NAMELIST /time_control/ auxhist21_end_s
NAMELIST /time_control/ auxhist21_end
NAMELIST /time_control/ io_form_auxhist21
NAMELIST /time_control/ frames_per_auxhist21
NAMELIST /time_control/ auxhist22_inname
NAMELIST /time_control/ auxhist22_outname
NAMELIST /time_control/ auxhist22_interval_y
NAMELIST /time_control/ auxhist22_interval_d
NAMELIST /time_control/ auxhist22_interval_h
NAMELIST /time_control/ auxhist22_interval_m
NAMELIST /time_control/ auxhist22_interval_s
NAMELIST /time_control/ auxhist22_interval
NAMELIST /time_control/ auxhist22_begin_y
NAMELIST /time_control/ auxhist22_begin_d
NAMELIST /time_control/ auxhist22_begin_h
NAMELIST /time_control/ auxhist22_begin_m
NAMELIST /time_control/ auxhist22_begin_s
NAMELIST /time_control/ auxhist22_begin
NAMELIST /time_control/ auxhist22_end_y
NAMELIST /time_control/ auxhist22_end_d
NAMELIST /time_control/ auxhist22_end_h
NAMELIST /time_control/ auxhist22_end_m
NAMELIST /time_control/ auxhist22_end_s
NAMELIST /time_control/ auxhist22_end
NAMELIST /time_control/ io_form_auxhist22
NAMELIST /time_control/ frames_per_auxhist22
NAMELIST /time_control/ auxhist23_inname
NAMELIST /time_control/ auxhist23_outname
NAMELIST /time_control/ auxhist23_interval_y
NAMELIST /time_control/ auxhist23_interval_d
NAMELIST /time_control/ auxhist23_interval_h
NAMELIST /time_control/ auxhist23_interval_m
NAMELIST /time_control/ auxhist23_interval_s
NAMELIST /time_control/ auxhist23_interval
NAMELIST /time_control/ auxhist23_begin_y
NAMELIST /time_control/ auxhist23_begin_d
NAMELIST /time_control/ auxhist23_begin_h
NAMELIST /time_control/ auxhist23_begin_m
NAMELIST /time_control/ auxhist23_begin_s
NAMELIST /time_control/ auxhist23_begin
NAMELIST /time_control/ auxhist23_end_y
NAMELIST /time_control/ auxhist23_end_d
NAMELIST /time_control/ auxhist23_end_h
NAMELIST /time_control/ auxhist23_end_m
NAMELIST /time_control/ auxhist23_end_s
NAMELIST /time_control/ auxhist23_end
NAMELIST /time_control/ io_form_auxhist23
NAMELIST /time_control/ frames_per_auxhist23
NAMELIST /time_control/ auxhist24_inname
NAMELIST /time_control/ auxhist24_outname
NAMELIST /time_control/ auxhist24_interval_y
NAMELIST /time_control/ auxhist24_interval_d
NAMELIST /time_control/ auxhist24_interval_h
NAMELIST /time_control/ auxhist24_interval_m
NAMELIST /time_control/ auxhist24_interval_s
NAMELIST /time_control/ auxhist24_interval
NAMELIST /time_control/ auxhist24_begin_y
NAMELIST /time_control/ auxhist24_begin_d
NAMELIST /time_control/ auxhist24_begin_h
NAMELIST /time_control/ auxhist24_begin_m
NAMELIST /time_control/ auxhist24_begin_s
NAMELIST /time_control/ auxhist24_begin
NAMELIST /time_control/ auxhist24_end_y
NAMELIST /time_control/ auxhist24_end_d
NAMELIST /time_control/ auxhist24_end_h
NAMELIST /time_control/ auxhist24_end_m
NAMELIST /time_control/ auxhist24_end_s
NAMELIST /time_control/ auxhist24_end
NAMELIST /time_control/ io_form_auxhist24
NAMELIST /time_control/ frames_per_auxhist24
NAMELIST /time_control/ auxinput1_outname
NAMELIST /time_control/ auxinput1_interval_y
NAMELIST /time_control/ auxinput1_interval_d
NAMELIST /time_control/ auxinput1_interval_h
NAMELIST /time_control/ auxinput1_interval_m
NAMELIST /time_control/ auxinput1_interval_s
NAMELIST /time_control/ auxinput1_interval
NAMELIST /time_control/ auxinput1_begin_y
NAMELIST /time_control/ auxinput1_begin_d
NAMELIST /time_control/ auxinput1_begin_h
NAMELIST /time_control/ auxinput1_begin_m
NAMELIST /time_control/ auxinput1_begin_s
NAMELIST /time_control/ auxinput1_begin
NAMELIST /time_control/ auxinput1_end_y
NAMELIST /time_control/ auxinput1_end_d
NAMELIST /time_control/ auxinput1_end_h
NAMELIST /time_control/ auxinput1_end_m
NAMELIST /time_control/ auxinput1_end_s
NAMELIST /time_control/ auxinput1_end
NAMELIST /time_control/ frames_per_auxinput1
NAMELIST /time_control/ auxinput2_inname
NAMELIST /time_control/ auxinput2_outname
NAMELIST /time_control/ auxinput2_interval_y
NAMELIST /time_control/ auxinput2_interval_d
NAMELIST /time_control/ auxinput2_interval_h
NAMELIST /time_control/ auxinput2_interval_m
NAMELIST /time_control/ auxinput2_interval_s
NAMELIST /time_control/ auxinput2_interval
NAMELIST /time_control/ auxinput2_begin_y
NAMELIST /time_control/ auxinput2_begin_d
NAMELIST /time_control/ auxinput2_begin_h
NAMELIST /time_control/ auxinput2_begin_m
NAMELIST /time_control/ auxinput2_begin_s
NAMELIST /time_control/ auxinput2_begin
NAMELIST /time_control/ auxinput2_end_y
NAMELIST /time_control/ auxinput2_end_d
NAMELIST /time_control/ auxinput2_end_h
NAMELIST /time_control/ auxinput2_end_m
NAMELIST /time_control/ auxinput2_end_s
NAMELIST /time_control/ auxinput2_end
NAMELIST /time_control/ io_form_auxinput2
NAMELIST /time_control/ frames_per_auxinput2
NAMELIST /time_control/ auxinput3_inname
NAMELIST /time_control/ auxinput3_outname
NAMELIST /time_control/ auxinput3_interval_y
NAMELIST /time_control/ auxinput3_interval_d
NAMELIST /time_control/ auxinput3_interval_h
NAMELIST /time_control/ auxinput3_interval_m
NAMELIST /time_control/ auxinput3_interval_s
NAMELIST /time_control/ auxinput3_interval
NAMELIST /time_control/ auxinput3_begin_y
NAMELIST /time_control/ auxinput3_begin_d
NAMELIST /time_control/ auxinput3_begin_h
NAMELIST /time_control/ auxinput3_begin_m
NAMELIST /time_control/ auxinput3_begin_s
NAMELIST /time_control/ auxinput3_begin
NAMELIST /time_control/ auxinput3_end_y
NAMELIST /time_control/ auxinput3_end_d
NAMELIST /time_control/ auxinput3_end_h
NAMELIST /time_control/ auxinput3_end_m
NAMELIST /time_control/ auxinput3_end_s
NAMELIST /time_control/ auxinput3_end
NAMELIST /time_control/ io_form_auxinput3
NAMELIST /time_control/ frames_per_auxinput3
NAMELIST /time_control/ auxinput4_inname
NAMELIST /time_control/ auxinput4_outname
NAMELIST /time_control/ auxinput4_interval_y
NAMELIST /time_control/ auxinput4_interval_d
NAMELIST /time_control/ auxinput4_interval_h
NAMELIST /time_control/ auxinput4_interval_m
NAMELIST /time_control/ auxinput4_interval_s
NAMELIST /time_control/ auxinput4_interval
NAMELIST /time_control/ auxinput4_begin_y
NAMELIST /time_control/ auxinput4_begin_d
NAMELIST /time_control/ auxinput4_begin_h
NAMELIST /time_control/ auxinput4_begin_m
NAMELIST /time_control/ auxinput4_begin_s
NAMELIST /time_control/ auxinput4_begin
NAMELIST /time_control/ auxinput4_end_y
NAMELIST /time_control/ auxinput4_end_d
NAMELIST /time_control/ auxinput4_end_h
NAMELIST /time_control/ auxinput4_end_m
NAMELIST /time_control/ auxinput4_end_s
NAMELIST /time_control/ auxinput4_end
NAMELIST /time_control/ io_form_auxinput4
NAMELIST /time_control/ frames_per_auxinput4
NAMELIST /time_control/ auxinput5_inname
NAMELIST /time_control/ auxinput5_outname
NAMELIST /time_control/ auxinput5_interval_y
NAMELIST /time_control/ auxinput5_interval_d
NAMELIST /time_control/ auxinput5_interval_h
NAMELIST /time_control/ auxinput5_interval_m
NAMELIST /time_control/ auxinput5_interval_s
NAMELIST /time_control/ auxinput5_interval
NAMELIST /time_control/ auxinput5_begin_y
NAMELIST /time_control/ auxinput5_begin_d
NAMELIST /time_control/ auxinput5_begin_h
NAMELIST /time_control/ auxinput5_begin_m
NAMELIST /time_control/ auxinput5_begin_s
NAMELIST /time_control/ auxinput5_begin
NAMELIST /time_control/ auxinput5_end_y
NAMELIST /time_control/ auxinput5_end_d
NAMELIST /time_control/ auxinput5_end_h
NAMELIST /time_control/ auxinput5_end_m
NAMELIST /time_control/ auxinput5_end_s
NAMELIST /time_control/ auxinput5_end
NAMELIST /time_control/ io_form_auxinput5
NAMELIST /time_control/ frames_per_auxinput5
NAMELIST /time_control/ auxinput6_inname
NAMELIST /time_control/ auxinput6_outname
NAMELIST /time_control/ auxinput6_interval_y
NAMELIST /time_control/ auxinput6_interval_d
NAMELIST /time_control/ auxinput6_interval_h
NAMELIST /time_control/ auxinput6_interval_m
NAMELIST /time_control/ auxinput6_interval_s
NAMELIST /time_control/ auxinput6_interval
NAMELIST /time_control/ auxinput6_begin_y
NAMELIST /time_control/ auxinput6_begin_d
NAMELIST /time_control/ auxinput6_begin_h
NAMELIST /time_control/ auxinput6_begin_m
NAMELIST /time_control/ auxinput6_begin_s
NAMELIST /time_control/ auxinput6_begin
NAMELIST /time_control/ auxinput6_end_y
NAMELIST /time_control/ auxinput6_end_d
NAMELIST /time_control/ auxinput6_end_h
NAMELIST /time_control/ auxinput6_end_m
NAMELIST /time_control/ auxinput6_end_s
NAMELIST /time_control/ auxinput6_end
NAMELIST /time_control/ io_form_auxinput6
NAMELIST /time_control/ frames_per_auxinput6
NAMELIST /time_control/ auxinput7_inname
NAMELIST /time_control/ auxinput7_outname
NAMELIST /time_control/ auxinput7_interval_y
NAMELIST /time_control/ auxinput7_interval_d
NAMELIST /time_control/ auxinput7_interval_h
NAMELIST /time_control/ auxinput7_interval_m
NAMELIST /time_control/ auxinput7_interval_s
NAMELIST /time_control/ auxinput7_interval
NAMELIST /time_control/ auxinput7_begin_y
NAMELIST /time_control/ auxinput7_begin_d
NAMELIST /time_control/ auxinput7_begin_h
NAMELIST /time_control/ auxinput7_begin_m
NAMELIST /time_control/ auxinput7_begin_s
NAMELIST /time_control/ auxinput7_begin
NAMELIST /time_control/ auxinput7_end_y
NAMELIST /time_control/ auxinput7_end_d
NAMELIST /time_control/ auxinput7_end_h
NAMELIST /time_control/ auxinput7_end_m
NAMELIST /time_control/ auxinput7_end_s
NAMELIST /time_control/ auxinput7_end
NAMELIST /time_control/ io_form_auxinput7
NAMELIST /time_control/ frames_per_auxinput7
NAMELIST /time_control/ auxinput8_inname
NAMELIST /time_control/ auxinput8_outname
NAMELIST /time_control/ auxinput8_interval_y
NAMELIST /time_control/ auxinput8_interval_d
NAMELIST /time_control/ auxinput8_interval_h
NAMELIST /time_control/ auxinput8_interval_m
NAMELIST /time_control/ auxinput8_interval_s
NAMELIST /time_control/ auxinput8_interval
NAMELIST /time_control/ auxinput8_begin_y
NAMELIST /time_control/ auxinput8_begin_d
NAMELIST /time_control/ auxinput8_begin_h
NAMELIST /time_control/ auxinput8_begin_m
NAMELIST /time_control/ auxinput8_begin_s
NAMELIST /time_control/ auxinput8_begin
NAMELIST /time_control/ auxinput8_end_y
NAMELIST /time_control/ auxinput8_end_d
NAMELIST /time_control/ auxinput8_end_h
NAMELIST /time_control/ auxinput8_end_m
NAMELIST /time_control/ auxinput8_end_s
NAMELIST /time_control/ auxinput8_end
NAMELIST /time_control/ io_form_auxinput8
NAMELIST /time_control/ frames_per_auxinput8
NAMELIST /time_control/ auxinput9_inname
NAMELIST /time_control/ auxinput9_outname
NAMELIST /time_control/ auxinput9_interval_y
NAMELIST /time_control/ auxinput9_interval_d
NAMELIST /time_control/ auxinput9_interval_h
NAMELIST /time_control/ auxinput9_interval_m
NAMELIST /time_control/ auxinput9_interval_s
NAMELIST /time_control/ auxinput9_interval
NAMELIST /time_control/ auxinput9_begin_y
NAMELIST /time_control/ auxinput9_begin_d
NAMELIST /time_control/ auxinput9_begin_h
NAMELIST /time_control/ auxinput9_begin_m
NAMELIST /time_control/ auxinput9_begin_s
NAMELIST /time_control/ auxinput9_begin
NAMELIST /time_control/ auxinput9_end_y
NAMELIST /time_control/ auxinput9_end_d
NAMELIST /time_control/ auxinput9_end_h
NAMELIST /time_control/ auxinput9_end_m
NAMELIST /time_control/ auxinput9_end_s
NAMELIST /time_control/ auxinput9_end
NAMELIST /time_control/ io_form_auxinput9
NAMELIST /time_control/ frames_per_auxinput9
NAMELIST /time_control/ auxinput10_inname
NAMELIST /time_control/ auxinput10_outname
NAMELIST /time_control/ auxinput10_interval_y
NAMELIST /time_control/ auxinput10_interval_d
NAMELIST /time_control/ auxinput10_interval_h
NAMELIST /time_control/ auxinput10_interval_m
NAMELIST /time_control/ auxinput10_interval_s
NAMELIST /time_control/ auxinput10_interval
NAMELIST /time_control/ auxinput10_begin_y
NAMELIST /time_control/ auxinput10_begin_d
NAMELIST /time_control/ auxinput10_begin_h
NAMELIST /time_control/ auxinput10_begin_m
NAMELIST /time_control/ auxinput10_begin_s
NAMELIST /time_control/ auxinput10_begin
NAMELIST /time_control/ auxinput10_end_y
NAMELIST /time_control/ auxinput10_end_d
NAMELIST /time_control/ auxinput10_end_h
NAMELIST /time_control/ auxinput10_end_m
NAMELIST /time_control/ auxinput10_end_s
NAMELIST /time_control/ auxinput10_end
NAMELIST /time_control/ io_form_auxinput10
NAMELIST /time_control/ frames_per_auxinput10
NAMELIST /time_control/ auxinput11_inname
NAMELIST /time_control/ auxinput11_outname
NAMELIST /time_control/ auxinput11_interval_y
NAMELIST /time_control/ auxinput11_interval_d
NAMELIST /time_control/ auxinput11_interval_h
NAMELIST /time_control/ auxinput11_interval_m
NAMELIST /time_control/ auxinput11_interval_s
NAMELIST /time_control/ auxinput11_interval
NAMELIST /time_control/ auxinput11_begin_y
NAMELIST /time_control/ auxinput11_begin_d
NAMELIST /time_control/ auxinput11_begin_h
NAMELIST /time_control/ auxinput11_begin_m
NAMELIST /time_control/ auxinput11_begin_s
NAMELIST /time_control/ auxinput11_begin
NAMELIST /time_control/ auxinput11_end_y
NAMELIST /time_control/ auxinput11_end_d
NAMELIST /time_control/ auxinput11_end_h
NAMELIST /time_control/ auxinput11_end_m
NAMELIST /time_control/ auxinput11_end_s
NAMELIST /time_control/ auxinput11_end
NAMELIST /time_control/ io_form_auxinput11
NAMELIST /time_control/ frames_per_auxinput11
NAMELIST /time_control/ auxinput12_inname
NAMELIST /time_control/ auxinput12_outname
NAMELIST /time_control/ auxinput12_interval_y
NAMELIST /time_control/ auxinput12_interval_d
NAMELIST /time_control/ auxinput12_interval_h
NAMELIST /time_control/ auxinput12_interval_m
NAMELIST /time_control/ auxinput12_interval_s
NAMELIST /time_control/ auxinput12_interval
NAMELIST /time_control/ auxinput12_begin_y
NAMELIST /time_control/ auxinput12_begin_d
NAMELIST /time_control/ auxinput12_begin_h
NAMELIST /time_control/ auxinput12_begin_m
NAMELIST /time_control/ auxinput12_begin_s
NAMELIST /time_control/ auxinput12_begin
NAMELIST /time_control/ auxinput12_end_y
NAMELIST /time_control/ auxinput12_end_d
NAMELIST /time_control/ auxinput12_end_h
NAMELIST /time_control/ auxinput12_end_m
NAMELIST /time_control/ auxinput12_end_s
NAMELIST /time_control/ auxinput12_end
NAMELIST /time_control/ io_form_auxinput12
NAMELIST /time_control/ frames_per_auxinput12
NAMELIST /time_control/ auxinput13_inname
NAMELIST /time_control/ auxinput13_outname
NAMELIST /time_control/ auxinput13_interval_y
NAMELIST /time_control/ auxinput13_interval_d
NAMELIST /time_control/ auxinput13_interval_h
NAMELIST /time_control/ auxinput13_interval_m
NAMELIST /time_control/ auxinput13_interval_s
NAMELIST /time_control/ auxinput13_interval
NAMELIST /time_control/ auxinput13_begin_y
NAMELIST /time_control/ auxinput13_begin_d
NAMELIST /time_control/ auxinput13_begin_h
NAMELIST /time_control/ auxinput13_begin_m
NAMELIST /time_control/ auxinput13_begin_s
NAMELIST /time_control/ auxinput13_begin
NAMELIST /time_control/ auxinput13_end_y
NAMELIST /time_control/ auxinput13_end_d
NAMELIST /time_control/ auxinput13_end_h
NAMELIST /time_control/ auxinput13_end_m
NAMELIST /time_control/ auxinput13_end_s
NAMELIST /time_control/ auxinput13_end
NAMELIST /time_control/ io_form_auxinput13
NAMELIST /time_control/ frames_per_auxinput13
NAMELIST /time_control/ auxinput14_inname
NAMELIST /time_control/ auxinput14_outname
NAMELIST /time_control/ auxinput14_interval_y
NAMELIST /time_control/ auxinput14_interval_d
NAMELIST /time_control/ auxinput14_interval_h
NAMELIST /time_control/ auxinput14_interval_m
NAMELIST /time_control/ auxinput14_interval_s
NAMELIST /time_control/ auxinput14_interval
NAMELIST /time_control/ auxinput14_begin_y
NAMELIST /time_control/ auxinput14_begin_d
NAMELIST /time_control/ auxinput14_begin_h
NAMELIST /time_control/ auxinput14_begin_m
NAMELIST /time_control/ auxinput14_begin_s
NAMELIST /time_control/ auxinput14_begin
NAMELIST /time_control/ auxinput14_end_y
NAMELIST /time_control/ auxinput14_end_d
NAMELIST /time_control/ auxinput14_end_h
NAMELIST /time_control/ auxinput14_end_m
NAMELIST /time_control/ auxinput14_end_s
NAMELIST /time_control/ auxinput14_end
NAMELIST /time_control/ io_form_auxinput14
NAMELIST /time_control/ frames_per_auxinput14
NAMELIST /time_control/ auxinput15_inname
NAMELIST /time_control/ auxinput15_outname
NAMELIST /time_control/ auxinput15_interval_y
NAMELIST /time_control/ auxinput15_interval_d
NAMELIST /time_control/ auxinput15_interval_h
NAMELIST /time_control/ auxinput15_interval_m
NAMELIST /time_control/ auxinput15_interval_s
NAMELIST /time_control/ auxinput15_interval
NAMELIST /time_control/ auxinput15_begin_y
NAMELIST /time_control/ auxinput15_begin_d
NAMELIST /time_control/ auxinput15_begin_h
NAMELIST /time_control/ auxinput15_begin_m
NAMELIST /time_control/ auxinput15_begin_s
NAMELIST /time_control/ auxinput15_begin
NAMELIST /time_control/ auxinput15_end_y
NAMELIST /time_control/ auxinput15_end_d
NAMELIST /time_control/ auxinput15_end_h
NAMELIST /time_control/ auxinput15_end_m
NAMELIST /time_control/ auxinput15_end_s
NAMELIST /time_control/ auxinput15_end
NAMELIST /time_control/ io_form_auxinput15
NAMELIST /time_control/ frames_per_auxinput15
NAMELIST /time_control/ auxinput16_inname
NAMELIST /time_control/ auxinput16_outname
NAMELIST /time_control/ auxinput16_interval_y
NAMELIST /time_control/ auxinput16_interval_d
NAMELIST /time_control/ auxinput16_interval_h
NAMELIST /time_control/ auxinput16_interval_m
NAMELIST /time_control/ auxinput16_interval_s
NAMELIST /time_control/ auxinput16_interval
NAMELIST /time_control/ auxinput16_begin_y
NAMELIST /time_control/ auxinput16_begin_d
NAMELIST /time_control/ auxinput16_begin_h
NAMELIST /time_control/ auxinput16_begin_m
NAMELIST /time_control/ auxinput16_begin_s
NAMELIST /time_control/ auxinput16_begin
NAMELIST /time_control/ auxinput16_end_y
NAMELIST /time_control/ auxinput16_end_d
NAMELIST /time_control/ auxinput16_end_h
NAMELIST /time_control/ auxinput16_end_m
NAMELIST /time_control/ auxinput16_end_s
NAMELIST /time_control/ auxinput16_end
NAMELIST /time_control/ io_form_auxinput16
NAMELIST /time_control/ frames_per_auxinput16
NAMELIST /time_control/ auxinput17_inname
NAMELIST /time_control/ auxinput17_outname
NAMELIST /time_control/ auxinput17_interval_y
NAMELIST /time_control/ auxinput17_interval_d
NAMELIST /time_control/ auxinput17_interval_h
NAMELIST /time_control/ auxinput17_interval_m
NAMELIST /time_control/ auxinput17_interval_s
NAMELIST /time_control/ auxinput17_interval
NAMELIST /time_control/ auxinput17_begin_y
NAMELIST /time_control/ auxinput17_begin_d
NAMELIST /time_control/ auxinput17_begin_h
NAMELIST /time_control/ auxinput17_begin_m
NAMELIST /time_control/ auxinput17_begin_s
NAMELIST /time_control/ auxinput17_begin
NAMELIST /time_control/ auxinput17_end_y
NAMELIST /time_control/ auxinput17_end_d
NAMELIST /time_control/ auxinput17_end_h
NAMELIST /time_control/ auxinput17_end_m
NAMELIST /time_control/ auxinput17_end_s
NAMELIST /time_control/ auxinput17_end
NAMELIST /time_control/ io_form_auxinput17
NAMELIST /time_control/ frames_per_auxinput17
NAMELIST /time_control/ auxinput18_inname
NAMELIST /time_control/ auxinput18_outname
NAMELIST /time_control/ auxinput18_interval_y
NAMELIST /time_control/ auxinput18_interval_d
NAMELIST /time_control/ auxinput18_interval_h
NAMELIST /time_control/ auxinput18_interval_m
NAMELIST /time_control/ auxinput18_interval_s
NAMELIST /time_control/ auxinput18_interval
NAMELIST /time_control/ auxinput18_begin_y
NAMELIST /time_control/ auxinput18_begin_d
NAMELIST /time_control/ auxinput18_begin_h
NAMELIST /time_control/ auxinput18_begin_m
NAMELIST /time_control/ auxinput18_begin_s
NAMELIST /time_control/ auxinput18_begin
NAMELIST /time_control/ auxinput18_end_y
NAMELIST /time_control/ auxinput18_end_d
NAMELIST /time_control/ auxinput18_end_h
NAMELIST /time_control/ auxinput18_end_m
NAMELIST /time_control/ auxinput18_end_s
NAMELIST /time_control/ auxinput18_end
NAMELIST /time_control/ io_form_auxinput18
NAMELIST /time_control/ frames_per_auxinput18
NAMELIST /time_control/ auxinput19_inname
NAMELIST /time_control/ auxinput19_outname
NAMELIST /time_control/ auxinput19_interval_y
NAMELIST /time_control/ auxinput19_interval_d
NAMELIST /time_control/ auxinput19_interval_h
NAMELIST /time_control/ auxinput19_interval_m
NAMELIST /time_control/ auxinput19_interval_s
NAMELIST /time_control/ auxinput19_interval
NAMELIST /time_control/ auxinput19_begin_y
NAMELIST /time_control/ auxinput19_begin_d
NAMELIST /time_control/ auxinput19_begin_h
NAMELIST /time_control/ auxinput19_begin_m
NAMELIST /time_control/ auxinput19_begin_s
NAMELIST /time_control/ auxinput19_begin
NAMELIST /time_control/ auxinput19_end_y
NAMELIST /time_control/ auxinput19_end_d
NAMELIST /time_control/ auxinput19_end_h
NAMELIST /time_control/ auxinput19_end_m
NAMELIST /time_control/ auxinput19_end_s
NAMELIST /time_control/ auxinput19_end
NAMELIST /time_control/ io_form_auxinput19
NAMELIST /time_control/ frames_per_auxinput19
NAMELIST /time_control/ auxinput20_inname
NAMELIST /time_control/ auxinput20_outname
NAMELIST /time_control/ auxinput20_interval_y
NAMELIST /time_control/ auxinput20_interval_d
NAMELIST /time_control/ auxinput20_interval_h
NAMELIST /time_control/ auxinput20_interval_m
NAMELIST /time_control/ auxinput20_interval_s
NAMELIST /time_control/ auxinput20_interval
NAMELIST /time_control/ auxinput20_begin_y
NAMELIST /time_control/ auxinput20_begin_d
NAMELIST /time_control/ auxinput20_begin_h
NAMELIST /time_control/ auxinput20_begin_m
NAMELIST /time_control/ auxinput20_begin_s
NAMELIST /time_control/ auxinput20_begin
NAMELIST /time_control/ auxinput20_end_y
NAMELIST /time_control/ auxinput20_end_d
NAMELIST /time_control/ auxinput20_end_h
NAMELIST /time_control/ auxinput20_end_m
NAMELIST /time_control/ auxinput20_end_s
NAMELIST /time_control/ auxinput20_end
NAMELIST /time_control/ io_form_auxinput20
NAMELIST /time_control/ frames_per_auxinput20
NAMELIST /time_control/ auxinput21_inname
NAMELIST /time_control/ auxinput21_outname
NAMELIST /time_control/ auxinput21_interval_y
NAMELIST /time_control/ auxinput21_interval_d
NAMELIST /time_control/ auxinput21_interval_h
NAMELIST /time_control/ auxinput21_interval_m
NAMELIST /time_control/ auxinput21_interval_s
NAMELIST /time_control/ auxinput21_interval
NAMELIST /time_control/ auxinput21_begin_y
NAMELIST /time_control/ auxinput21_begin_d
NAMELIST /time_control/ auxinput21_begin_h
NAMELIST /time_control/ auxinput21_begin_m
NAMELIST /time_control/ auxinput21_begin_s
NAMELIST /time_control/ auxinput21_begin
NAMELIST /time_control/ auxinput21_end_y
NAMELIST /time_control/ auxinput21_end_d
NAMELIST /time_control/ auxinput21_end_h
NAMELIST /time_control/ auxinput21_end_m
NAMELIST /time_control/ auxinput21_end_s
NAMELIST /time_control/ auxinput21_end
NAMELIST /time_control/ io_form_auxinput21
NAMELIST /time_control/ frames_per_auxinput21
NAMELIST /time_control/ auxinput22_inname
NAMELIST /time_control/ auxinput22_outname
NAMELIST /time_control/ auxinput22_interval_y
NAMELIST /time_control/ auxinput22_interval_d
NAMELIST /time_control/ auxinput22_interval_h
NAMELIST /time_control/ auxinput22_interval_m
NAMELIST /time_control/ auxinput22_interval_s
NAMELIST /time_control/ auxinput22_interval
NAMELIST /time_control/ auxinput22_begin_y
NAMELIST /time_control/ auxinput22_begin_d
NAMELIST /time_control/ auxinput22_begin_h
NAMELIST /time_control/ auxinput22_begin_m
NAMELIST /time_control/ auxinput22_begin_s
NAMELIST /time_control/ auxinput22_begin
NAMELIST /time_control/ auxinput22_end_y
NAMELIST /time_control/ auxinput22_end_d
NAMELIST /time_control/ auxinput22_end_h
NAMELIST /time_control/ auxinput22_end_m
NAMELIST /time_control/ auxinput22_end_s
NAMELIST /time_control/ auxinput22_end
NAMELIST /time_control/ io_form_auxinput22
NAMELIST /time_control/ frames_per_auxinput22
NAMELIST /time_control/ auxinput23_inname
NAMELIST /time_control/ auxinput23_outname
NAMELIST /time_control/ auxinput23_interval_y
NAMELIST /time_control/ auxinput23_interval_d
NAMELIST /time_control/ auxinput23_interval_h
NAMELIST /time_control/ auxinput23_interval_m
NAMELIST /time_control/ auxinput23_interval_s
NAMELIST /time_control/ auxinput23_interval
NAMELIST /time_control/ auxinput23_begin_y
NAMELIST /time_control/ auxinput23_begin_d
NAMELIST /time_control/ auxinput23_begin_h
NAMELIST /time_control/ auxinput23_begin_m
NAMELIST /time_control/ auxinput23_begin_s
NAMELIST /time_control/ auxinput23_begin
NAMELIST /time_control/ auxinput23_end_y
NAMELIST /time_control/ auxinput23_end_d
NAMELIST /time_control/ auxinput23_end_h
NAMELIST /time_control/ auxinput23_end_m
NAMELIST /time_control/ auxinput23_end_s
NAMELIST /time_control/ auxinput23_end
NAMELIST /time_control/ io_form_auxinput23
NAMELIST /time_control/ frames_per_auxinput23
NAMELIST /time_control/ auxinput24_inname
NAMELIST /time_control/ auxinput24_outname
NAMELIST /time_control/ auxinput24_interval_y
NAMELIST /time_control/ auxinput24_interval_d
NAMELIST /time_control/ auxinput24_interval_h
NAMELIST /time_control/ auxinput24_interval_m
NAMELIST /time_control/ auxinput24_interval_s
NAMELIST /time_control/ auxinput24_interval
NAMELIST /time_control/ auxinput24_begin_y
NAMELIST /time_control/ auxinput24_begin_d
NAMELIST /time_control/ auxinput24_begin_h
NAMELIST /time_control/ auxinput24_begin_m
NAMELIST /time_control/ auxinput24_begin_s
NAMELIST /time_control/ auxinput24_begin
NAMELIST /time_control/ auxinput24_end_y
NAMELIST /time_control/ auxinput24_end_d
NAMELIST /time_control/ auxinput24_end_h
NAMELIST /time_control/ auxinput24_end_m
NAMELIST /time_control/ auxinput24_end_s
NAMELIST /time_control/ auxinput24_end
NAMELIST /time_control/ io_form_auxinput24
NAMELIST /time_control/ frames_per_auxinput24
NAMELIST /time_control/ history_interval
NAMELIST /time_control/ frames_per_outfile
NAMELIST /time_control/ restart
NAMELIST /time_control/ restart_interval
NAMELIST /time_control/ io_form_input
NAMELIST /time_control/ io_form_history
NAMELIST /time_control/ io_form_restart
NAMELIST /time_control/ io_form_boundary
NAMELIST /time_control/ debug_level
NAMELIST /time_control/ self_test_domain
NAMELIST /time_control/ history_outname
NAMELIST /time_control/ history_inname
NAMELIST /time_control/ use_netcdf_classic
NAMELIST /time_control/ history_interval_d
NAMELIST /time_control/ history_interval_h
NAMELIST /time_control/ history_interval_m
NAMELIST /time_control/ history_interval_s
NAMELIST /time_control/ inputout_interval_d
NAMELIST /time_control/ inputout_interval_h
NAMELIST /time_control/ inputout_interval_m
NAMELIST /time_control/ inputout_interval_s
NAMELIST /time_control/ inputout_interval
NAMELIST /time_control/ restart_interval_d
NAMELIST /time_control/ restart_interval_h
NAMELIST /time_control/ restart_interval_m
NAMELIST /time_control/ restart_interval_s
NAMELIST /time_control/ history_begin_y
NAMELIST /time_control/ history_begin_d
NAMELIST /time_control/ history_begin_h
NAMELIST /time_control/ history_begin_m
NAMELIST /time_control/ history_begin_s
NAMELIST /time_control/ history_begin
NAMELIST /time_control/ inputout_begin_y
NAMELIST /time_control/ inputout_begin_d
NAMELIST /time_control/ inputout_begin_h
NAMELIST /time_control/ inputout_begin_m
NAMELIST /time_control/ inputout_begin_s
NAMELIST /time_control/ restart_begin_y
NAMELIST /time_control/ restart_begin_d
NAMELIST /time_control/ restart_begin_h
NAMELIST /time_control/ restart_begin_m
NAMELIST /time_control/ restart_begin_s
NAMELIST /time_control/ restart_begin
NAMELIST /time_control/ history_end_y
NAMELIST /time_control/ history_end_d
NAMELIST /time_control/ history_end_h
NAMELIST /time_control/ history_end_m
NAMELIST /time_control/ history_end_s
NAMELIST /time_control/ history_end
NAMELIST /time_control/ inputout_end_y
NAMELIST /time_control/ inputout_end_d
NAMELIST /time_control/ inputout_end_h
NAMELIST /time_control/ inputout_end_m
NAMELIST /time_control/ inputout_end_s
NAMELIST /time_control/ reset_simulation_start
NAMELIST /domains/ sr_x
NAMELIST /domains/ sr_y
NAMELIST /fdda/ sgfdda_inname
NAMELIST /fdda/ gfdda_inname
NAMELIST /fdda/ sgfdda_interval_d
NAMELIST /fdda/ sgfdda_interval_h
NAMELIST /fdda/ sgfdda_interval_m
NAMELIST /fdda/ sgfdda_interval_s
NAMELIST /fdda/ sgfdda_interval_y
NAMELIST /fdda/ sgfdda_interval
NAMELIST /fdda/ gfdda_interval_d
NAMELIST /fdda/ gfdda_interval_h
NAMELIST /fdda/ gfdda_interval_m
NAMELIST /fdda/ gfdda_interval_s
NAMELIST /fdda/ gfdda_interval_y
NAMELIST /fdda/ gfdda_interval
NAMELIST /fdda/ sgfdda_begin_y
NAMELIST /fdda/ sgfdda_begin_d
NAMELIST /fdda/ sgfdda_begin_h
NAMELIST /fdda/ sgfdda_begin_m
NAMELIST /fdda/ sgfdda_begin_s
NAMELIST /fdda/ gfdda_begin_y
NAMELIST /fdda/ gfdda_begin_d
NAMELIST /fdda/ gfdda_begin_h
NAMELIST /fdda/ gfdda_begin_m
NAMELIST /fdda/ gfdda_begin_s
NAMELIST /fdda/ sgfdda_end_y
NAMELIST /fdda/ sgfdda_end_d
NAMELIST /fdda/ sgfdda_end_h
NAMELIST /fdda/ sgfdda_end_m
NAMELIST /fdda/ sgfdda_end_s
NAMELIST /fdda/ gfdda_end_y
NAMELIST /fdda/ gfdda_end_d
NAMELIST /fdda/ gfdda_end_h
NAMELIST /fdda/ gfdda_end_m
NAMELIST /fdda/ gfdda_end_s
NAMELIST /fdda/ io_form_sgfdda
NAMELIST /fdda/ io_form_gfdda
NAMELIST /time_control/ iofields_filename
NAMELIST /time_control/ ignore_iofields_warning
NAMELIST /time_control/ ncd_nofill
NAMELIST /wrfvar1/ update_sfcdiags
NAMELIST /wrfvar1/ use_wrf_sfcinfo
NAMELIST /wrfvar1/ use_background_errors
NAMELIST /wrfvar1/ write_increments
NAMELIST /wrfvar1/ var4d
NAMELIST /wrfvar1/ var4d_bin
NAMELIST /wrfvar1/ var4d_bin_rain
NAMELIST /wrfvar1/ var4d_lbc
NAMELIST /wrfvar1/ multi_inc
NAMELIST /wrfvar1/ print_detail_radar
NAMELIST /wrfvar1/ print_detail_rain
NAMELIST /wrfvar1/ print_detail_rad
NAMELIST /wrfvar1/ print_detail_xa
NAMELIST /wrfvar1/ print_detail_xb
NAMELIST /wrfvar1/ print_detail_obs
NAMELIST /wrfvar1/ print_detail_f_obs
NAMELIST /wrfvar1/ print_detail_map
NAMELIST /wrfvar1/ print_detail_grad
NAMELIST /wrfvar1/ print_detail_regression
NAMELIST /wrfvar1/ print_detail_spectral
NAMELIST /wrfvar1/ print_detail_testing
NAMELIST /wrfvar1/ print_detail_parallel
NAMELIST /wrfvar1/ print_detail_be
NAMELIST /wrfvar1/ print_detail_outerloop
NAMELIST /wrfvar1/ check_max_iv_print
NAMELIST /wrfvar1/ check_buddy_print
NAMELIST /wrfvar2/ analysis_accu
NAMELIST /wrfvar2/ calc_w_increment
NAMELIST /wrfvar2/ dt_cloud_model
NAMELIST /wrfvar2/ write_mod_filtered_obs
NAMELIST /wrfvar2/ wind_sd
NAMELIST /wrfvar2/ wind_sd_buoy
NAMELIST /wrfvar2/ wind_sd_synop
NAMELIST /wrfvar2/ wind_sd_ships
NAMELIST /wrfvar2/ wind_sd_metar
NAMELIST /wrfvar2/ wind_sd_sound
NAMELIST /wrfvar2/ wind_sd_pilot
NAMELIST /wrfvar2/ wind_sd_airep
NAMELIST /wrfvar2/ wind_sd_qscat
NAMELIST /wrfvar2/ wind_sd_tamdar
NAMELIST /wrfvar2/ wind_sd_geoamv
NAMELIST /wrfvar2/ wind_sd_mtgirs
NAMELIST /wrfvar2/ wind_sd_polaramv
NAMELIST /wrfvar2/ wind_sd_profiler
NAMELIST /wrfvar2/ wind_stats_sd
NAMELIST /wrfvar2/ qc_rej_both
NAMELIST /wrfvar3/ fg_format
NAMELIST /wrfvar3/ ob_format
NAMELIST /wrfvar3/ ob_format_gpsro
NAMELIST /wrfvar3/ num_fgat_time
NAMELIST /wrfvar4/ thin_conv
NAMELIST /wrfvar4/ thin_conv_ascii
NAMELIST /wrfvar4/ thin_mesh_conv
NAMELIST /wrfvar4/ thin_rainobs
NAMELIST /wrfvar4/ use_synopobs
NAMELIST /wrfvar4/ use_shipsobs
NAMELIST /wrfvar4/ use_metarobs
NAMELIST /wrfvar4/ use_soundobs
NAMELIST /wrfvar4/ use_mtgirsobs
NAMELIST /wrfvar4/ use_tamdarobs
NAMELIST /wrfvar4/ use_pilotobs
NAMELIST /wrfvar4/ use_airepobs
NAMELIST /wrfvar4/ use_geoamvobs
NAMELIST /wrfvar4/ use_polaramvobs
NAMELIST /wrfvar4/ use_bogusobs
NAMELIST /wrfvar4/ use_buoyobs
NAMELIST /wrfvar4/ use_profilerobs
NAMELIST /wrfvar4/ use_satemobs
NAMELIST /wrfvar4/ use_gpsztdobs
NAMELIST /wrfvar4/ use_gpspwobs
NAMELIST /wrfvar4/ use_gpsrefobs
NAMELIST /wrfvar4/ top_km_gpsro
NAMELIST /wrfvar4/ bot_km_gpsro
NAMELIST /wrfvar4/ use_ssmiretrievalobs
NAMELIST /wrfvar4/ use_ssmitbobs
NAMELIST /wrfvar4/ use_ssmt1obs
NAMELIST /wrfvar4/ use_ssmt2obs
NAMELIST /wrfvar4/ use_qscatobs
NAMELIST /wrfvar4/ use_radarobs
NAMELIST /wrfvar4/ use_radar_rv
NAMELIST /wrfvar4/ use_radar_rf
NAMELIST /wrfvar4/ use_radar_rqv
NAMELIST /wrfvar4/ use_radar_rhv
NAMELIST /wrfvar4/ use_3dvar_phy
NAMELIST /wrfvar4/ use_rainobs
NAMELIST /wrfvar4/ use_hirs2obs
NAMELIST /wrfvar4/ use_hirs3obs
NAMELIST /wrfvar4/ use_hirs4obs
NAMELIST /wrfvar4/ use_mhsobs
NAMELIST /wrfvar4/ use_msuobs
NAMELIST /wrfvar4/ use_amsuaobs
NAMELIST /wrfvar4/ use_amsubobs
NAMELIST /wrfvar4/ use_airsobs
NAMELIST /wrfvar4/ use_airsretobs
NAMELIST /wrfvar4/ use_eos_amsuaobs
NAMELIST /wrfvar4/ use_hsbobs
NAMELIST /wrfvar4/ use_ssmisobs
NAMELIST /wrfvar4/ use_iasiobs
NAMELIST /wrfvar4/ use_seviriobs
NAMELIST /wrfvar4/ use_amsr2obs
NAMELIST /wrfvar4/ use_kma1dvar
NAMELIST /wrfvar4/ use_filtered_rad
NAMELIST /wrfvar4/ use_obs_errfac
NAMELIST /wrfvar4/ use_atmsobs
NAMELIST /wrfvar4/ use_mwtsobs
NAMELIST /wrfvar4/ use_mwhsobs
NAMELIST /wrfvar5/ check_max_iv
NAMELIST /wrfvar5/ max_error_t
NAMELIST /wrfvar5/ max_error_uv
NAMELIST /wrfvar5/ max_error_spd
NAMELIST /wrfvar5/ max_error_dir
NAMELIST /wrfvar5/ max_omb_spd
NAMELIST /wrfvar5/ max_omb_dir
NAMELIST /wrfvar5/ max_error_pw
NAMELIST /wrfvar5/ max_error_ref
NAMELIST /wrfvar5/ max_error_rh
NAMELIST /wrfvar5/ max_error_q
NAMELIST /wrfvar5/ max_error_p
NAMELIST /wrfvar5/ max_error_tb
NAMELIST /wrfvar5/ max_error_thickness
NAMELIST /wrfvar5/ max_error_rv
NAMELIST /wrfvar5/ max_error_rf
NAMELIST /wrfvar5/ max_error_rain
NAMELIST /wrfvar5/ max_error_buv
NAMELIST /wrfvar5/ max_error_bt
NAMELIST /wrfvar5/ max_error_bq
NAMELIST /wrfvar5/ max_error_slp
NAMELIST /wrfvar5/ check_buddy
NAMELIST /wrfvar5/ put_rand_seed
NAMELIST /wrfvar5/ omb_set_rand
NAMELIST /wrfvar5/ omb_add_noise
NAMELIST /wrfvar5/ position_lev_dependant
NAMELIST /wrfvar5/ obs_qc_pointer
NAMELIST /wrfvar5/ qmarker_retain
NAMELIST /wrfvar5/ max_sound_input
NAMELIST /wrfvar5/ max_mtgirs_input
NAMELIST /wrfvar5/ max_tamdar_input
NAMELIST /wrfvar5/ max_synop_input
NAMELIST /wrfvar5/ max_geoamv_input
NAMELIST /wrfvar5/ max_polaramv_input
NAMELIST /wrfvar5/ max_airep_input
NAMELIST /wrfvar5/ max_satem_input
NAMELIST /wrfvar5/ max_pilot_input
NAMELIST /wrfvar5/ max_radar_input
NAMELIST /wrfvar5/ max_rain_input
NAMELIST /wrfvar5/ max_metar_input
NAMELIST /wrfvar5/ max_gpspw_input
NAMELIST /wrfvar5/ max_ships_input
NAMELIST /wrfvar5/ max_profiler_input
NAMELIST /wrfvar5/ max_bogus_input
NAMELIST /wrfvar5/ max_buoy_input
NAMELIST /wrfvar5/ max_ssmi_rv_input
NAMELIST /wrfvar5/ max_ssmi_tb_input
NAMELIST /wrfvar5/ max_ssmt1_input
NAMELIST /wrfvar5/ max_ssmt2_input
NAMELIST /wrfvar5/ max_qscat_input
NAMELIST /wrfvar5/ max_gpsref_input
NAMELIST /wrfvar5/ max_airsr_input
NAMELIST /wrfvar5/ max_tovs_input
NAMELIST /wrfvar5/ max_ssmis_input
NAMELIST /wrfvar5/ report_start
NAMELIST /wrfvar5/ report_end
NAMELIST /wrfvar5/ tovs_start
NAMELIST /wrfvar5/ tovs_end
NAMELIST /wrfvar5/ gpsref_thinning
NAMELIST /wrfvar6/ outer_loop_restart
NAMELIST /wrfvar6/ max_ext_its
NAMELIST /wrfvar6/ ntmax
NAMELIST /wrfvar6/ nsave
NAMELIST /wrfvar6/ write_interval
NAMELIST /wrfvar6/ eps
NAMELIST /wrfvar6/ precondition_cg
NAMELIST /wrfvar6/ precondition_factor
NAMELIST /wrfvar6/ use_lanczos
NAMELIST /wrfvar6/ read_lanczos
NAMELIST /wrfvar6/ write_lanczos
NAMELIST /wrfvar6/ orthonorm_gradient
NAMELIST /wrfvar7/ cv_options
NAMELIST /wrfvar7/ cloud_cv_options
NAMELIST /wrfvar7/ as1
NAMELIST /wrfvar7/ as2
NAMELIST /wrfvar7/ as3
NAMELIST /wrfvar7/ as4
NAMELIST /wrfvar7/ as5
NAMELIST /wrfvar7/ do_normalize
NAMELIST /wrfvar7/ use_rf
NAMELIST /wrfvar7/ rf_passes
NAMELIST /wrfvar7/ var_scaling1
NAMELIST /wrfvar7/ var_scaling2
NAMELIST /wrfvar7/ var_scaling3
NAMELIST /wrfvar7/ var_scaling4
NAMELIST /wrfvar7/ var_scaling5
NAMELIST /wrfvar7/ var_scaling6
NAMELIST /wrfvar7/ var_scaling7
NAMELIST /wrfvar7/ var_scaling8
NAMELIST /wrfvar7/ var_scaling9
NAMELIST /wrfvar7/ var_scaling10
NAMELIST /wrfvar7/ var_scaling11
NAMELIST /wrfvar7/ len_scaling1
NAMELIST /wrfvar7/ len_scaling2
NAMELIST /wrfvar7/ len_scaling3
NAMELIST /wrfvar7/ len_scaling4
NAMELIST /wrfvar7/ len_scaling5
NAMELIST /wrfvar7/ len_scaling6
NAMELIST /wrfvar7/ len_scaling7
NAMELIST /wrfvar7/ len_scaling8
NAMELIST /wrfvar7/ len_scaling9
NAMELIST /wrfvar7/ len_scaling10
NAMELIST /wrfvar7/ len_scaling11
NAMELIST /wrfvar7/ je_factor
NAMELIST /wrfvar7/ power_truncation
NAMELIST /wrfvar8/ def_sub_domain
NAMELIST /wrfvar8/ x_start_sub_domain
NAMELIST /wrfvar8/ y_start_sub_domain
NAMELIST /wrfvar8/ x_end_sub_domain
NAMELIST /wrfvar8/ y_end_sub_domain
NAMELIST /wrfvar9/ stdout
NAMELIST /wrfvar9/ stderr
NAMELIST /wrfvar9/ trace_unit
NAMELIST /wrfvar9/ trace_pe
NAMELIST /wrfvar9/ trace_repeat_head
NAMELIST /wrfvar9/ trace_repeat_body
NAMELIST /wrfvar9/ trace_max_depth
NAMELIST /wrfvar9/ trace_use
NAMELIST /wrfvar9/ trace_use_frequent
NAMELIST /wrfvar9/ trace_use_dull
NAMELIST /wrfvar9/ trace_memory
NAMELIST /wrfvar9/ trace_all_pes
NAMELIST /wrfvar9/ trace_csv
NAMELIST /wrfvar9/ use_html
NAMELIST /wrfvar9/ warnings_are_fatal
NAMELIST /wrfvar10/ test_transforms
NAMELIST /wrfvar10/ test_gradient
NAMELIST /wrfvar10/ test_statistics
NAMELIST /wrfvar10/ interpolate_stats
NAMELIST /wrfvar10/ be_eta
NAMELIST /wrfvar10/ test_dm_exact
NAMELIST /wrfvar11/ cv_options_hum
NAMELIST /wrfvar11/ check_rh
NAMELIST /wrfvar11/ set_omb_rand_fac
NAMELIST /wrfvar11/ seed_array1
NAMELIST /wrfvar11/ seed_array2
NAMELIST /wrfvar11/ sfc_assi_options
NAMELIST /wrfvar11/ psfc_from_slp
NAMELIST /wrfvar11/ calculate_cg_cost_fn
NAMELIST /wrfvar11/ lat_stats_option
NAMELIST /wrfvar11/ interp_option
NAMELIST /wrfvar12/ balance_type
NAMELIST /wrfvar12/ use_wpec
NAMELIST /wrfvar12/ wpec_factor
NAMELIST /wrfvar13/ vert_corr
NAMELIST /wrfvar13/ vertical_ip
NAMELIST /wrfvar13/ vert_evalue
NAMELIST /wrfvar13/ max_vert_var1
NAMELIST /wrfvar13/ max_vert_var2
NAMELIST /wrfvar13/ max_vert_var3
NAMELIST /wrfvar13/ max_vert_var4
NAMELIST /wrfvar13/ max_vert_var5
NAMELIST /wrfvar13/ max_vert_var6
NAMELIST /wrfvar13/ max_vert_var7
NAMELIST /wrfvar13/ max_vert_var8
NAMELIST /wrfvar13/ max_vert_var9
NAMELIST /wrfvar13/ max_vert_var10
NAMELIST /wrfvar13/ max_vert_var11
NAMELIST /wrfvar13/ max_vert_var_alpha
NAMELIST /wrfvar13/ psi_chi_factor
NAMELIST /wrfvar13/ psi_t_factor
NAMELIST /wrfvar13/ psi_ps_factor
NAMELIST /wrfvar13/ psi_rh_factor
NAMELIST /wrfvar13/ chi_u_t_factor
NAMELIST /wrfvar13/ chi_u_ps_factor
NAMELIST /wrfvar13/ chi_u_rh_factor
NAMELIST /wrfvar13/ t_u_rh_factor
NAMELIST /wrfvar13/ ps_u_rh_factor
NAMELIST /wrfvar14/ rttov_emis_atlas_ir
NAMELIST /wrfvar14/ rttov_emis_atlas_mw
NAMELIST /wrfvar14/ rtminit_print
NAMELIST /wrfvar14/ rtminit_nsensor
NAMELIST /wrfvar14/ rtminit_platform
NAMELIST /wrfvar14/ rtminit_satid
NAMELIST /wrfvar14/ rtminit_sensor
NAMELIST /wrfvar14/ rad_monitoring
NAMELIST /wrfvar14/ thinning_mesh
NAMELIST /wrfvar14/ thinning
NAMELIST /wrfvar14/ read_biascoef
NAMELIST /wrfvar14/ biascorr
NAMELIST /wrfvar14/ biasprep
NAMELIST /wrfvar14/ rttov_scatt
NAMELIST /wrfvar14/ write_profile
NAMELIST /wrfvar14/ write_jacobian
NAMELIST /wrfvar14/ qc_rad
NAMELIST /wrfvar14/ write_iv_rad_ascii
NAMELIST /wrfvar14/ write_oa_rad_ascii
NAMELIST /wrfvar14/ write_filtered_rad
NAMELIST /wrfvar14/ use_error_factor_rad
NAMELIST /wrfvar14/ use_landem
NAMELIST /wrfvar14/ use_antcorr
NAMELIST /wrfvar14/ use_mspps_emis
NAMELIST /wrfvar14/ use_mspps_ts
NAMELIST /wrfvar14/ mw_emis_sea
NAMELIST /wrfvar14/ tovs_min_transfer
NAMELIST /wrfvar14/ tovs_batch
NAMELIST /wrfvar14/ rtm_option
NAMELIST /wrfvar14/ use_crtm_kmatrix
NAMELIST /wrfvar14/ use_rttov_kmatrix
NAMELIST /wrfvar14/ crtm_cloud
NAMELIST /wrfvar14/ only_sea_rad
NAMELIST /wrfvar14/ use_pseudo_rad
NAMELIST /wrfvar14/ pseudo_rad_platid
NAMELIST /wrfvar14/ pseudo_rad_satid
NAMELIST /wrfvar14/ pseudo_rad_senid
NAMELIST /wrfvar14/ pseudo_rad_ichan
NAMELIST /wrfvar14/ pseudo_rad_lat
NAMELIST /wrfvar14/ pseudo_rad_lon
NAMELIST /wrfvar14/ pseudo_rad_inv
NAMELIST /wrfvar14/ pseudo_rad_err
NAMELIST /wrfvar14/ use_simulated_rad
NAMELIST /wrfvar14/ simulated_rad_io
NAMELIST /wrfvar14/ simulated_rad_ngrid
NAMELIST /wrfvar14/ use_varbc
NAMELIST /wrfvar14/ freeze_varbc
NAMELIST /wrfvar14/ varbc_factor
NAMELIST /wrfvar14/ varbc_nbgerr
NAMELIST /wrfvar14/ varbc_nobsmin
NAMELIST /wrfvar14/ use_clddet_mmr
NAMELIST /wrfvar14/ use_clddet_ecmwf
NAMELIST /wrfvar14/ airs_warmest_fov
NAMELIST /wrfvar14/ use_satcv
NAMELIST /wrfvar14/ use_blacklist_rad
NAMELIST /wrfvar14/ calc_weightfunc
NAMELIST /wrfvar14/ crtm_coef_path
NAMELIST /wrfvar14/ crtm_irwater_coef
NAMELIST /wrfvar14/ crtm_mwwater_coef
NAMELIST /wrfvar14/ crtm_irland_coef
NAMELIST /wrfvar14/ crtm_visland_coef
NAMELIST /wrfvar15/ num_pseudo
NAMELIST /wrfvar15/ pseudo_x
NAMELIST /wrfvar15/ pseudo_y
NAMELIST /wrfvar15/ pseudo_z
NAMELIST /wrfvar15/ pseudo_val
NAMELIST /wrfvar15/ pseudo_err
NAMELIST /wrfvar16/ alphacv_method
NAMELIST /wrfvar16/ ensdim_alpha
NAMELIST /wrfvar16/ alpha_truncation
NAMELIST /wrfvar16/ alpha_corr_type
NAMELIST /wrfvar16/ alpha_corr_scale
NAMELIST /wrfvar16/ alpha_std_dev
NAMELIST /wrfvar16/ alpha_vertloc
NAMELIST /wrfvar16/ alpha_hydrometeors
NAMELIST /wrfvar16/ hybrid_dual_res
NAMELIST /wrfvar16/ dual_res_upscale_opt
NAMELIST /wrfvar17/ analysis_type
NAMELIST /wrfvar17/ sensitivity_option
NAMELIST /wrfvar17/ adj_sens
NAMELIST /wrfvar18/ analysis_date
NAMELIST /wrfvar19/ pseudo_var
NAMELIST /wrfvar20/ documentation_url
NAMELIST /wrfvar21/ time_window_min
NAMELIST /wrfvar22/ time_window_max
NAMELIST /perturbation/ jcdfi_use
NAMELIST /perturbation/ jcdfi_diag
NAMELIST /perturbation/ jcdfi_penalty
NAMELIST /perturbation/ enable_identity
NAMELIST /perturbation/ trajectory_io
NAMELIST /perturbation/ var4d_detail_out
NAMELIST /perturbation/ var4d_run
NAMELIST /physics/ mp_physics_ad
NAMELIST /physics/ chem_opt


      esmf_coupling = .FALSE.

      quilting_enabled = .FALSE.
      IF ( disable_quilt ) RETURN

      DO i = 1,int_num_handles
        okay_to_write(i) = .FALSE.
        int_handle_in_use(i) = .FALSE.
        server_for_handle(i) = 0 
        int_num_bytes_to_write(i) = 0
      ENDDO

      CALL MPI_INITIALIZED( mpi_inited, ierr )
      IF ( .NOT. mpi_inited ) THEN
        CALL wrf_error_fatal3("<stdin>",6433,&
"module_io_quilt_old.F : MPI not init'd" )
      ENDIF
      CALL wrf_get_dm_quilt_comm( mpi_comm_here )   

      CALL MPI_Comm_rank( mpi_comm_here, mytask, ierr ) ;
      CALL Mpi_Comm_Size (  mpi_comm_here, ntasks, ierr  ) ;

      IF ( mytask .EQ. 0 ) THEN
        OPEN ( unit=27, file="namelist.input", form="formatted", status="old" )
        nio_groups = 1
        nio_tasks_per_group  = 0
        poll_servers = .false.
        READ ( 27 , NML = namelist_quilt, IOSTAT=io_status )
        IF (io_status .NE. 0) THEN
          CALL wrf_error_fatal3("<stdin>",6448,&
"ERROR reading namelist namelist_quilt" )
        ENDIF
        REWIND(27)
        nproc_x = -1
        nproc_y = -1
        READ ( UNIT = 27 , NML = domains , IOSTAT=io_status )
        IF (io_status .NE. 0) THEN
          CALL wrf_error_fatal3("<stdin>",6456,&
"ERROR reading namelist domains" )
        ENDIF
        CLOSE ( 27 )
        IF ( esmf_coupling ) THEN
          IF ( nio_tasks_per_group > 0 ) THEN
            CALL wrf_error_fatal3("<stdin>",6462,&
"frame/module_io_quilt.F: cannot use "// &
                                 "ESMF coupling with quilt tasks") ;
          ENDIF
        ENDIF
        if(poll_servers) then
           temp_poll=1
        else
           temp_poll=0
        endif
      ENDIF

      CALL mpi_bcast( nio_tasks_per_group  , 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nio_groups , 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( temp_poll , 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nproc_x , 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nproc_y , 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )

      poll_servers = (temp_poll == 1)

      CALL setup_quilt_servers( nio_tasks_per_group,            &
                                mytask,               &
                                ntasks,               &
                                nproc_x,              &
                                nproc_y,              &
                                nio_groups,           &
                                nio_tasks_in_group,   &
                                mpi_comm_here,        &
                                mpi_comm_local,       &
                                mpi_comm_io_groups)

      call init_module_wrf_error(on_io_server=.true.)

       
       IF ( compute_node ) THEN
          mpi_comm_allcompute = mpi_comm_local
          CALL wrf_set_dm_communicator( mpi_comm_local )
       ELSE
          CALL quilt    
       ENDIF
      RETURN
    END SUBROUTINE init_module_wrf_quilt




END MODULE module_wrf_quilt







SUBROUTINE disable_quilting




  USE module_wrf_quilt
  disable_quilt = .TRUE.
  RETURN
END SUBROUTINE disable_quilting

SUBROUTINE quilting_disabled( reslt )




  USE module_wrf_quilt
  LOGICAL, INTENT(OUT) :: reslt
  reslt = disable_quilt
write(0,*)"module_io_quilt_old.F",2931,disable_quilt
  RETURN
END SUBROUTINE quilting_disabled

LOGICAL FUNCTION  use_output_servers_for(ioform)







  USE module_wrf_quilt
  integer, intent(in) :: ioform
  use_output_servers_for = quilting_enabled
  use_output_servers_for = ( use_output_servers_for .and. ioform<100 )
  RETURN
END FUNCTION use_output_servers_for

LOGICAL FUNCTION  use_output_servers()




  USE module_wrf_quilt
  use_output_servers = quilting_enabled
  RETURN
END FUNCTION use_output_servers

LOGICAL FUNCTION  use_input_servers()




  USE module_wrf_quilt
  use_input_servers = .FALSE.
  RETURN
END FUNCTION use_input_servers

SUBROUTINE wrf_quilt_open_for_write_begin( FileName , gridid, Comm_compute, Comm_io, SysDepInfo, &
                                     DataHandle , io_form_arg, Status )





  USE module_wrf_quilt
  USE module_state_description, ONLY: IO_PNETCDF
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  CHARACTER *(*), INTENT(IN)  :: FileName
  INTEGER ,       INTENT(IN)  :: gridid
  INTEGER ,       INTENT(IN)  :: Comm_compute , Comm_io
  CHARACTER *(*), INTENT(IN)  :: SysDepInfo
  INTEGER ,       INTENT(OUT) :: DataHandle
  INTEGER ,       INTENT(IN)  :: io_form_arg
  INTEGER ,       INTENT(OUT) :: Status

  CHARACTER*132   :: locFileName, locSysDepInfo
  INTEGER i, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy
  INTEGER, EXTERNAL :: use_package

  CALL wrf_debug ( 50, 'in wrf_quilt_open_for_write_begin' ) 
  CALL int_get_fresh_handle(i)
  okay_to_write(i) = .false.
  DataHandle = i

  locFileName = FileName
  locSysDepInfo = SysDepInfo

  CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )

  SELECT CASE(use_package(io_form_arg))

  CASE DEFAULT

     IF ( wrf_dm_on_monitor() ) THEN
        CALL int_gen_ofwb_header( hdrbuf, hdrbufsize, itypesize, &
                                  locFileName,locSysDepInfo,io_form_arg,DataHandle )
     ELSE
        CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
     ENDIF

  END SELECT

  iserver = get_server_id ( DataHandle )
  CALL get_mpi_comm_io_groups( comm_io_group , iserver )

  CALL Mpi_Comm_Size (  comm_io_group , tasks_in_group , ierr  )


  
  reduced = 0
  reduced(1) = hdrbufsize 
  IF ( wrf_dm_on_monitor() )  reduced(2) = i 
  CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )


  
  CALL collect_on_comm_debug("module_io_quilt_old.F",3047, comm_io_group,            &
                        onebyte,                       &
                        hdrbuf, hdrbufsize , &
                        dummy, 0 )

  Status = 0


  RETURN  
END SUBROUTINE wrf_quilt_open_for_write_begin

SUBROUTINE wrf_quilt_open_for_write_commit( DataHandle , Status )







  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN ) :: DataHandle
  INTEGER ,       INTENT(OUT) :: Status
  INTEGER i, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy

  CALL wrf_debug ( 50, 'in wrf_quilt_open_for_write_commit' ) 
  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      okay_to_write( DataHandle ) = .true.
    ENDIF
  ENDIF

  CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )


  IF ( wrf_dm_on_monitor() ) THEN
     CALL int_gen_handle_header( hdrbuf, hdrbufsize, itypesize, &
                                 DataHandle, int_open_for_write_commit )
  ELSE
     CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
  ENDIF

  iserver = get_server_id ( DataHandle )
  CALL get_mpi_comm_io_groups( comm_io_group , iserver )

  CALL Mpi_Comm_Size (  comm_io_group , tasks_in_group , ierr  )


  
  reduced = 0
  reduced(1) = hdrbufsize 
  IF ( wrf_dm_on_monitor() )  reduced(2) = DataHandle
  CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )


  
  CALL collect_on_comm_debug("module_io_quilt_old.F",3123, comm_io_group,            &
                        onebyte,                       &
                        hdrbuf, hdrbufsize , &
                        dummy, 0 )

  Status = 0

  RETURN  
END SUBROUTINE wrf_quilt_open_for_write_commit

SUBROUTINE wrf_quilt_open_for_read ( FileName , Comm_compute, Comm_io, SysDepInfo, &
                               DataHandle , Status )





  IMPLICIT NONE
  CHARACTER *(*), INTENT(IN)  :: FileName
  INTEGER ,       INTENT(IN)  :: Comm_compute , Comm_io
  CHARACTER *(*), INTENT(IN)  :: SysDepInfo
  INTEGER ,       INTENT(OUT) :: DataHandle
  INTEGER ,       INTENT(OUT) :: Status

  CALL wrf_debug ( 50, 'in wrf_quilt_open_for_read' ) 
  DataHandle = -1
  Status = -1
  CALL wrf_error_fatal3("<stdin>",6787,&
"frame/module_io_quilt.F: wrf_quilt_open_for_read not yet supported" )
  RETURN  
END SUBROUTINE wrf_quilt_open_for_read

SUBROUTINE wrf_quilt_inquire_opened ( DataHandle, FileName , FileStatus, Status )





  USE module_wrf_quilt
  IMPLICIT NONE
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101
      integer, parameter  :: WRF_FILE_OPENED_FOR_WRITE            = 102
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103
      integer, parameter  :: WRF_REAL                             = 105
      integer, parameter  :: WRF_DOUBLE                           = 105
      integer, parameter  :: WRF_FLOAT=WRF_REAL
      integer, parameter  :: WRF_INTEGER                          = 106
      integer, parameter  :: WRF_LOGICAL                          = 107
      integer, parameter  :: WRF_COMPLEX                          = 108
      integer, parameter  :: WRF_DOUBLE_COMPLEX                   = 109
      integer, parameter  :: WRF_FILE_OPENED_FOR_UPDATE           = 110


      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER *(*), INTENT(IN)  :: FileName
  INTEGER ,       INTENT(OUT) :: FileStatus
  INTEGER ,       INTENT(OUT) :: Status

  Status = 0

  CALL wrf_debug ( 50, 'in wrf_quilt_inquire_opened' ) 
  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      IF ( okay_to_write( DataHandle ) ) THEN
        FileStatus = WRF_FILE_OPENED_FOR_WRITE
      ENDIF
    ENDIF
  ENDIF
  Status = 0
  
  RETURN
END SUBROUTINE wrf_quilt_inquire_opened

SUBROUTINE wrf_quilt_inquire_filename ( DataHandle, FileName , FileStatus, Status )










  USE module_wrf_quilt
  IMPLICIT NONE
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101
      integer, parameter  :: WRF_FILE_OPENED_FOR_WRITE            = 102
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103
      integer, parameter  :: WRF_REAL                             = 105
      integer, parameter  :: WRF_DOUBLE                           = 105
      integer, parameter  :: WRF_FLOAT=WRF_REAL
      integer, parameter  :: WRF_INTEGER                          = 106
      integer, parameter  :: WRF_LOGICAL                          = 107
      integer, parameter  :: WRF_COMPLEX                          = 108
      integer, parameter  :: WRF_DOUBLE_COMPLEX                   = 109
      integer, parameter  :: WRF_FILE_OPENED_FOR_UPDATE           = 110


      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER *(*), INTENT(OUT) :: FileName
  INTEGER ,       INTENT(OUT) :: FileStatus
  INTEGER ,       INTENT(OUT) :: Status
  CALL wrf_debug ( 50, 'in wrf_quilt_inquire_filename' ) 
  Status = 0
  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      IF ( okay_to_write( DataHandle ) ) THEN
        FileStatus = WRF_FILE_OPENED_FOR_WRITE
      ELSE
        FileStatus = WRF_FILE_OPENED_NOT_COMMITTED
      ENDIF
    ELSE
        FileStatus = WRF_FILE_NOT_OPENED
    ENDIF
    Status = 0
    FileName = "bogusfornow"
  ELSE
    Status = -1
  ENDIF
  RETURN
END SUBROUTINE wrf_quilt_inquire_filename

SUBROUTINE wrf_quilt_iosync ( DataHandle, Status )



















  USE module_wrf_quilt
  IMPLICIT NONE
  include "mpif.h"
  INTEGER ,       INTENT(IN)  :: DataHandle
  INTEGER ,       INTENT(OUT) :: Status

  INTEGER locsize , itypesize
  INTEGER ierr, tasks_in_group, comm_io_group, dummy, i

  CALL wrf_debug ( 50, 'in wrf_quilt_iosync' ) 


  IF ( associated ( int_local_output_buffer ) ) THEN

    iserver = get_server_id ( DataHandle )
    CALL get_mpi_comm_io_groups( comm_io_group , iserver )

    CALL Mpi_Comm_Size (  comm_io_group , tasks_in_group , ierr  )

    locsize = int_num_bytes_to_write(DataHandle)


    
    reduced = 0
    reduced(1) = locsize 
    IF ( wrf_dm_on_monitor() )  reduced(2) = DataHandle
    CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )


    
    CALL collect_on_comm_debug("module_io_quilt_old.F",3290, comm_io_group,            &
                          onebyte,                       &
                          int_local_output_buffer, locsize , &
                          dummy, 0 )


    int_local_output_cursor = 1

    DEALLOCATE ( int_local_output_buffer )
    NULLIFY ( int_local_output_buffer )
  ELSE
    CALL wrf_message ("frame/module_io_quilt.F: wrf_quilt_iosync: no buffer allocated")
  ENDIF

  Status = 0
  RETURN
END SUBROUTINE wrf_quilt_iosync

SUBROUTINE wrf_quilt_ioclose ( DataHandle, Status )







  USE module_wrf_quilt
  USE module_timing
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  INTEGER ,       INTENT(OUT) :: Status
  INTEGER i, itypesize, tasks_in_group, comm_io_group, ierr
  REAL dummy


  CALL wrf_debug ( 50, 'in wrf_quilt_ioclose' ) 
  CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )




  IF ( wrf_dm_on_monitor() ) THEN
     CALL int_gen_handle_header( hdrbuf, hdrbufsize, itypesize, &
                                 DataHandle , int_ioclose )
  ELSE
     CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
  ENDIF

  iserver = get_server_id ( DataHandle )
  CALL get_mpi_comm_io_groups( comm_io_group , iserver )

  CALL Mpi_Comm_Size (  comm_io_group , tasks_in_group , ierr  )


  
  reduced = 0
  IF ( wrf_dm_on_monitor() )  reduced(2) = DataHandle
  CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )
!!JMTIMING   CALL end_timing("MPI_Reduce in ioclose")


  int_handle_in_use(DataHandle) = .false.
  CALL set_server_id( DataHandle, 0 ) 
  okay_to_write(DataHandle) = .false.
  okay_to_commit(DataHandle) = .false.
  int_local_output_cursor = 1
  int_num_bytes_to_write(DataHandle) = 0
  IF ( associated ( int_local_output_buffer ) ) THEN
    DEALLOCATE ( int_local_output_buffer )
    NULLIFY ( int_local_output_buffer )
  ENDIF

  Status = 0
!!JMTIMING   CALL end_timing( "wrf_quilt_ioclose" )

  RETURN
END SUBROUTINE wrf_quilt_ioclose

SUBROUTINE wrf_quilt_ioexit( Status )





  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(OUT) :: Status
  INTEGER                     :: DataHandle, actual_iserver
  INTEGER i, itypesize, tasks_in_group, comm_io_group, me, ierr 
  REAL dummy

  CALL wrf_debug ( 50, 'in wrf_quilt_ioexit' ) 
  CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )




  IF ( wrf_dm_on_monitor() ) THEN
     CALL int_gen_handle_header( hdrbuf, hdrbufsize, itypesize, &
                                 DataHandle , int_ioexit )  
  ELSE
     CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
  ENDIF

  DO iserver = 1, nio_groups
    if(poll_servers) then
       
       
       

       call wrf_quilt_find_server(actual_iserver)

       
       
       
    else
       
       actual_iserver=iserver
    endif
    CALL get_mpi_comm_io_groups( comm_io_group , actual_iserver )

    CALL Mpi_Comm_Size (  comm_io_group , tasks_in_group , ierr  )
    CALL mpi_comm_rank( comm_io_group , me , ierr )


    hdrbufsize = -100 
    reduced = 0
    IF ( me .eq. 0 ) reduced(1) = hdrbufsize 
    CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )

  ENDDO
  Status = 0

  RETURN  
END SUBROUTINE wrf_quilt_ioexit

SUBROUTINE wrf_quilt_get_next_time ( DataHandle, DateStr, Status )





  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*)               :: DateStr
  INTEGER                     :: Status
  RETURN
END SUBROUTINE wrf_quilt_get_next_time

SUBROUTINE wrf_quilt_get_previous_time ( DataHandle, DateStr, Status )





  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*)               :: DateStr
  INTEGER                     :: Status
  RETURN
END SUBROUTINE wrf_quilt_get_previous_time

SUBROUTINE wrf_quilt_set_time ( DataHandle, Data,  Status )





  USE module_wrf_quilt
  USE module_state_description, ONLY: IO_PNETCDF
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Data
  INTEGER                     :: Status
  INTEGER i, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy
  INTEGER                 :: Count
  INTEGER, EXTERNAL       :: use_package

  CALL wrf_debug ( 50, 'in wrf_quilt_set_time' )

  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
      Count = 0   



      IF ( wrf_dm_on_monitor() ) THEN
         CALL int_gen_ti_header_char( hdrbuf, hdrbufsize, itypesize, &
                                      DataHandle, "TIMESTAMP", "", Data, int_set_time )
      ELSE
         CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
      ENDIF

      iserver = get_server_id ( DataHandle )
      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      CALL Mpi_Comm_Size (  comm_io_group , tasks_in_group , ierr  )

      
      reduced = 0
      reduced(1) = hdrbufsize 
      IF ( wrf_dm_on_monitor() )  reduced(2) = DataHandle
      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )
      
      CALL collect_on_comm_debug("module_io_quilt_old.F",3562, comm_io_group,            &
                            onebyte,                       &
                            hdrbuf, hdrbufsize , &
                            dummy, 0 )
    ENDIF
  ENDIF

RETURN
END SUBROUTINE wrf_quilt_set_time

SUBROUTINE wrf_quilt_get_next_var ( DataHandle, VarName, Status )






  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*)               :: VarName
  INTEGER                     :: Status
  RETURN
END SUBROUTINE wrf_quilt_get_next_var

SUBROUTINE wrf_quilt_get_dom_ti_real ( DataHandle,Element,   Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  REAL,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Outcount
  INTEGER                     :: Status
  CALL wrf_message('wrf_quilt_get_dom_ti_real not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_dom_ti_real 

SUBROUTINE wrf_quilt_put_dom_ti_real ( DataHandle,Element,   Data, Count,  Status )








  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  REAL ,          INTENT(IN)  :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status

  CHARACTER*132   :: locElement
  INTEGER i, typesize, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy


  CALL wrf_debug ( 50, 'in wrf_quilt_put_dom_ti_real' ) 
  CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
  locElement = Element

  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
      CALL MPI_TYPE_SIZE( MPI_REAL, typesize, ierr )

      IF ( wrf_dm_on_monitor() ) THEN
         CALL int_gen_ti_header( hdrbuf, hdrbufsize, itypesize, typesize, &
                                 DataHandle, locElement, Data, Count, int_dom_ti_real )
      ELSE
         CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
      ENDIF

      iserver = get_server_id ( DataHandle )
      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      CALL Mpi_Comm_Size (  comm_io_group , tasks_in_group , ierr  )


      
      reduced = 0
      reduced(1) = hdrbufsize 
      IF ( wrf_dm_on_monitor() )  reduced(2) = DataHandle
      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )

      
      CALL collect_on_comm_debug("module_io_quilt_old.F",3680, comm_io_group,            &
                            onebyte,                       &
                            hdrbuf, hdrbufsize , &
                            dummy, 0 )
    ENDIF
  ENDIF

  Status = 0

RETURN
END SUBROUTINE wrf_quilt_put_dom_ti_real 

SUBROUTINE wrf_quilt_get_dom_ti_double ( DataHandle,Element,   Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  real*8                      :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",7410,&
'wrf_quilt_get_dom_ti_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_dom_ti_double 

SUBROUTINE wrf_quilt_put_dom_ti_double ( DataHandle,Element,   Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  REAL*8 ,        INTENT(IN)  :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",7432,&
'wrf_quilt_put_dom_ti_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_put_dom_ti_double 

SUBROUTINE wrf_quilt_get_dom_ti_integer ( DataHandle,Element,   Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  integer                     :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status
  CALL wrf_message('wrf_quilt_get_dom_ti_integer not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_dom_ti_integer 

SUBROUTINE wrf_quilt_put_dom_ti_integer ( DataHandle,Element,   Data, Count,  Status )








  USE module_wrf_quilt
  USE module_state_description, ONLY: IO_PNETCDF
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  INTEGER ,       INTENT(IN)  :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status

  CHARACTER*132   :: locElement
  INTEGER i, typesize, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy
  INTEGER, EXTERNAL :: use_package



  locElement = Element

  CALL wrf_debug ( 50, 'in wrf_quilt_put_dom_ti_integer' ) 

  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )
      CALL MPI_TYPE_SIZE( MPI_INTEGER, typesize, ierr )



      IF ( wrf_dm_on_monitor() ) THEN
         CALL int_gen_ti_header( hdrbuf, hdrbufsize, itypesize, typesize, &
                                 DataHandle, locElement, Data, Count,     &
                                 int_dom_ti_integer )
      ELSE
         CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
      ENDIF

      iserver = get_server_id ( DataHandle )
      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      CALL Mpi_Comm_Size (  comm_io_group , tasks_in_group , ierr  )


      
      reduced = 0
      reduced(1) = hdrbufsize 
      IF ( wrf_dm_on_monitor() )  reduced(2) = DataHandle
      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )


      
      CALL collect_on_comm_debug("module_io_quilt_old.F",3840, comm_io_group,            &
                            onebyte,                       &
                            hdrbuf, hdrbufsize , &
                            dummy, 0 )
    ENDIF
  ENDIF
  CALL wrf_debug ( 50, 'returning from wrf_quilt_put_dom_ti_integer' ) 


RETURN
END SUBROUTINE wrf_quilt_put_dom_ti_integer 

SUBROUTINE wrf_quilt_get_dom_ti_logical ( DataHandle,Element,   Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  logical                     :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status

RETURN
END SUBROUTINE wrf_quilt_get_dom_ti_logical 

SUBROUTINE wrf_quilt_put_dom_ti_logical ( DataHandle,Element,   Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  logical ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status

  INTEGER i
  INTEGER one_or_zero(Count)

  DO i = 1, Count
    IF ( Data(i) ) THEN
      one_or_zero(i) = 1
    ELSE
      one_or_zero(i) = 0
    ENDIF
  ENDDO

  CALL wrf_quilt_put_dom_ti_integer ( DataHandle,Element,   one_or_zero, Count,  Status )
RETURN
END SUBROUTINE wrf_quilt_put_dom_ti_logical 

SUBROUTINE wrf_quilt_get_dom_ti_char ( DataHandle,Element,   Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*)               :: Data
  INTEGER                     :: Status
  CALL wrf_message('wrf_quilt_get_dom_ti_char not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_dom_ti_char 

SUBROUTINE wrf_quilt_put_dom_ti_char ( DataHandle, Element,  Data,  Status )








  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: Data
  INTEGER                     :: Status
  INTEGER i, itypesize, tasks_in_group, ierr, comm_io_group, me
  REAL dummy


  CALL wrf_debug ( 50, 'in wrf_quilt_put_dom_ti_char' ) 

  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )



      IF ( wrf_dm_on_monitor() ) THEN
         CALL int_gen_ti_header_char( hdrbuf, hdrbufsize, itypesize, &
                                      DataHandle, Element, "", Data, int_dom_ti_char )
      ELSE
         CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
      ENDIF

      iserver = get_server_id ( DataHandle )

      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      CALL Mpi_Comm_Size (  comm_io_group , tasks_in_group , ierr  )
      







      
      reduced_dummy = 0 
      reduced = 0
      reduced(1) = hdrbufsize 
      IF ( wrf_dm_on_monitor() )  reduced(2) = DataHandle


      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )


      


      CALL collect_on_comm_debug("module_io_quilt_old.F",4011, comm_io_group,            &
                            onebyte,                       &
                            hdrbuf, hdrbufsize , &
                            dummy, 0 )

    ENDIF
  ENDIF


RETURN
END SUBROUTINE wrf_quilt_put_dom_ti_char 

SUBROUTINE wrf_quilt_get_dom_td_real ( DataHandle,Element, DateStr,  Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  real                        :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_dom_td_real 

SUBROUTINE wrf_quilt_put_dom_td_real ( DataHandle,Element, DateStr,  Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  real ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_dom_td_real 

SUBROUTINE wrf_quilt_get_dom_td_double ( DataHandle,Element, DateStr,  Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  real*8                          :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",7810,&
'wrf_quilt_get_dom_td_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_dom_td_double 

SUBROUTINE wrf_quilt_put_dom_td_double ( DataHandle,Element, DateStr,  Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  real*8 ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",7833,&
'wrf_quilt_put_dom_td_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_put_dom_td_double 

SUBROUTINE wrf_quilt_get_dom_td_integer ( DataHandle,Element, DateStr,  Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  integer                          :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_dom_td_integer 

SUBROUTINE wrf_quilt_put_dom_td_integer ( DataHandle,Element, DateStr,  Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  integer ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_dom_td_integer 

SUBROUTINE wrf_quilt_get_dom_td_logical ( DataHandle,Element, DateStr,  Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  logical                          :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_dom_td_logical 

SUBROUTINE wrf_quilt_put_dom_td_logical ( DataHandle,Element, DateStr,  Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  logical ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_dom_td_logical 

SUBROUTINE wrf_quilt_get_dom_td_char ( DataHandle,Element, DateStr,  Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*)               :: Data
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_dom_td_char 

SUBROUTINE wrf_quilt_put_dom_td_char ( DataHandle,Element, DateStr,  Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN) :: Data
  INTEGER                          :: Status
RETURN
END SUBROUTINE wrf_quilt_put_dom_td_char 

SUBROUTINE wrf_quilt_get_var_ti_real ( DataHandle,Element,  Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real                          :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_ti_real 

SUBROUTINE wrf_quilt_put_var_ti_real ( DataHandle,Element,  Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_ti_real 

SUBROUTINE wrf_quilt_get_var_ti_double ( DataHandle,Element,  Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real*8                      :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",8030,&
'wrf_quilt_get_var_ti_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_var_ti_double 

SUBROUTINE wrf_quilt_put_var_ti_double ( DataHandle,Element,  Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real*8 ,        INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",8053,&
'wrf_quilt_put_var_ti_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_put_var_ti_double 

SUBROUTINE wrf_quilt_get_var_ti_integer ( DataHandle,Element,  Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  integer                     :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_ti_integer 

SUBROUTINE wrf_quilt_put_var_ti_integer ( DataHandle,Element,  Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  integer ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_ti_integer 

SUBROUTINE wrf_quilt_get_var_ti_logical ( DataHandle,Element,  Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  logical                     :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_ti_logical 

SUBROUTINE wrf_quilt_put_var_ti_logical ( DataHandle,Element,  Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  logical ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_ti_logical 

SUBROUTINE wrf_quilt_get_var_ti_char ( DataHandle,Element,  Varname, Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  CHARACTER*(*)               :: Data
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_ti_char 

SUBROUTINE wrf_quilt_put_var_ti_char ( DataHandle,Element,  Varname, Data,  Status )









  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER ::  int_ioexit			=  	     10
  INTEGER, PARAMETER ::  int_open_for_write_begin	=  	     20
  INTEGER, PARAMETER ::  int_open_for_write_commit	=  	     30
  INTEGER, PARAMETER ::  int_open_for_read 		=  	     40
  INTEGER, PARAMETER ::  int_inquire_opened 		=  	     60
  INTEGER, PARAMETER ::  int_inquire_filename 		=  	     70
  INTEGER, PARAMETER ::  int_iosync 			=  	     80
  INTEGER, PARAMETER ::  int_ioclose 			=  	     90
  INTEGER, PARAMETER ::  int_next_time 			=  	    100
  INTEGER, PARAMETER ::  int_set_time 			=  	    110
  INTEGER, PARAMETER ::  int_next_var 			=  	    120
  INTEGER, PARAMETER ::  int_dom_ti_real 		=  	    140
  INTEGER, PARAMETER ::  int_dom_ti_double 		=  	    160
  INTEGER, PARAMETER ::  int_dom_ti_integer 		=  	    180
  INTEGER, PARAMETER ::  int_dom_ti_logical 		=  	    200
  INTEGER, PARAMETER ::  int_dom_ti_char 		=  	    220
  INTEGER, PARAMETER ::  int_dom_td_real 		=  	    240
  INTEGER, PARAMETER ::  int_dom_td_double 		=  	    260
  INTEGER, PARAMETER ::  int_dom_td_integer 		=  	    280
  INTEGER, PARAMETER ::  int_dom_td_logical 		=  	    300
  INTEGER, PARAMETER ::  int_dom_td_char 		=  	    320
  INTEGER, PARAMETER ::  int_var_ti_real 		=  	    340
  INTEGER, PARAMETER ::  int_var_ti_double 		=  	    360
  INTEGER, PARAMETER ::  int_var_ti_integer 		=  	    380
  INTEGER, PARAMETER ::  int_var_ti_logical 		=  	    400
  INTEGER, PARAMETER ::  int_var_ti_char 		=  	    420
  INTEGER, PARAMETER ::  int_var_td_real 		=  	    440
  INTEGER, PARAMETER ::  int_var_td_double 		=  	    460
  INTEGER, PARAMETER ::  int_var_td_integer 		=  	    480
  INTEGER, PARAMETER ::  int_var_td_logical 		=  	    500
  INTEGER, PARAMETER ::  int_var_td_char 		=  	    520
  INTEGER, PARAMETER ::  int_field 			=  	    530
  INTEGER, PARAMETER ::  int_var_info 			=  	    540
  INTEGER, PARAMETER ::  int_noop 			=  	    550
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  CHARACTER*(*) , INTENT(IN)  :: Data
  INTEGER                     :: Status
  INTEGER i, itypesize, tasks_in_group, ierr, comm_io_group
  REAL dummy



  CALL wrf_debug ( 50, 'in wrf_quilt_put_var_ti_char' ) 

  IF ( DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles ) THEN
    IF ( int_handle_in_use( DataHandle ) ) THEN
      CALL MPI_TYPE_SIZE( MPI_INTEGER, itypesize, ierr )

      IF ( wrf_dm_on_monitor() ) THEN
         CALL int_gen_ti_header_char( hdrbuf, hdrbufsize, itypesize, &
                                      DataHandle, TRIM(Element),     &
                                      TRIM(VarName), TRIM(Data), int_var_ti_char )
      ELSE
         CALL int_gen_noop_header( hdrbuf, hdrbufsize, itypesize )
      ENDIF

      iserver = get_server_id ( DataHandle )
      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      CALL Mpi_Comm_Size (  comm_io_group , tasks_in_group , ierr  )


      
      reduced = 0
      reduced(1) = hdrbufsize 
      IF ( wrf_dm_on_monitor() )  reduced(2) = DataHandle
      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )

      
      CALL collect_on_comm_debug("module_io_quilt_old.F",4543, comm_io_group,            &
                            onebyte,                       &
                            hdrbuf, hdrbufsize , &
                            dummy, 0 )
    ENDIF
  ENDIF


RETURN
END SUBROUTINE wrf_quilt_put_var_ti_char 

SUBROUTINE wrf_quilt_get_var_td_real ( DataHandle,Element,  DateStr,Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real                        :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_td_real 

SUBROUTINE wrf_quilt_put_var_td_real ( DataHandle,Element,  DateStr,Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_td_real 

SUBROUTINE wrf_quilt_get_var_td_double ( DataHandle,Element,  DateStr,Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real*8                      :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",8327,&
'wrf_quilt_get_var_td_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_get_var_td_double 

SUBROUTINE wrf_quilt_put_var_td_double ( DataHandle,Element,  DateStr,Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  real*8 ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
  CALL wrf_error_fatal3("<stdin>",8351,&
'wrf_quilt_put_var_td_double not supported yet')
RETURN
END SUBROUTINE wrf_quilt_put_var_td_double 

SUBROUTINE wrf_quilt_get_var_td_integer ( DataHandle,Element,  DateStr,Varname, Data, Count, Outcount,Status)











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  integer                     :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_td_integer 

SUBROUTINE wrf_quilt_put_var_td_integer ( DataHandle,Element,  DateStr,Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  integer ,       INTENT(IN)  :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_td_integer 

SUBROUTINE wrf_quilt_get_var_td_logical ( DataHandle,Element,  DateStr,Varname, Data, Count, Outcount, Status )











  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  logical                          :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                      :: OutCount
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_td_logical 

SUBROUTINE wrf_quilt_put_var_td_logical ( DataHandle,Element,  DateStr,Varname, Data, Count,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  logical ,            INTENT(IN) :: Data(*)
  INTEGER ,       INTENT(IN)  :: Count
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_td_logical 

SUBROUTINE wrf_quilt_get_var_td_char ( DataHandle,Element,  DateStr,Varname, Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  CHARACTER*(*)               :: Data
  INTEGER                     :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_td_char 

SUBROUTINE wrf_quilt_put_var_td_char ( DataHandle,Element,  DateStr,Varname, Data,  Status )










  IMPLICIT NONE
  INTEGER ,       INTENT(IN)  :: DataHandle
  CHARACTER*(*) , INTENT(IN)  :: Element
  CHARACTER*(*) , INTENT(IN)  :: DateStr
  CHARACTER*(*) , INTENT(IN)  :: VarName 
  CHARACTER*(*) , INTENT(IN) :: Data
  INTEGER                    :: Status
RETURN
END SUBROUTINE wrf_quilt_put_var_td_char 

SUBROUTINE wrf_quilt_read_field ( DataHandle , DateStr , VarName , Field , FieldType , Comm , IOComm, &
                            DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                            DomainStart , DomainEnd ,                                    &
                            MemoryStart , MemoryEnd ,                                    &
                            PatchStart , PatchEnd ,                                      &
                            Status )







  IMPLICIT NONE
  INTEGER ,       INTENT(IN)    :: DataHandle 
  CHARACTER*(*) , INTENT(INOUT) :: DateStr
  CHARACTER*(*) , INTENT(INOUT) :: VarName
  INTEGER ,       INTENT(INOUT) :: Field(*)
  integer                       ,intent(in)    :: FieldType
  integer                       ,intent(inout) :: Comm
  integer                       ,intent(inout) :: IOComm
  integer                       ,intent(in)    :: DomainDesc
  character*(*)                 ,intent(in)    :: MemoryOrder
  character*(*)                 ,intent(in)    :: Stagger
  character*(*) , dimension (*) ,intent(in)    :: DimNames
  integer ,dimension(*)         ,intent(in)    :: DomainStart, DomainEnd
  integer ,dimension(*)         ,intent(in)    :: MemoryStart, MemoryEnd
  integer ,dimension(*)         ,intent(in)    :: PatchStart,  PatchEnd
  integer                       ,intent(out)   :: Status
  Status = 0
RETURN
END SUBROUTINE wrf_quilt_read_field

SUBROUTINE wrf_quilt_write_field ( DataHandle , DateStr , VarName , Field , FieldType , Comm , IOComm,  &
                             DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                             DomainStart , DomainEnd ,                                    &
                             MemoryStart , MemoryEnd ,                                    &
                             PatchStart , PatchEnd ,                                      &
                             Status )


















  USE module_state_description
  USE module_wrf_quilt
  IMPLICIT NONE
  INCLUDE 'mpif.h'
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101
      integer, parameter  :: WRF_FILE_OPENED_FOR_WRITE            = 102
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103
      integer, parameter  :: WRF_REAL                             = 105
      integer, parameter  :: WRF_DOUBLE                           = 105
      integer, parameter  :: WRF_FLOAT=WRF_REAL
      integer, parameter  :: WRF_INTEGER                          = 106
      integer, parameter  :: WRF_LOGICAL                          = 107
      integer, parameter  :: WRF_COMPLEX                          = 108
      integer, parameter  :: WRF_DOUBLE_COMPLEX                   = 109
      integer, parameter  :: WRF_FILE_OPENED_FOR_UPDATE           = 110


      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102
  INTEGER ,       INTENT(IN)    :: DataHandle 
  CHARACTER*(*) , INTENT(IN)    :: DateStr
  CHARACTER*(*) , INTENT(IN)    :: VarName

  integer                       ,intent(in)    :: FieldType
  integer                       ,intent(inout) :: Comm
  integer                       ,intent(inout) :: IOComm
  integer                       ,intent(in)    :: DomainDesc
  character*(*)                 ,intent(in)    :: MemoryOrder
  character*(*)                 ,intent(in)    :: Stagger
  character*(*) , dimension (*) ,intent(in)    :: DimNames
  integer ,dimension(*)         ,intent(in)    :: DomainStart, DomainEnd
  integer ,dimension(*)         ,intent(in)    :: MemoryStart, MemoryEnd
  integer ,dimension(*)         ,intent(in)    :: PatchStart,  PatchEnd
  integer                       ,intent(out)   :: Status

  integer ii,jj,kk,myrank

  REAL, DIMENSION( MemoryStart(1):MemoryEnd(1), &
                   MemoryStart(2):MemoryEnd(2), &
                   MemoryStart(3):MemoryEnd(3) ) :: Field
  INTEGER locsize , typesize, itypesize
  INTEGER ierr, tasks_in_group, comm_io_group, dummy, i
  INTEGER, EXTERNAL :: use_package


  CALL wrf_debug ( 50, 'in wrf_quilt_write_field' ) 

  IF ( .NOT. (DataHandle .GE. 1 .AND. DataHandle .LE. int_num_handles) ) THEN
    CALL wrf_error_fatal3("<stdin>",8595,&
"frame/module_io_quilt.F: wrf_quilt_write_field: invalid data handle" )
  ENDIF
  IF ( .NOT. int_handle_in_use( DataHandle ) ) THEN
    CALL wrf_error_fatal3("<stdin>",8599,&
"frame/module_io_quilt.F: wrf_quilt_write_field: DataHandle not opened" )
  ENDIF

  locsize = (PatchEnd(1)-PatchStart(1)+1)* &
            (PatchEnd(2)-PatchStart(2)+1)* &
            (PatchEnd(3)-PatchStart(3)+1)

  CALL mpi_type_size( MPI_INTEGER, itypesize, ierr )
  
  
  IF ( FieldType .EQ. WRF_DOUBLE ) THEN
    CALL mpi_type_size( MPI_DOUBLE_PRECISION, typesize, ierr )
  ELSE IF ( FieldType .EQ. WRF_FLOAT ) THEN
    CALL mpi_type_size( MPI_REAL, typesize, ierr )
  ELSE IF ( FieldType .EQ. WRF_INTEGER ) THEN
    CALL mpi_type_size( MPI_INTEGER, typesize, ierr )
  ELSE IF ( FieldType .EQ. WRF_LOGICAL ) THEN
    CALL mpi_type_size( MPI_LOGICAL, typesize, ierr )
  ENDIF

  IF ( .NOT. okay_to_write( DataHandle ) ) THEN

      
      
      
      

      CALL int_gen_write_field_header ( hdrbuf, hdrbufsize, itypesize, typesize,           &
                               DataHandle , DateStr , VarName , Field , FieldType , Comm , IOComm,  &
                               333933         , MemoryOrder , Stagger , DimNames ,              &   
                               DomainStart , DomainEnd ,                                    &
                               MemoryStart , MemoryEnd ,                                    &
                               PatchStart , PatchEnd )

      int_num_bytes_to_write(DataHandle) = int_num_bytes_to_write(DataHandle) + locsize * typesize + hdrbufsize

      

      iserver = get_server_id ( DataHandle )

      CALL get_mpi_comm_io_groups( comm_io_group , iserver )
      

      CALL Mpi_Comm_Size (  comm_io_group , tasks_in_group , ierr  )




      
      reduced = 0
      reduced(1) = hdrbufsize 
      IF ( wrf_dm_on_monitor() )  reduced(2) = DataHandle
      CALL MPI_Reduce( reduced, reduced_dummy, 2, MPI_INTEGER, MPI_SUM, tasks_in_group-1, comm_io_group, ierr )

      

      CALL collect_on_comm_debug("module_io_quilt_old.F",4965, comm_io_group,                   &
                            onebyte,                          &
                            hdrbuf, hdrbufsize ,                 &
                            dummy, 0 )

  ELSE

    IF ( .NOT. associated( int_local_output_buffer ) ) THEN
      ALLOCATE ( int_local_output_buffer( (int_num_bytes_to_write( DataHandle )+1)/itypesize ), Stat=ierr )
      IF(ierr /= 0)THEN
         CALL wrf_error_fatal3("<stdin>",8666,&
"frame/module_io_quilt.F: wrf_quilt_write_field: allocate of int_local_output_buffer failed" )
      END IF
      int_local_output_cursor = 1
    ENDIF
      iserver = get_server_id ( DataHandle )


    
    CALL int_gen_write_field_header ( hdrbuf, hdrbufsize, itypesize, typesize,           &
                             DataHandle , DateStr , VarName , Field , FieldType , Comm , IOComm,  &
                             0          , MemoryOrder , Stagger , DimNames ,              &   
                             DomainStart , DomainEnd ,                                    &
                             MemoryStart , MemoryEnd ,                                    &
                             PatchStart , PatchEnd )

    
    
    CALL int_pack_data ( hdrbuf , hdrbufsize , int_local_output_buffer, int_local_output_cursor )

    
    
    CALL int_pack_data ( Field(PatchStart(1):PatchEnd(1),PatchStart(2):PatchEnd(2),PatchStart(3):PatchEnd(3) ), &
                                  locsize * typesize , int_local_output_buffer, int_local_output_cursor )

  ENDIF
  Status = 0


  RETURN
END SUBROUTINE wrf_quilt_write_field

SUBROUTINE wrf_quilt_get_var_info ( DataHandle , VarName , NDim , MemoryOrder , Stagger , &
                              DomainStart , DomainEnd , Status )







  IMPLICIT NONE
  integer               ,intent(in)     :: DataHandle
  character*(*)         ,intent(in)     :: VarName
  integer                               :: NDim
  character*(*)                         :: MemoryOrder
  character*(*)                         :: Stagger
  integer ,dimension(*)                 :: DomainStart, DomainEnd
  integer                               :: Status
RETURN
END SUBROUTINE wrf_quilt_get_var_info

subroutine wrf_quilt_find_server(iserver)

  
  
  

  
  
  

  use module_wrf_quilt, only : in_avail, mpi_comm_avail, mpi_comm_local

  implicit none
  INCLUDE 'mpif.h'
  integer, intent(inout) :: iserver
  integer :: ierr
  character(255) :: message

  call wrf_message('Polling I/O servers...')

  if(in_avail) then
     call mpi_recv(iserver,1,MPI_INTEGER,MPI_ANY_SOURCE,0,mpi_comm_avail,MPI_STATUS_IGNORE,ierr)
     if(ierr/=0) then
        call wrf_error_fatal3("<stdin>",8741,&
'mpi_recv failed in wrf_quilt_find_server')
     endif
  endif

  call mpi_bcast(iserver,1,MPI_INTEGER,0,mpi_comm_local,ierr)
  if(ierr/=0) then
     call wrf_error_fatal3("<stdin>",8748,&
'mpi_bcast failed in wrf_quilt_find_server')
  endif

  write(message,'("I/O server ",I0," is ready for operations.")') iserver
  call wrf_message(message)


end subroutine wrf_quilt_find_server
subroutine wrf_quilt_server_ready()

  
  
  
  

  
  
  
  

  use module_wrf_quilt, only : mpi_comm_local, in_avail, availrank, mpi_comm_avail

  implicit none
  INCLUDE 'mpif.h'
  integer :: ierr
  character*255 :: message

  write(message,*) 'Entering wrf_quilt_server_ready.'
  call wrf_debug(1,message)

  call mpi_barrier(mpi_comm_local,ierr)
  if(ierr/=0) then
     call wrf_error_fatal3("<stdin>",8781,&
'mpi_barrier failed in wrf_quilt_server_ready')
  endif

  if(in_avail) then
     write(message,'("mpi_ssend ioserver=",I0," in wrf_quilt_server_ready")') availrank
     call wrf_debug(1,message)
     call mpi_ssend(availrank,1,MPI_INTEGER,0,0,mpi_comm_avail,ierr)
     if(ierr/=0) then
        call wrf_error_fatal3("<stdin>",8790,&
'mpi_ssend failed in wrf_quilt_server_ready')
     endif
  endif

  call mpi_barrier(mpi_comm_local,ierr)
  if(ierr/=0) then
     call wrf_error_fatal3("<stdin>",8797,&
'mpi_barrier failed in wrf_quilt_server_ready')
  endif

  write(message,*) 'Leaving wrf_quilt_server_ready.'
  call wrf_debug(1,message)

end subroutine wrf_quilt_server_ready

SUBROUTINE get_mpi_comm_io_groups( retval, isrvr )





      USE module_wrf_quilt
      IMPLICIT NONE
      INTEGER, INTENT(IN ) :: isrvr
      INTEGER, INTENT(OUT) :: retval
      retval = mpi_comm_io_groups(isrvr)
      RETURN
END SUBROUTINE get_mpi_comm_io_groups

SUBROUTINE get_nio_tasks_in_group( retval )





      USE module_wrf_quilt
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: retval
      retval = nio_tasks_in_group
      RETURN
END SUBROUTINE get_nio_tasks_in_group

SUBROUTINE collect_on_comm_debug(file,line, comm_io_group,   &
                        sze,                                 &
                        hdrbuf, hdrbufsize ,                 &
                        outbuf, outbufsize                   )
  IMPLICIT NONE
  CHARACTER*(*) file
  INTEGER line
  INTEGER comm_io_group
  INTEGER sze
  INTEGER hdrbuf(*), outbuf(*)
  INTEGER hdrbufsize, outbufsize 

  CALL collect_on_comm( comm_io_group,                       &
                        sze,                                 &
                        hdrbuf, hdrbufsize ,                 &
                        outbuf, outbufsize                   )
  RETURN
END


SUBROUTINE collect_on_comm_debug2(file,line,var,tag,sz,hdr_rec_size, &
                        comm_io_group,                       &
                        sze,                                 &
                        hdrbuf, hdrbufsize ,                 &
                        outbuf, outbufsize                   )
  IMPLICIT NONE
  CHARACTER*(*) file,var
  INTEGER line,tag,sz,hdr_rec_size
  INTEGER comm_io_group
  INTEGER sze
  INTEGER hdrbuf(*), outbuf(*)
  INTEGER hdrbufsize, outbufsize

  CALL collect_on_comm( comm_io_group,                       &
                        sze,                                 &
                        hdrbuf, hdrbufsize ,                 &
                        outbuf, outbufsize                   )
  RETURN
END

