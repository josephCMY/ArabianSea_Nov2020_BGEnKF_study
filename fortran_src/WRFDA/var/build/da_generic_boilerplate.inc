!
! WRFVAR generic type macro file
!
! This file is used to generate a series of simple boiler-plate calls 
! to support residual generic types for bitwise-exact testing.  
! It contains M4 macros and then
! a series of invocations of the macros to generate the subroutine
! definitions, which are then included in another source file.  
!

! $1 = specific ob name, $2 = specific ob type, $3 = ob counter



















!--- sound sound num_sound

SUBROUTINE da_y_type_ex_sound( iv, re, slice )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Extract all sound obs from y and place them in generic 
!           object slice.  
!           Call da_y_facade_free() to deallocate memory allocated here.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (iv_type),       INTENT(IN   ) :: iv     ! Innovation vector
   type (y_type),        INTENT(IN   ) :: re     ! all residual obs
   type (y_facade_type), INTENT(INOUT) :: slice  ! selected residual obs
   ! Local declarations
   INTEGER :: n

   CALL da_y_facade_create( slice, iv%info(sound)%nlocal, iv%info(sound)%ntotal )
   DO n=1, slice%num_obs
stop
!     CALL da_res_generic_set_info( slice%obs(n),                     &
!                                     iv%sound(n)%loc%proc_domain,      &
!                                     iv%sound(n)%loc%obs_global_index )
!     CALL da_res_sound_to_generic( re%sound(n), iv%sound(n)%info%levels, &
!                                     slice%obs(n) )
   ENDDO

END SUBROUTINE da_y_type_ex_sound  


!--- sound sound num_sound

SUBROUTINE da_y_type_ins_sound_global( slice_glob, re_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert obs from generic object slice_glob into 
!           globally-scoped y_type re_glob.  re_glob is 
!           allocated minimally.  Caller must deallocate.  
!           Memory for global object slice_glob is deallocated here.  
!           Do not use slice_glob after this call.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(INOUT) :: slice_glob ! generic
   type (y_type),        INTENT(INOUT) :: re_glob    ! selected residual obs
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize obs
   ! deallocation is done in free_global_sound()
   ALLOCATE( re_glob%sound(slice_glob%num_obs) )
   DO n=1, slice_glob%num_obs
     CALL da_res_sound_from_generic( slice_glob%obs(n), re_glob%sound(n) )
   ENDDO
   re_glob%nlocal(sound) = slice_glob%num_obs  ! duplication!
   CALL da_y_facade_free( slice_glob )

END SUBROUTINE da_y_type_ins_sound_global 


!--- sound sound num_sound

SUBROUTINE da_iv_type_ins_sound_global( slice_glob, iv_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert meta-data from generic object slice_glob into 
!           globally-scoped iv_type iv_glob.  iv_glob is 
!           allocated minimally.  Caller must deallocate.  
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(IN   ) :: slice_glob ! selected residual obs
   type (iv_type),       INTENT(INOUT) :: iv_glob    ! partial Innovation vector
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize needed bits of iv_glob (ugly)
   iv_glob%info(sound)%nlocal  = slice_glob%num_obs
   iv_glob%info(sound)%ntotal = slice_glob%num_obs_glo
   ! deallocation is done in free_global_sound()
   ALLOCATE( iv_glob%sound(iv_glob%info(sound)%nlocal) )
   DO n=1, iv_glob%info(sound)%nlocal
stop
!     iv_glob%sound(n)%loc%proc_domain = slice_glob%obs(n)%proc_domain
!     iv_glob%sound(n)%loc%obs_global_index = &
!                                        slice_glob%obs(n)%obs_global_index
!     IF ( da_res_generic_has_vector( slice_glob%obs(n) ) ) THEN
!       iv_glob%sound(n)%info%levels = SIZE(slice_glob%obs(n)%values(1)%ptr)
!     ENDIF
   ENDDO

END SUBROUTINE da_iv_type_ins_sound_global  


!--- sound sound num_sound

!------------------------------------------------------------------------------
! PURPOSE:  Collect local arrays of residual_sound_type objects into 
!           global arrays in serial-code storage order.  This is used to 
!           perform bitwise-exact floating-point summations in 
!           serial-code-order during bitwise-exact testing of 
!           distributed-memory parallel configurations.  
!
! METHOD:   Indices stowed away during Read_Obs() are used to restore serial 
!           storage order.  Memory for global objects is allocated here.  
!           Global objects are minimally allocated to save memory.  
!           Memory is deallocated in free_global_sound().  
!------------------------------------------------------------------------------
  SUBROUTINE da_to_global_sound( iv,      re,      jo_grad_y, &
                                  iv_glob, re_glob, jo_grad_y_glob )

    IMPLICIT NONE

    ! task-local objects
    type (iv_type), INTENT( IN) :: iv             ! Innovation vector
    type (y_type),  INTENT( IN) :: re             ! residual vector
    type (y_type),  INTENT( IN) :: jo_grad_y      ! Grad_y(Jo)
    ! task-global objects
    type (iv_type), INTENT(OUT) :: iv_glob        ! Innovation vector
    type (y_type),  INTENT(OUT) :: re_glob        ! residual vector
    type (y_type),  INTENT(OUT) :: jo_grad_y_glob ! Grad_y(Jo)

    ! Local declarations
    type (y_facade_type) :: re_slice, re_glob_slice
    type (y_facade_type) :: jo_grad_y_slice, jo_grad_y_glob_slice
    type (residual_template_type) :: template  ! allocation info

    ! create process-local generic objects from specific objects
    CALL da_y_type_ex_sound( iv, re,        re_slice )
    CALL da_y_type_ex_sound( iv, jo_grad_y, jo_grad_y_slice )

    ! create global versions of generic objects from process-local objects
    ! and destroy process-local generic objects
    CALL da_res_sound_create_template( template )  ! use template in case 
                                                     ! some tasks have no obs
    CALL da_y_facade_to_global( re_slice,        template, re_glob_slice )
    CALL da_y_facade_to_global( jo_grad_y_slice, template, jo_grad_y_glob_slice )

    ! create global versions of specific objects
    ! and destroy global versions of generic objects
    ! iv first
    CALL da_iv_type_ins_sound_global( re_glob_slice, iv_glob )
    ! then y_types
    CALL da_y_type_ins_sound_global( re_glob_slice,        re_glob )
    CALL da_y_type_ins_sound_global( jo_grad_y_glob_slice, jo_grad_y_glob )
    ! global versions of specific objects are destroyed in 
    ! free_global_sound()

    RETURN

  END SUBROUTINE da_to_global_sound  


!--- sonde_sfc synop num_sound

SUBROUTINE da_y_type_ex_sonde_sfc( iv, re, slice )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Extract all sonde_sfc obs from y and place them in generic 
!           object slice.  
!           Call da_y_facade_free() to deallocate memory allocated here.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (iv_type),       INTENT(IN   ) :: iv     ! Innovation vector
   type (y_type),        INTENT(IN   ) :: re     ! all residual obs
   type (y_facade_type), INTENT(INOUT) :: slice  ! selected residual obs
   ! Local declarations
   INTEGER :: n

   CALL da_y_facade_create( slice, iv%info(sonde_sfc)%nlocal, iv%info(sonde_sfc)%ntotal )
   DO n=1, slice%num_obs
stop
!     CALL da_res_generic_set_info( slice%obs(n),                     &
!                                     iv%sonde_sfc(n)%loc%proc_domain,      &
!                                     iv%sonde_sfc(n)%loc%obs_global_index )
!     CALL da_res_synop_to_generic( re%sonde_sfc(n), iv%sonde_sfc(n)%info%levels, &
!                                     slice%obs(n) )
   ENDDO

END SUBROUTINE da_y_type_ex_sonde_sfc  


!--- sonde_sfc synop num_sound

SUBROUTINE da_y_type_ins_sonde_sfc_global( slice_glob, re_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert obs from generic object slice_glob into 
!           globally-scoped y_type re_glob.  re_glob is 
!           allocated minimally.  Caller must deallocate.  
!           Memory for global object slice_glob is deallocated here.  
!           Do not use slice_glob after this call.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(INOUT) :: slice_glob ! generic
   type (y_type),        INTENT(INOUT) :: re_glob    ! selected residual obs
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize obs
   ! deallocation is done in free_global_sonde_sfc()
   ALLOCATE( re_glob%sonde_sfc(slice_glob%num_obs) )
   DO n=1, slice_glob%num_obs
     CALL da_res_synop_from_generic( slice_glob%obs(n), re_glob%sonde_sfc(n) )
   ENDDO
   re_glob%nlocal(sonde_sfc) = slice_glob%num_obs  ! duplication!
   CALL da_y_facade_free( slice_glob )

END SUBROUTINE da_y_type_ins_sonde_sfc_global 


!--- sonde_sfc synop num_sound

SUBROUTINE da_iv_type_ins_sonde_sfc_global( slice_glob, iv_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert meta-data from generic object slice_glob into 
!           globally-scoped iv_type iv_glob.  iv_glob is 
!           allocated minimally.  Caller must deallocate.  
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(IN   ) :: slice_glob ! selected residual obs
   type (iv_type),       INTENT(INOUT) :: iv_glob    ! partial Innovation vector
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize needed bits of iv_glob (ugly)
   iv_glob%info(sonde_sfc)%nlocal  = slice_glob%num_obs
   iv_glob%info(sonde_sfc)%ntotal = slice_glob%num_obs_glo
   ! deallocation is done in free_global_sonde_sfc()
   ALLOCATE( iv_glob%sonde_sfc(iv_glob%info(sonde_sfc)%nlocal) )
   DO n=1, iv_glob%info(sonde_sfc)%nlocal
stop
!     iv_glob%sonde_sfc(n)%loc%proc_domain = slice_glob%obs(n)%proc_domain
!     iv_glob%sonde_sfc(n)%loc%obs_global_index = &
!                                        slice_glob%obs(n)%obs_global_index
!     IF ( da_res_generic_has_vector( slice_glob%obs(n) ) ) THEN
!       iv_glob%sonde_sfc(n)%info%levels = SIZE(slice_glob%obs(n)%values(1)%ptr)
!     ENDIF
   ENDDO

END SUBROUTINE da_iv_type_ins_sonde_sfc_global  


!--- sonde_sfc synop num_sound

!------------------------------------------------------------------------------
! PURPOSE:  Collect local arrays of residual_synop_type objects into 
!           global arrays in serial-code storage order.  This is used to 
!           perform bitwise-exact floating-point summations in 
!           serial-code-order during bitwise-exact testing of 
!           distributed-memory parallel configurations.  
!
! METHOD:   Indices stowed away during Read_Obs() are used to restore serial 
!           storage order.  Memory for global objects is allocated here.  
!           Global objects are minimally allocated to save memory.  
!           Memory is deallocated in free_global_sonde_sfc().  
!------------------------------------------------------------------------------
  SUBROUTINE da_to_global_sonde_sfc( iv,      re,      jo_grad_y, &
                                  iv_glob, re_glob, jo_grad_y_glob )

    IMPLICIT NONE

    ! task-local objects
    type (iv_type), INTENT( IN) :: iv             ! Innovation vector
    type (y_type),  INTENT( IN) :: re             ! residual vector
    type (y_type),  INTENT( IN) :: jo_grad_y      ! Grad_y(Jo)
    ! task-global objects
    type (iv_type), INTENT(OUT) :: iv_glob        ! Innovation vector
    type (y_type),  INTENT(OUT) :: re_glob        ! residual vector
    type (y_type),  INTENT(OUT) :: jo_grad_y_glob ! Grad_y(Jo)

    ! Local declarations
    type (y_facade_type) :: re_slice, re_glob_slice
    type (y_facade_type) :: jo_grad_y_slice, jo_grad_y_glob_slice
    type (residual_template_type) :: template  ! allocation info

    ! create process-local generic objects from specific objects
    CALL da_y_type_ex_sonde_sfc( iv, re,        re_slice )
    CALL da_y_type_ex_sonde_sfc( iv, jo_grad_y, jo_grad_y_slice )

    ! create global versions of generic objects from process-local objects
    ! and destroy process-local generic objects
    CALL da_res_synop_create_template( template )  ! use template in case 
                                                     ! some tasks have no obs
    CALL da_y_facade_to_global( re_slice,        template, re_glob_slice )
    CALL da_y_facade_to_global( jo_grad_y_slice, template, jo_grad_y_glob_slice )

    ! create global versions of specific objects
    ! and destroy global versions of generic objects
    ! iv first
    CALL da_iv_type_ins_sonde_sfc_global( re_glob_slice, iv_glob )
    ! then y_types
    CALL da_y_type_ins_sonde_sfc_global( re_glob_slice,        re_glob )
    CALL da_y_type_ins_sonde_sfc_global( jo_grad_y_glob_slice, jo_grad_y_glob )
    ! global versions of specific objects are destroyed in 
    ! free_global_sonde_sfc()

    RETURN

  END SUBROUTINE da_to_global_sonde_sfc  


!--- synop synop num_synop

SUBROUTINE da_y_type_ex_synop( iv, re, slice )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Extract all synop obs from y and place them in generic 
!           object slice.  
!           Call da_y_facade_free() to deallocate memory allocated here.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (iv_type),       INTENT(IN   ) :: iv     ! Innovation vector
   type (y_type),        INTENT(IN   ) :: re     ! all residual obs
   type (y_facade_type), INTENT(INOUT) :: slice  ! selected residual obs
   ! Local declarations
   INTEGER :: n

   CALL da_y_facade_create( slice, iv%info(synop)%nlocal, iv%info(synop)%ntotal )
   DO n=1, slice%num_obs
stop
!     CALL da_res_generic_set_info( slice%obs(n),                     &
!                                     iv%synop(n)%loc%proc_domain,      &
!                                     iv%synop(n)%loc%obs_global_index )
!     CALL da_res_synop_to_generic( re%synop(n), iv%synop(n)%info%levels, &
!                                     slice%obs(n) )
   ENDDO

END SUBROUTINE da_y_type_ex_synop  


!--- synop synop num_synop

SUBROUTINE da_y_type_ins_synop_global( slice_glob, re_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert obs from generic object slice_glob into 
!           globally-scoped y_type re_glob.  re_glob is 
!           allocated minimally.  Caller must deallocate.  
!           Memory for global object slice_glob is deallocated here.  
!           Do not use slice_glob after this call.
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(INOUT) :: slice_glob ! generic
   type (y_type),        INTENT(INOUT) :: re_glob    ! selected residual obs
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize obs
   ! deallocation is done in free_global_synop()
   ALLOCATE( re_glob%synop(slice_glob%num_obs) )
   DO n=1, slice_glob%num_obs
     CALL da_res_synop_from_generic( slice_glob%obs(n), re_glob%synop(n) )
   ENDDO
   re_glob%nlocal(synop) = slice_glob%num_obs  ! duplication!
   CALL da_y_facade_free( slice_glob )

END SUBROUTINE da_y_type_ins_synop_global 


!--- synop synop num_synop

SUBROUTINE da_iv_type_ins_synop_global( slice_glob, iv_glob )

!------------------------------------------------------------------------------
! PURPOSE:  Eliminate redundant code for many obs types.  
!
! METHOD:   Insert meta-data from generic object slice_glob into 
!           globally-scoped iv_type iv_glob.  iv_glob is 
!           allocated minimally.  Caller must deallocate.  
!------------------------------------------------------------------------------
   IMPLICIT NONE

   type (y_facade_type), INTENT(IN   ) :: slice_glob ! selected residual obs
   type (iv_type),       INTENT(INOUT) :: iv_glob    ! partial Innovation vector
   ! Local declarations
   INTEGER :: n

   ! allocate and initialize needed bits of iv_glob (ugly)
   iv_glob%info(synop)%nlocal  = slice_glob%num_obs
   iv_glob%info(synop)%ntotal = slice_glob%num_obs_glo
   ! deallocation is done in free_global_synop()
   ALLOCATE( iv_glob%synop(iv_glob%info(synop)%nlocal) )
   DO n=1, iv_glob%info(synop)%nlocal
stop
!     iv_glob%synop(n)%loc%proc_domain = slice_glob%obs(n)%proc_domain
!     iv_glob%synop(n)%loc%obs_global_index = &
!                                        slice_glob%obs(n)%obs_global_index
!     IF ( da_res_generic_has_vector( slice_glob%obs(n) ) ) THEN
!       iv_glob%synop(n)%info%levels = SIZE(slice_glob%obs(n)%values(1)%ptr)
!     ENDIF
   ENDDO

END SUBROUTINE da_iv_type_ins_synop_global  


!--- synop synop num_synop

!------------------------------------------------------------------------------
! PURPOSE:  Collect local arrays of residual_synop_type objects into 
!           global arrays in serial-code storage order.  This is used to 
!           perform bitwise-exact floating-point summations in 
!           serial-code-order during bitwise-exact testing of 
!           distributed-memory parallel configurations.  
!
! METHOD:   Indices stowed away during Read_Obs() are used to restore serial 
!           storage order.  Memory for global objects is allocated here.  
!           Global objects are minimally allocated to save memory.  
!           Memory is deallocated in free_global_synop().  
!------------------------------------------------------------------------------
  SUBROUTINE da_to_global_synop( iv,      re,      jo_grad_y, &
                                  iv_glob, re_glob, jo_grad_y_glob )

    IMPLICIT NONE

    ! task-local objects
    type (iv_type), INTENT( IN) :: iv             ! Innovation vector
    type (y_type),  INTENT( IN) :: re             ! residual vector
    type (y_type),  INTENT( IN) :: jo_grad_y      ! Grad_y(Jo)
    ! task-global objects
    type (iv_type), INTENT(OUT) :: iv_glob        ! Innovation vector
    type (y_type),  INTENT(OUT) :: re_glob        ! residual vector
    type (y_type),  INTENT(OUT) :: jo_grad_y_glob ! Grad_y(Jo)

    ! Local declarations
    type (y_facade_type) :: re_slice, re_glob_slice
    type (y_facade_type) :: jo_grad_y_slice, jo_grad_y_glob_slice
    type (residual_template_type) :: template  ! allocation info

    ! create process-local generic objects from specific objects
    CALL da_y_type_ex_synop( iv, re,        re_slice )
    CALL da_y_type_ex_synop( iv, jo_grad_y, jo_grad_y_slice )

    ! create global versions of generic objects from process-local objects
    ! and destroy process-local generic objects
    CALL da_res_synop_create_template( template )  ! use template in case 
                                                     ! some tasks have no obs
    CALL da_y_facade_to_global( re_slice,        template, re_glob_slice )
    CALL da_y_facade_to_global( jo_grad_y_slice, template, jo_grad_y_glob_slice )

    ! create global versions of specific objects
    ! and destroy global versions of generic objects
    ! iv first
    CALL da_iv_type_ins_synop_global( re_glob_slice, iv_glob )
    ! then y_types
    CALL da_y_type_ins_synop_global( re_glob_slice,        re_glob )
    CALL da_y_type_ins_synop_global( jo_grad_y_glob_slice, jo_grad_y_glob )
    ! global versions of specific objects are destroyed in 
    ! free_global_synop()

    RETURN

  END SUBROUTINE da_to_global_synop  

