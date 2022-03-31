

























MODULE module_internal_header_util
































INTERFACE int_get_ti_header
   MODULE PROCEDURE int_get_ti_header_integer, int_get_ti_header_real
END INTERFACE
INTERFACE int_gen_ti_header
   MODULE PROCEDURE int_gen_ti_header_integer, int_gen_ti_header_real
END INTERFACE
INTERFACE int_get_td_header
   MODULE PROCEDURE int_get_td_header_integer, int_get_td_header_real
END INTERFACE
INTERFACE int_gen_td_header
   MODULE PROCEDURE int_gen_td_header_integer, int_gen_td_header_real
END INTERFACE

PRIVATE :: int_pack_string, int_unpack_string

CONTAINS


INTEGER FUNCTION get_hdr_tag( hdrbuf )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: hdrbuf(*)
  get_hdr_tag = hdrbuf(2)
  RETURN
END FUNCTION get_hdr_tag

INTEGER FUNCTION get_hdr_rec_size( hdrbuf )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: hdrbuf(*)
  get_hdr_rec_size = hdrbuf(1)
  RETURN
END FUNCTION get_hdr_rec_size

SUBROUTINE int_gen_write_field_header ( hdrbuf, hdrbufsize, itypesize, ftypesize, &
                                        DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm, &
                                        DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                        DomainStart , DomainEnd ,                                    &
                                        MemoryStart , MemoryEnd ,                                    &
                                        PatchStart , PatchEnd )




















































  IMPLICIT NONE
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
  INTEGER,       INTENT(INOUT)  ::  hdrbuf(*)
  INTEGER,       INTENT(INOUT)  ::  hdrbufsize
  INTEGER,       INTENT(INOUT)  ::  itypesize, ftypesize
  INTEGER ,      INTENT(IN)     :: DataHandle
  CHARACTER*(*), INTENT(IN)  :: DateStr
  CHARACTER*(*), INTENT(IN)  :: VarName
  REAL, DIMENSION(*)            :: Dummy
  INTEGER                       ,intent(in)    :: FieldType
  INTEGER                       ,intent(inout) :: Comm
  INTEGER                       ,intent(inout) :: IOComm
  INTEGER                       ,intent(in)    :: DomainDesc
  CHARACTER*(*)                 ,intent(in)    :: MemoryOrder
  CHARACTER*(*)                 ,intent(in)    :: Stagger
  CHARACTER*(*) , dimension (*) ,intent(in)    :: DimNames
  INTEGER ,dimension(*)         ,intent(in)    :: DomainStart, DomainEnd
  INTEGER ,dimension(*)         ,intent(in)    :: MemoryStart, MemoryEnd
  INTEGER ,dimension(*)         ,intent(in)    :: PatchStart,  PatchEnd

  INTEGER i, n


  hdrbuf(1) = 0 
  hdrbuf(2) = int_field
  hdrbuf(3) = ftypesize

  i = 4
  hdrbuf(i) = DataHandle      ; i = i+1
  call int_pack_string( DateStr, hdrbuf(i), n ) ; i = i + n
  call int_pack_string( VarName, hdrbuf(i), n ) ; i = i + n
  hdrbuf(i) = FieldType       ; i = i+1
  call int_pack_string( MemoryOrder, hdrbuf(i), n ) ; i = i + n
  call int_pack_string( Stagger,     hdrbuf(i), n ) ; i = i + n
  call int_pack_string( DimNames(1), hdrbuf(i), n ) ; i = i + n
  call int_pack_string( DimNames(2), hdrbuf(i), n ) ; i = i + n
  call int_pack_string( DimNames(3), hdrbuf(i), n ) ; i = i + n
  hdrbuf(i) = DomainStart(1)     ; i = i+1
  hdrbuf(i) = DomainStart(2)     ; i = i+1
  hdrbuf(i) = DomainStart(3)     ; i = i+1
  hdrbuf(i) = DomainEnd(1)       ; i = i+1
  hdrbuf(i) = DomainEnd(2)       ; i = i+1
  hdrbuf(i) = DomainEnd(3)       ; i = i+1
  hdrbuf(i) = PatchStart(1)     ; i = i+1
  hdrbuf(i) = PatchStart(2)     ; i = i+1
  hdrbuf(i) = PatchStart(3)     ; i = i+1
  hdrbuf(i) = PatchEnd(1)       ; i = i+1
  hdrbuf(i) = PatchEnd(2)       ; i = i+1
  hdrbuf(i) = PatchEnd(3)       ; i = i+1
  hdrbuf(i) = DomainDesc        ; i = i+1

  hdrbufsize = (i-1) * itypesize  
  hdrbuf(1) = hdrbufsize

  RETURN
END SUBROUTINE int_gen_write_field_header


SUBROUTINE int_get_write_field_header ( hdrbuf, hdrbufsize, itypesize, ftypesize, &
                                        DataHandle , DateStr , VarName , Dummy , FieldType , Comm , IOComm,  &
                                        DomainDesc , MemoryOrder , Stagger , DimNames ,              &
                                        DomainStart , DomainEnd ,                                    &
                                        MemoryStart , MemoryEnd ,                                    &
                                        PatchStart , PatchEnd )






  IMPLICIT NONE
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
  INTEGER,       INTENT(INOUT)  ::  hdrbuf(*)
  INTEGER,       INTENT(OUT)    ::  hdrbufsize
  INTEGER,       INTENT(INOUT)  ::  itypesize, ftypesize
  INTEGER ,      INTENT(OUT)    :: DataHandle
  CHARACTER*(*), INTENT(INOUT)  :: DateStr
  CHARACTER*(*), INTENT(INOUT)  :: VarName
  REAL, DIMENSION(*)            :: Dummy
  INTEGER                                       :: FieldType
  INTEGER                                       :: Comm
  INTEGER                                       :: IOComm
  INTEGER                                       :: DomainDesc
  CHARACTER*(*)                                 :: MemoryOrder
  CHARACTER*(*)                                 :: Stagger
  CHARACTER*(*) , dimension (*)                 :: DimNames
  INTEGER ,dimension(*)                         :: DomainStart, DomainEnd
  INTEGER ,dimension(*)                         :: MemoryStart, MemoryEnd
  INTEGER ,dimension(*)                         :: PatchStart,  PatchEnd

  CHARACTER*132 mess
  INTEGER i, n

  hdrbufsize = hdrbuf(1)
  IF ( hdrbuf(2) .NE. int_field ) THEN
    write(mess,*)'int_get_write_field_header: hdrbuf(2) ne int_field ',hdrbuf(2),int_field
    CALL wrf_error_fatal3("<stdin>",299,&
mess )
  ENDIF
  ftypesize = hdrbuf(3)

   i = 4
   DataHandle = hdrbuf(i)     ; i = i+1
  call int_unpack_string( DateStr, hdrbuf(i), n )     ; i = i+n
  call int_unpack_string( VarName, hdrbuf(i), n )     ; i = i+n
   FieldType = hdrbuf(i)      ; i = i+1
  call int_unpack_string( MemoryOrder, hdrbuf(i), n ) ; i = i+n
  call int_unpack_string( Stagger, hdrbuf(i), n )     ; i = i+n
  call int_unpack_string( DimNames(1), hdrbuf(i), n ) ; i = i+n
  call int_unpack_string( DimNames(2), hdrbuf(i), n ) ; i = i+n
  call int_unpack_string( DimNames(3), hdrbuf(i), n ) ; i = i+n
   DomainStart(1) = hdrbuf(i)    ; i = i+1
   DomainStart(2) = hdrbuf(i)    ; i = i+1
   DomainStart(3) = hdrbuf(i)    ; i = i+1
   DomainEnd(1) = hdrbuf(i)       ; i = i+1
   DomainEnd(2) = hdrbuf(i)       ; i = i+1
   DomainEnd(3) = hdrbuf(i)       ; i = i+1
   PatchStart(1) = hdrbuf(i)     ; i = i+1
   PatchStart(2) = hdrbuf(i)     ; i = i+1
   PatchStart(3) = hdrbuf(i)     ; i = i+1
   PatchEnd(1) = hdrbuf(i)       ; i = i+1
   PatchEnd(2) = hdrbuf(i)       ; i = i+1
   PatchEnd(3) = hdrbuf(i)       ; i = i+1
   DomainDesc = hdrbuf(i)       ; i = i+1

  RETURN
END SUBROUTINE int_get_write_field_header




SUBROUTINE int_gen_ofr_header( hdrbuf, hdrbufsize, itypesize, &
                                FileName, SysDepInfo, DataHandle )



























  IMPLICIT NONE
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
  INTEGER,       INTENT(INOUT) ::  hdrbuf(*)
  INTEGER,       INTENT(OUT)   ::  hdrbufsize
  INTEGER,       INTENT(INOUT) ::  itypesize
  INTEGER ,      INTENT(IN)    :: DataHandle
  CHARACTER*(*), INTENT(INOUT) :: FileName
  CHARACTER*(*), INTENT(INOUT) :: SysDepInfo

  INTEGER i, n, i1

  hdrbuf(1) = 0  
  hdrbuf(2) = int_open_for_read
  i = 3
  hdrbuf(i) = DataHandle     ; i = i+1

  call int_pack_string( TRIM(FileName), hdrbuf(i), n )   ; i = i + n
  call int_pack_string( TRIM(SysDepInfo), hdrbuf(i), n ) ; i = i + n
  hdrbufsize = (i-1) * itypesize  
  hdrbuf(1) = hdrbufsize
  RETURN
END SUBROUTINE int_gen_ofr_header


SUBROUTINE int_get_ofr_header( hdrbuf, hdrbufsize, itypesize, &
                                FileName, SysDepInfo, DataHandle )






  IMPLICIT NONE
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
  INTEGER,       INTENT(INOUT) ::  hdrbuf(*)
  INTEGER,       INTENT(OUT)   ::  hdrbufsize
  INTEGER,       INTENT(INOUT) ::  itypesize
  INTEGER ,      INTENT(OUT)   :: DataHandle
  CHARACTER*(*), INTENT(INOUT) :: FileName
  CHARACTER*(*), INTENT(INOUT) :: SysDepInfo

  INTEGER i, n

  hdrbufsize = hdrbuf(1)



  i = 3
  DataHandle = hdrbuf(i)    ; i = i+1
  call int_unpack_string( FileName, hdrbuf(i), n ) ; i = i+n
  call int_unpack_string( SysDepInfo, hdrbuf(i), n ) ; i = i+n
  RETURN
END SUBROUTINE int_get_ofr_header




SUBROUTINE int_gen_ofwb_header( hdrbuf, hdrbufsize, itypesize, &
                                FileName, SysDepInfo, io_form, DataHandle )





























  IMPLICIT NONE
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
  INTEGER,       INTENT(INOUT) :: hdrbuf(*)
  INTEGER,       INTENT(OUT)   :: hdrbufsize
  INTEGER,       INTENT(INOUT) :: itypesize
  INTEGER ,      INTENT(IN)    :: io_form
  INTEGER ,      INTENT(IN)    :: DataHandle
  CHARACTER*(*), INTENT(INOUT) :: FileName
  CHARACTER*(*), INTENT(INOUT) :: SysDepInfo

  INTEGER i, n, j

  hdrbuf(1) = 0  
  hdrbuf(2) = int_open_for_write_begin
  i = 3
  hdrbuf(i) = DataHandle     ; i = i+1
  hdrbuf(i) = io_form        ; i = i+1

  call int_pack_string( FileName, hdrbuf(i), n )   ; i = i + n


  call int_pack_string( SysDepInfo, hdrbuf(i), n ) ; i = i + n

  hdrbufsize = (i-1) * itypesize  
  hdrbuf(1) = hdrbufsize

  RETURN
END SUBROUTINE int_gen_ofwb_header


SUBROUTINE int_get_ofwb_header( hdrbuf, hdrbufsize, itypesize, &
                                FileName, SysDepInfo, io_form, DataHandle )






  IMPLICIT NONE
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
  INTEGER,       INTENT(INOUT)  :: hdrbuf(*)
  INTEGER,       INTENT(OUT)    :: hdrbufsize
  INTEGER,       INTENT(INOUT)  :: itypesize
  INTEGER ,      INTENT(OUT)    :: DataHandle
  INTEGER ,      INTENT(OUT)    :: io_form
  CHARACTER*(*), INTENT (INOUT) :: FileName
  CHARACTER*(*), INTENT (INOUT) :: SysDepInfo

  INTEGER i, n, j

  hdrbufsize = hdrbuf(1)




  i = 3
  DataHandle = hdrbuf(i)    ; i = i+1

  io_form    = hdrbuf(i)    ; i = i+1



  call int_unpack_string( FileName, hdrbuf(i), n ) ; i = i+n



  call int_unpack_string( SysDepInfo, hdrbuf(i), n ) ; i = i+n



  RETURN
END SUBROUTINE int_get_ofwb_header



SUBROUTINE int_gen_handle_header( hdrbuf, hdrbufsize, itypesize, &
                                DataHandle , code )




















  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT) ::  hdrbuf(*)
  INTEGER, INTENT(OUT)   ::  hdrbufsize
  INTEGER, INTENT(INOUT) ::  itypesize
  INTEGER ,INTENT(IN)    :: DataHandle, code

  INTEGER i

  hdrbuf(1) = 0  
  hdrbuf(2) = code
  i = 3
  hdrbuf(i) = DataHandle     ; i = i+1
  hdrbufsize = (i-1) * itypesize  
  hdrbuf(1) = hdrbufsize
  RETURN
END SUBROUTINE int_gen_handle_header

SUBROUTINE int_get_handle_header( hdrbuf, hdrbufsize, itypesize, &
                                DataHandle , code )






  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT) ::  hdrbuf(*)
  INTEGER, INTENT(OUT)   ::  hdrbufsize
  INTEGER, INTENT(INOUT) ::  itypesize
  INTEGER ,INTENT(OUT)   :: DataHandle, code

  INTEGER i

  hdrbufsize = hdrbuf(1)
  code       = hdrbuf(2)
  i = 3
  DataHandle = hdrbuf(i)    ; i = i+1
  RETURN
END SUBROUTINE int_get_handle_header



SUBROUTINE int_gen_ti_header_integer( hdrbuf, hdrbufsize, itypesize, typesize, &
                                      DataHandle, Element, Data, Count, code )






























  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize, typesize
  CHARACTER*(*), INTENT(INOUT) ::  Element
  INTEGER, INTENT(IN)          ::  Data(*)
  INTEGER, INTENT(IN)          ::  DataHandle, Count, code

  INTEGER i, n

  CALL int_gen_ti_header_c ( hdrbuf, hdrbufsize, itypesize, typesize, &
                             DataHandle, Data, Count, code )
  i = hdrbufsize/itypesize + 1 ;

  CALL int_pack_string ( Element, hdrbuf( i ), n ) ; i = i + n
  hdrbufsize = n * itypesize + hdrbufsize 
  hdrbuf(1) = hdrbufsize
  RETURN
END SUBROUTINE int_gen_ti_header_integer

SUBROUTINE int_gen_ti_header_integer_varna( hdrbuf, hdrbufsize, itypesize, typesize, &
                                      DataHandle, Element, VarName, Data, Count, code )


































  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize, typesize
  CHARACTER*(*), INTENT(IN)    ::  Element, VarName
  INTEGER, INTENT(IN)          ::  Data(*)
  INTEGER, INTENT(IN)          ::  DataHandle, Count, code

  INTEGER i, n

  CALL int_gen_ti_header_c ( hdrbuf, hdrbufsize, itypesize, typesize, &
                             DataHandle, Data, Count, code )
  i = hdrbufsize/itypesize + 1 ;

  CALL int_pack_string ( Element, hdrbuf( i ), n ) ; i = i + n
  CALL int_pack_string ( VarName, hdrbuf( i ), n ) ; i = i + n
  hdrbufsize = i * itypesize + hdrbufsize 
  hdrbuf(1) = hdrbufsize
  RETURN
END SUBROUTINE int_gen_ti_header_integer_varna

SUBROUTINE int_gen_ti_header_real( hdrbuf, hdrbufsize, itypesize, typesize, &
                                   DataHandle, Element, Data, Count, code )





  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize, typesize
  CHARACTER*(*), INTENT(INOUT) ::  Element
  REAL, INTENT(IN)             ::  Data(*)
  INTEGER, INTENT(IN)          ::  DataHandle, Count, code

  INTEGER i, n

  CALL int_gen_ti_header_c ( hdrbuf, hdrbufsize, itypesize, typesize, &
                             DataHandle, Data, Count, code )
  i = hdrbufsize/itypesize + 1 ;

  CALL int_pack_string ( Element, hdrbuf( i ), n ) ; i = i + n
  hdrbufsize = n * itypesize + hdrbufsize 
  hdrbuf(1) = hdrbufsize
  RETURN
END SUBROUTINE int_gen_ti_header_real

SUBROUTINE int_get_ti_header_integer_varna( hdrbuf, hdrbufsize, itypesize, typesize, &
                              DataHandle, Element, VarName, Data, Count, code)






  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize, typesize
  CHARACTER*(*), INTENT(INOUT) ::  Element, VarName
  INTEGER, INTENT(OUT)         ::  Data(*)
  INTEGER, INTENT(OUT)         ::  DataHandle, Count, code

  INTEGER i, n


  CALL int_get_ti_header_c ( hdrbuf, hdrbufsize, n, itypesize, typesize, &
                           DataHandle, Data, Count, code )
  i = n/itypesize + 1
  CALL int_unpack_string ( Element, hdrbuf( i ), n ) ; i=i+n;
  CALL int_unpack_string ( VarName, hdrbuf( i ), n ) ; i = i + n


  hdrbufsize = hdrbuf(1)
  RETURN
END SUBROUTINE int_get_ti_header_integer_varna

SUBROUTINE int_get_ti_header_integer( hdrbuf, hdrbufsize, itypesize, typesize, &
                              DataHandle, Element, Data, Count, code )






  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize, typesize
  CHARACTER*(*), INTENT(INOUT) ::  Element
  INTEGER, INTENT(OUT)         ::  Data(*)
  INTEGER, INTENT(OUT)         ::  DataHandle, Count, code

  INTEGER i, n


  CALL int_get_ti_header_c ( hdrbuf, hdrbufsize, n, itypesize, typesize, &
                           DataHandle, Data, Count, code )
  i = 1 
  CALL int_unpack_string ( Element, hdrbuf( n/itypesize + 1 ), n ) ;

  hdrbufsize = hdrbuf(1)
  RETURN
END SUBROUTINE int_get_ti_header_integer

SUBROUTINE int_get_ti_header_real( hdrbuf, hdrbufsize, itypesize, typesize, &
                              DataHandle, Element, Data, Count, code )






  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize, typesize
  CHARACTER*(*), INTENT(INOUT) ::  Element
  REAL, INTENT(OUT)            ::  Data(*)
  INTEGER, INTENT(OUT)         ::  DataHandle, Count, code

  INTEGER i, n


  CALL int_get_ti_header_c ( hdrbuf, hdrbufsize, n, itypesize, typesize, &
                           DataHandle, Data, Count, code )
  i = 1
  CALL int_unpack_string ( Element, hdrbuf( n/itypesize + 1 ), n ) ;

  hdrbufsize = hdrbuf(1)
  RETURN
END SUBROUTINE int_get_ti_header_real



SUBROUTINE int_gen_ti_header_char( hdrbuf, hdrbufsize, itypesize, &
                              DataHandle, Element, VarName, Data, code )

































  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize
  CHARACTER*(*), INTENT(IN)    :: Element, Data, VarName
  INTEGER, INTENT(IN)          ::  DataHandle, code

  INTEGER                      ::  DummyData
  INTEGER i, n, Count, DummyCount

  DummyCount = 0
  CALL int_gen_ti_header_c ( hdrbuf, hdrbufsize, itypesize, 1, &
                             DataHandle, DummyData, DummyCount, code )
  i = hdrbufsize/itypesize+1 ;
  CALL int_pack_string ( Element, hdrbuf( i ), n ) ; i = i + n
  CALL int_pack_string ( Data   , hdrbuf( i ), n ) ; i = i + n
  CALL int_pack_string ( VarName   , hdrbuf( i ), n ) ; i = i + n
  hdrbufsize = (i-1) * itypesize + hdrbufsize 
  hdrbuf(1) = hdrbufsize
  RETURN
END SUBROUTINE int_gen_ti_header_char

SUBROUTINE int_get_ti_header_char( hdrbuf, hdrbufsize, itypesize, &
                              DataHandle, Element, VarName, Data, code )






  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize
  CHARACTER*(*), INTENT(INOUT) ::  Element, Data, VarName
  INTEGER, INTENT(OUT)         ::  DataHandle, code

  INTEGER i, n, DummyCount, typesize
  CHARACTER * 132  dummyData

  CALL int_get_ti_header_c ( hdrbuf, hdrbufsize, n, itypesize, typesize, &
                           DataHandle, dummyData, DummyCount, code )
  i = n/itypesize+1 ;
  CALL int_unpack_string ( Element, hdrbuf( i ), n ) ; i = i + n
  CALL int_unpack_string ( Data   , hdrbuf( i ), n ) ; i = i + n
  CALL int_unpack_string ( VarName  , hdrbuf( i ), n ) ; i = i + n
  hdrbufsize = hdrbuf(1)

  RETURN
END SUBROUTINE int_get_ti_header_char




SUBROUTINE int_gen_td_header_char( hdrbuf, hdrbufsize, itypesize, &
                              DataHandle, DateStr, Element, Data, code )































  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize
  CHARACTER*(*), INTENT(INOUT) ::  DateStr, Element, Data
  INTEGER, INTENT(IN)          ::  DataHandle, code

  INTEGER i, n, DummyCount, DummyData

  DummyCount = 0

  CALL int_gen_ti_header_c ( hdrbuf, hdrbufsize, itypesize, 1, &
                           DataHandle, DummyData, DummyCount, code )
  i = hdrbufsize/itypesize + 1 ;
  CALL int_pack_string ( Element, hdrbuf( i ), n ) ; i = i + n
  CALL int_pack_string ( DateStr, hdrbuf( i ), n ) ; i = i + n
  CALL int_pack_string ( Data   , hdrbuf( i ), n ) ; i = i + n
  hdrbufsize = (i-1) * itypesize + hdrbufsize 
  hdrbuf(1) = hdrbufsize
  RETURN
END SUBROUTINE int_gen_td_header_char

SUBROUTINE int_get_td_header_char( hdrbuf, hdrbufsize, itypesize, &
                              DataHandle, DateStr, Element, Data, code )






  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize
  CHARACTER*(*), INTENT(INOUT) ::  DateStr, Element, Data
  INTEGER, INTENT(OUT)         ::  DataHandle, code

  INTEGER i, n, Count, typesize


  CALL int_get_ti_header_c ( hdrbuf, hdrbufsize, n, itypesize, typesize, &
                           DataHandle, Data, Count, code )
  i = n/itypesize + 1 ;
  CALL int_unpack_string ( Element, hdrbuf( i ), n ) ; i = i + n ;
  CALL int_unpack_string ( DateStr, hdrbuf( i ), n ) ; i = i + n ;
  CALL int_unpack_string ( Data   , hdrbuf( i ), n ) ; i = i + n ;
  hdrbufsize = hdrbuf(1)
  RETURN
END SUBROUTINE int_get_td_header_char

SUBROUTINE int_gen_td_header_integer( hdrbuf, hdrbufsize, itypesize, typesize, &
                                      DataHandle, DateStr, Element, Data, Count, code )
































  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize, typesize
  CHARACTER*(*), INTENT(INOUT) ::  DateStr, Element
  INTEGER, INTENT(IN)          ::  Data(*)
  INTEGER, INTENT(IN)          ::  DataHandle, Count, code

  INTEGER i, n


  CALL int_gen_ti_header_c ( hdrbuf, hdrbufsize, itypesize, typesize, &
                           DataHandle, Data, Count, code )
  i = hdrbufsize/itypesize + 1 ;
  CALL int_pack_string ( DateStr, hdrbuf( i ), n ) ; i = i + n
  CALL int_pack_string ( Element, hdrbuf( i ), n ) ; i = i + n
  hdrbufsize = (i-1) * itypesize + hdrbufsize 
  hdrbuf(1) = hdrbufsize
  RETURN
END SUBROUTINE int_gen_td_header_integer

SUBROUTINE int_gen_td_header_real( hdrbuf, hdrbufsize, itypesize, typesize, &
                                   DataHandle, DateStr, Element, Data, Count, code )





  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize, typesize
  CHARACTER*(*), INTENT(INOUT) ::  DateStr, Element
  REAL, INTENT(IN)             ::  Data(*)
  INTEGER, INTENT(IN)          ::  DataHandle, Count, code

  INTEGER i, n


  CALL int_gen_ti_header_c ( hdrbuf, hdrbufsize, itypesize, typesize, &
                           DataHandle, Data, Count, code )
  i = hdrbufsize/itypesize + 1 ;
  CALL int_pack_string ( DateStr, hdrbuf( i ), n ) ; i = i + n
  CALL int_pack_string ( Element, hdrbuf( i ), n ) ; i = i + n
  hdrbufsize = (i-1) * itypesize + hdrbufsize 
  hdrbuf(1) = hdrbufsize
  RETURN
END SUBROUTINE int_gen_td_header_real

SUBROUTINE int_get_td_header_integer( hdrbuf, hdrbufsize, itypesize, typesize, &
                              DataHandle, DateStr, Element, Data, Count, code )






  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize, typesize
  CHARACTER*(*), INTENT(INOUT) ::  DateStr, Element
  INTEGER, INTENT(OUT)         ::  Data(*)
  INTEGER, INTENT(OUT)         ::  DataHandle, Count, code

  INTEGER i, n


  CALL int_get_ti_header_c ( hdrbuf, hdrbufsize, n, itypesize, typesize, &
                           DataHandle, Data, Count, code )
  i = n/itypesize + 1 ;
  CALL int_unpack_string ( DateStr, hdrbuf( i ), n ) ; i = i + n ;
  CALL int_unpack_string ( Element, hdrbuf( i ), n ) ; i = i + n ;
  hdrbufsize = hdrbuf(1)
  RETURN
END SUBROUTINE int_get_td_header_integer

SUBROUTINE int_get_td_header_real( hdrbuf, hdrbufsize, itypesize, typesize, &
                              DataHandle, DateStr, Element, Data, Count, code )






  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT)       ::  hdrbuf(*)
  INTEGER, INTENT(OUT)         ::  hdrbufsize
  INTEGER, INTENT(IN)          ::  itypesize, typesize
  CHARACTER*(*), INTENT(INOUT) ::  DateStr, Element
  REAL , INTENT(OUT)           ::  Data(*)
  INTEGER, INTENT(OUT)         ::  DataHandle, Count, code

  INTEGER i, n


  CALL int_get_ti_header_c ( hdrbuf, hdrbufsize, n, itypesize, typesize, &
                           DataHandle, Data, Count, code )
  i = n/itypesize + 1 ;
  CALL int_unpack_string ( DateStr, hdrbuf( i ), n ) ; i = i + n ;
  CALL int_unpack_string ( Element, hdrbuf( i ), n ) ; i = i + n ;
  hdrbufsize = hdrbuf(1)
  RETURN
END SUBROUTINE int_get_td_header_real



SUBROUTINE int_gen_noop_header ( hdrbuf, hdrbufsize, itypesize )
  IMPLICIT NONE

















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
  INTEGER, INTENT(INOUT) ::  hdrbuf(*)
  INTEGER, INTENT(OUT)   ::  hdrbufsize
  INTEGER, INTENT(INOUT) ::  itypesize

  INTEGER i

  hdrbuf(1) = 0  
  hdrbuf(2) = int_noop
  i = 3
  hdrbufsize = (i-1) * itypesize  
  hdrbuf(1) = hdrbufsize
  RETURN
END SUBROUTINE int_gen_noop_header

SUBROUTINE int_get_noop_header( hdrbuf, hdrbufsize, itypesize )






  IMPLICIT NONE
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
  INTEGER, INTENT(INOUT) ::  hdrbuf(*)
  INTEGER, INTENT(OUT)   ::  hdrbufsize
  INTEGER, INTENT(INOUT) ::  itypesize

  INTEGER i

  hdrbufsize = hdrbuf(1)
  IF ( hdrbuf(2) .NE. int_noop ) THEN
    CALL wrf_error_fatal3("<stdin>",1920,&
"int_get_noop_header: hdrbuf ne int_noop")
  ENDIF
  i = 3
  RETURN
END SUBROUTINE int_get_noop_header



SUBROUTINE int_pack_string ( str, buf, n )
  IMPLICIT NONE






  CHARACTER*(*), INTENT(IN)          :: str
  INTEGER, INTENT(OUT)               :: n    
  INTEGER, INTENT(OUT), DIMENSION(*) :: buf

  INTEGER i

  n = 1
  buf(n) = LEN(TRIM(str))
  n = n+1
  DO i = 1, LEN(TRIM(str))
    buf(n) = ichar(str(i:i))
    n = n+1
  ENDDO
  n = n - 1
END SUBROUTINE int_pack_string

SUBROUTINE int_unpack_string ( str, buf, n )
  IMPLICIT NONE






  CHARACTER*(*), INTENT(OUT)        :: str
  INTEGER, INTENT(OUT)              :: n       
  INTEGER, INTENT(IN), DIMENSION(*) :: buf

  INTEGER i
  INTEGER strlen

  strlen = buf(1)
  str = ""
  DO i = 1, strlen
    str(i:i) = char(buf(i+1))
  ENDDO
  n = strlen + 1
END SUBROUTINE int_unpack_string

END MODULE module_internal_header_util

