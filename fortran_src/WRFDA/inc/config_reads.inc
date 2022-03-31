!STARTOFREGISTRYGENERATEDINCLUDE 'inc/config_reads.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
! Contains namelist statements for module_config.F.
#ifndef NAMELIST_READ_UNIT
#  define NAMELIST_READ_UNIT nml_read_unit
#endif
#ifndef NAMELIST_WRITE_UNIT
#  define NAMELIST_WRITE_UNIT nml_write_unit
#endif
!
 nml_read_error = .FALSE.
 NML_LOOP : DO i=1,35
    REWIND ( UNIT = NAMELIST_READ_UNIT )
    SELECT CASE ( i )
       CASE ( 1 ) 
          nml_name = "time_control"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = time_control , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = time_control )
#endif
          CYCLE NML_LOOP
       CASE ( 2 ) 
          nml_name = "pio_control"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = pio_control , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = pio_control )
#endif
          CYCLE NML_LOOP
       CASE ( 3 ) 
          nml_name = "dfi_control"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = dfi_control , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = dfi_control )
#endif
          CYCLE NML_LOOP
       CASE ( 4 ) 
          nml_name = "domains"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = domains , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = domains )
#endif
          CYCLE NML_LOOP
       CASE ( 5 ) 
          nml_name = "tc"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = tc , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = tc )
#endif
          CYCLE NML_LOOP
       CASE ( 6 ) 
          nml_name = "physics"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = physics , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = physics )
#endif
          CYCLE NML_LOOP
       CASE ( 7 ) 
          nml_name = "noah_mp"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = noah_mp , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = noah_mp )
#endif
          CYCLE NML_LOOP
       CASE ( 8 ) 
          nml_name = "fdda"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = fdda , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = fdda )
#endif
          CYCLE NML_LOOP
       CASE ( 9 ) 
          nml_name = "scm"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = scm , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = scm )
#endif
          CYCLE NML_LOOP
       CASE ( 10 ) 
          nml_name = "dynamics"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = dynamics , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = dynamics )
#endif
          CYCLE NML_LOOP
       CASE ( 11 ) 
          nml_name = "bdy_control"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = bdy_control , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = bdy_control )
#endif
          CYCLE NML_LOOP
       CASE ( 12 ) 
          nml_name = "grib2"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = grib2 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = grib2 )
#endif
          CYCLE NML_LOOP
       CASE ( 13 ) 
          nml_name = "wrfvar1"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar1 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar1 )
#endif
          CYCLE NML_LOOP
       CASE ( 14 ) 
          nml_name = "wrfvar2"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar2 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar2 )
#endif
          CYCLE NML_LOOP
       CASE ( 15 ) 
          nml_name = "wrfvar3"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar3 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar3 )
#endif
          CYCLE NML_LOOP
       CASE ( 16 ) 
          nml_name = "wrfvar4"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar4 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar4 )
#endif
          CYCLE NML_LOOP
       CASE ( 17 ) 
          nml_name = "wrfvar5"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar5 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar5 )
#endif
          CYCLE NML_LOOP
       CASE ( 18 ) 
          nml_name = "wrfvar6"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar6 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar6 )
#endif
          CYCLE NML_LOOP
       CASE ( 19 ) 
          nml_name = "wrfvar7"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar7 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar7 )
#endif
          CYCLE NML_LOOP
       CASE ( 20 ) 
          nml_name = "wrfvar8"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar8 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar8 )
#endif
          CYCLE NML_LOOP
       CASE ( 21 ) 
          nml_name = "wrfvar9"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar9 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar9 )
#endif
          CYCLE NML_LOOP
       CASE ( 22 ) 
          nml_name = "wrfvar10"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar10 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar10 )
#endif
          CYCLE NML_LOOP
       CASE ( 23 ) 
          nml_name = "wrfvar11"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar11 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar11 )
#endif
          CYCLE NML_LOOP
       CASE ( 24 ) 
          nml_name = "wrfvar12"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar12 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar12 )
#endif
          CYCLE NML_LOOP
       CASE ( 25 ) 
          nml_name = "wrfvar13"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar13 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar13 )
#endif
          CYCLE NML_LOOP
       CASE ( 26 ) 
          nml_name = "wrfvar14"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar14 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar14 )
#endif
          CYCLE NML_LOOP
       CASE ( 27 ) 
          nml_name = "wrfvar15"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar15 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar15 )
#endif
          CYCLE NML_LOOP
       CASE ( 28 ) 
          nml_name = "wrfvar16"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar16 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar16 )
#endif
          CYCLE NML_LOOP
       CASE ( 29 ) 
          nml_name = "wrfvar17"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar17 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar17 )
#endif
          CYCLE NML_LOOP
       CASE ( 30 ) 
          nml_name = "wrfvar18"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar18 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar18 )
#endif
          CYCLE NML_LOOP
       CASE ( 31 ) 
          nml_name = "wrfvar19"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar19 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar19 )
#endif
          CYCLE NML_LOOP
       CASE ( 32 ) 
          nml_name = "wrfvar20"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar20 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar20 )
#endif
          CYCLE NML_LOOP
       CASE ( 33 ) 
          nml_name = "wrfvar21"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar21 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar21 )
#endif
          CYCLE NML_LOOP
       CASE ( 34 ) 
          nml_name = "wrfvar22"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = wrfvar22 , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = wrfvar22 )
#endif
          CYCLE NML_LOOP
       CASE ( 35 ) 
          nml_name = "perturbation"
          READ   ( UNIT = NAMELIST_READ_UNIT , NML = perturbation , ERR=9201, END=9202 )
#ifndef NO_NAMELIST_PRINT
          WRITE ( UNIT = NAMELIST_WRITE_UNIT, NML = perturbation )
#endif
          CYCLE NML_LOOP
    END SELECT
9201 CALL wrf_message("  ------ ERROR while reading namelist "//TRIM(nml_name)//" ------")
    nml_read_error = .TRUE.
    CALL wrf_alt_nml_obsolete(nml_read_unit, TRIM(nml_name))
    CYCLE NML_LOOP
9202 CALL wrf_debug(1,"Namelist "//TRIM(nml_name)//" not found in namelist.input.")
     CALL wrf_debug(1," --> Using registry defaults for variables in "//TRIM(nml_name))
 END DO NML_LOOP
 
 IF ( nml_read_error ) CALL wrf_error_fatal("ERRORS while reading one or more namelists from namelist.input.")
!ENDOFREGISTRYGENERATEDINCLUDE
