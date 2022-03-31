












program gen_be_etkf

















   use da_control, only : trace_use, stdout,filename_len

   use da_etkf, only : da_solve_etkf, da_solve_letkf, rand_filter
   use da_reporting, only : da_error, message

   implicit none

!     NetCDF-3.
!
! netcdf version 3 fortran interface:
!

!
! external netcdf data types:
!
      integer nf_byte
      integer nf_int1
      integer nf_char
      integer nf_short
      integer nf_int2
      integer nf_int
      integer nf_float
      integer nf_real
      integer nf_double
      integer nf_ubyte
      integer nf_ushort
      integer nf_uint
      integer nf_int64
      integer nf_uint64

      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)

!
! default fill values:
!
      integer           nf_fill_byte
      integer           nf_fill_int1
      integer           nf_fill_char
      integer           nf_fill_short
      integer           nf_fill_int2
      integer           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double

      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690e+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690d+36)

!
! mode flags for opening and creating a netcdf dataset:
!
      integer nf_nowrite
      integer nf_write
      integer nf_clobber
      integer nf_noclobber
      integer nf_fill
      integer nf_nofill
      integer nf_lock
      integer nf_share
      integer nf_64bit_offset
      integer nf_64bit_data
      integer nf_cdf5
      integer nf_sizehint_default
      integer nf_align_chunk
      integer nf_format_classic
      integer nf_format_64bit
      integer nf_format_64bit_offset
      integer nf_format_64bit_data
      integer nf_format_cdf5
      integer nf_diskless
      integer nf_mmap

      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)

!
! size argument for defining an unlimited dimension:
!
      integer nf_unlimited
      parameter (nf_unlimited = 0)

!
! global attribute id:
!
      integer nf_global
      parameter (nf_global = 0)

!
! implementation limits:
!
      integer nf_max_dims
      integer nf_max_attrs
      integer nf_max_vars
      integer nf_max_name
      integer nf_max_var_dims

      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)

!
! error codes:
!
      integer nf_noerr
      integer nf_ebadid
      integer nf_eexist
      integer nf_einval
      integer nf_eperm
      integer nf_enotindefine
      integer nf_eindefine
      integer nf_einvalcoords
      integer nf_emaxdims
      integer nf_enameinuse
      integer nf_enotatt
      integer nf_emaxatts
      integer nf_ebadtype
      integer nf_ebaddim
      integer nf_eunlimpos
      integer nf_emaxvars
      integer nf_enotvar
      integer nf_eglobal
      integer nf_enotnc
      integer nf_ests
      integer nf_emaxname
      integer nf_eunlimit
      integer nf_enorecvars
      integer nf_echar
      integer nf_eedge
      integer nf_estride
      integer nf_ebadname
      integer nf_erange
      integer nf_enomem
      integer nf_evarsize
      integer nf_edimsize
      integer nf_etrunc

      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
!
! error handling modes:
!
      integer  nf_fatal
      integer nf_verbose

      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)

!
! miscellaneous routines:
!
      character*80   nf_inq_libvers
      external       nf_inq_libvers

      character*80   nf_strerror
!                         (integer             ncerr)
      external       nf_strerror

      logical        nf_issyserr
!                         (integer             ncerr)
      external       nf_issyserr

!
! control routines:
!
      integer         nf_inq_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_inq_base_pe

      integer         nf_set_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_set_base_pe

      integer         nf_create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             ncid)
      external        nf_create

      integer         nf__create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create

      integer         nf__create_mp
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create_mp

      integer         nf_open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             ncid)
      external        nf_open

      integer         nf__open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open

      integer         nf__open_mp
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open_mp

      integer         nf_set_fill
!                         (integer             ncid,
!                          integer             fillmode,
!                          integer             old_mode)
      external        nf_set_fill

      integer         nf_set_default_format
!                          (integer             format,
!                          integer             old_format)
      external        nf_set_default_format

      integer         nf_redef
!                         (integer             ncid)
      external        nf_redef

      integer         nf_enddef
!                         (integer             ncid)
      external        nf_enddef

      integer         nf__enddef
!                         (integer             ncid,
!                          integer             h_minfree,
!                          integer             v_align,
!                          integer             v_minfree,
!                          integer             r_align)
      external        nf__enddef

      integer         nf_sync
!                         (integer             ncid)
      external        nf_sync

      integer         nf_abort
!                         (integer             ncid)
      external        nf_abort

      integer         nf_close
!                         (integer             ncid)
      external        nf_close

      integer         nf_delete
!                         (character*(*)       ncid)
      external        nf_delete

!
! general inquiry routines:
!

      integer         nf_inq
!                         (integer             ncid,
!                          integer             ndims,
!                          integer             nvars,
!                          integer             ngatts,
!                          integer             unlimdimid)
      external        nf_inq

! new inquire path

      integer nf_inq_path
      external nf_inq_path

      integer         nf_inq_ndims
!                         (integer             ncid,
!                          integer             ndims)
      external        nf_inq_ndims

      integer         nf_inq_nvars
!                         (integer             ncid,
!                          integer             nvars)
      external        nf_inq_nvars

      integer         nf_inq_natts
!                         (integer             ncid,
!                          integer             ngatts)
      external        nf_inq_natts

      integer         nf_inq_unlimdim
!                         (integer             ncid,
!                          integer             unlimdimid)
      external        nf_inq_unlimdim

      integer         nf_inq_format
!                         (integer             ncid,
!                          integer             format)
      external        nf_inq_format

!
! dimension routines:
!

      integer         nf_def_dim
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             len,
!                          integer             dimid)
      external        nf_def_dim

      integer         nf_inq_dimid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             dimid)
      external        nf_inq_dimid

      integer         nf_inq_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_dim

      integer         nf_inq_dimname
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_inq_dimname

      integer         nf_inq_dimlen
!                         (integer             ncid,
!                          integer             dimid,
!                          integer             len)
      external        nf_inq_dimlen

      integer         nf_rename_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_rename_dim

!
! general attribute routines:
!

      integer         nf_inq_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len)
      external        nf_inq_att

      integer         nf_inq_attid
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             attnum)
      external        nf_inq_attid

      integer         nf_inq_atttype
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype)
      external        nf_inq_atttype

      integer         nf_inq_attlen
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_attlen

      integer         nf_inq_attname
!                         (integer             ncid,
!                          integer             varid,
!                          integer             attnum,
!                          character(*)        name)
      external        nf_inq_attname

      integer         nf_copy_att
!                         (integer             ncid_in,
!                          integer             varid_in,
!                          character(*)        name,
!                          integer             ncid_out,
!                          integer             varid_out)
      external        nf_copy_att

      integer         nf_rename_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        curname,
!                          character(*)        newname)
      external        nf_rename_att

      integer         nf_del_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_del_att

!
! attribute put/get routines:
!

      integer         nf_put_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len,
!                          character(*)        text)
      external        nf_put_att_text

      integer         nf_get_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          character(*)        text)
      external        nf_get_att_text

      integer         nf_put_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int1_t           i1vals(1))
      external        nf_put_att_int1

      integer         nf_get_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int1_t           i1vals(1))
      external        nf_get_att_int1

      integer         nf_put_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int2_t           i2vals(1))
      external        nf_put_att_int2

      integer         nf_get_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int2_t           i2vals(1))
      external        nf_get_att_int2

      integer         nf_put_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          integer             ivals(1))
      external        nf_put_att_int

      integer         nf_get_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             ivals(1))
      external        nf_get_att_int

      integer         nf_put_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int8_t           i8vals(1))
      external        nf_put_att_int64

      integer         nf_get_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int8_t           i8vals(1))
      external        nf_get_att_int64

      integer         nf_put_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          real                rvals(1))
      external        nf_put_att_real

      integer         nf_get_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          real                rvals(1))
      external        nf_get_att_real

      integer         nf_put_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          double              dvals(1))
      external        nf_put_att_double

      integer         nf_get_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          double              dvals(1))
      external        nf_get_att_double

!
! general variable routines:
!

      integer         nf_def_var
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             varid)
      external        nf_def_var

      integer         nf_inq_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             natts)
      external        nf_inq_var

      integer         nf_inq_varid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             varid)
      external        nf_inq_varid

      integer         nf_inq_varname
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_inq_varname

      integer         nf_inq_vartype
!                         (integer             ncid,
!                          integer             varid,
!                          integer             xtype)
      external        nf_inq_vartype

      integer         nf_inq_varndims
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ndims)
      external        nf_inq_varndims

      integer         nf_inq_vardimid
!                         (integer             ncid,
!                          integer             varid,
!                          integer             dimids(1))
      external        nf_inq_vardimid

      integer         nf_inq_varnatts
!                         (integer             ncid,
!                          integer             varid,
!                          integer             natts)
      external        nf_inq_varnatts

      integer         nf_rename_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_rename_var

      integer         nf_copy_var
!                         (integer             ncid_in,
!                          integer             varid,
!                          integer             ncid_out)
      external        nf_copy_var

!
! entire variable put/get routines:
!

      integer         nf_put_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_put_var_text

      integer         nf_get_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_get_var_text

      integer         nf_put_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_put_var_int1

      integer         nf_get_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_get_var_int1

      integer         nf_put_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_put_var_int2

      integer         nf_get_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_get_var_int2

      integer         nf_put_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_put_var_int

      integer         nf_get_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_get_var_int

      integer         nf_put_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_put_var_real

      integer         nf_get_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_get_var_real

      integer         nf_put_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_put_var_double

      integer         nf_get_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_get_var_double

!
! single variable put/get routines:
!

      integer         nf_put_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_put_var1_text

      integer         nf_get_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_get_var1_text

      integer         nf_put_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_put_var1_int1

      integer         nf_get_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_get_var1_int1

      integer         nf_put_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_put_var1_int2

      integer         nf_get_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_get_var1_int2

      integer         nf_put_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_put_var1_int

      integer         nf_get_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_get_var1_int

      integer         nf_put_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_put_var1_real

      integer         nf_get_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_get_var1_real

      integer         nf_put_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_put_var1_double

      integer         nf_get_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_get_var1_double

!
! variable array put/get routines:
!

      integer         nf_put_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_put_vara_text

      integer         nf_get_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_get_vara_text

      integer         nf_put_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vara_int1

      integer         nf_get_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vara_int1

      integer         nf_put_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vara_int2

      integer         nf_get_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vara_int2

      integer         nf_put_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_put_vara_int

      integer         nf_get_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_get_vara_int

      integer         nf_put_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_put_vara_real

      integer         nf_get_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_get_vara_real

      integer         nf_put_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vara_double

      integer         nf_get_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vara_double

!
! strided variable put/get routines:
!

      integer         nf_put_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_put_vars_text

      integer         nf_get_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_get_vars_text

      integer         nf_put_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vars_int1

      integer         nf_get_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vars_int1

      integer         nf_put_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vars_int2

      integer         nf_get_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vars_int2

      integer         nf_put_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_put_vars_int

      integer         nf_get_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_get_vars_int

      integer         nf_put_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_put_vars_real

      integer         nf_get_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_get_vars_real

      integer         nf_put_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vars_double

      integer         nf_get_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vars_double

!
! mapped variable put/get routines:
!

      integer         nf_put_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_put_varm_text

      integer         nf_get_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_get_varm_text

      integer         nf_put_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_varm_int1

      integer         nf_get_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_varm_int1

      integer         nf_put_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_varm_int2

      integer         nf_get_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_varm_int2

      integer         nf_put_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_put_varm_int

      integer         nf_get_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_get_varm_int

      integer         nf_put_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_put_varm_real

      integer         nf_get_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_get_varm_real

      integer         nf_put_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_put_varm_double

      integer         nf_get_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_get_varm_double

!     64-bit int functions.
      integer nf_put_var1_int64
      external nf_put_var1_int64
      integer nf_put_vara_int64
      external nf_put_vara_int64
      integer nf_put_vars_int64
      external nf_put_vars_int64
      integer nf_put_varm_int64
      external nf_put_varm_int64
      integer nf_put_var_int64
      external nf_put_var_int64
      integer nf_get_var1_int64
      external nf_get_var1_int64
      integer nf_get_vara_int64
      external nf_get_vara_int64
      integer nf_get_vars_int64
      external nf_get_vars_int64
      integer nf_get_varm_int64
      external nf_get_varm_int64
      integer nf_get_var_int64
      external nf_get_var_int64


!     NetCDF-4.
!     This is part of netCDF-4. Copyright 2006, UCAR, See COPYRIGHT
!     file for distribution information.

!     Netcdf version 4 fortran interface.

!     $Id: netcdf4.inc,v 1.28 2010/05/25 13:53:02 ed Exp $

!     New netCDF-4 types.
      integer nf_string
      integer nf_vlen
      integer nf_opaque
      integer nf_enum
      integer nf_compound

      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)

!     New netCDF-4 fill values.
      integer           nf_fill_ubyte
      integer           nf_fill_ushort
!      real              nf_fill_uint
!      real              nf_fill_int64
!      real              nf_fill_uint64
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)

!     New constants.
      integer nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)

      integer nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)

      integer nf_netcdf4
      parameter (nf_netcdf4 = 4096)

      integer nf_classic_model
      parameter (nf_classic_model = 256)

      integer nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)

      integer nf_endian_native
      parameter (nf_endian_native = 0)
      integer nf_endian_little
      parameter (nf_endian_little = 1)
      integer nf_endian_big
      parameter (nf_endian_big = 2)

!     For NF_DEF_VAR_CHUNKING
      integer nf_chunked
      parameter (nf_chunked = 0)
      integer nf_contiguous
      parameter (nf_contiguous = 1)
      integer nf_compact
      parameter (nf_compact = 2)

!     For NF_DEF_VAR_FLETCHER32
      integer nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer nf_fletcher32
      parameter (nf_fletcher32 = 1)

!     For NF_DEF_VAR_DEFLATE
      integer nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer nf_shuffle
      parameter (nf_shuffle = 1)

!     For NF_DEF_VAR_SZIP
      integer nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)

!     For parallel I/O.
      integer nf_mpiio      
      parameter (nf_mpiio = 8192)
      integer nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer nf_pnetcdf
      parameter (nf_pnetcdf = 32768)

!     For NF_VAR_PAR_ACCESS.
      integer nf_independent
      parameter (nf_independent = 0)
      integer nf_collective
      parameter (nf_collective = 1)

!     New error codes.
      integer nf_ehdferr        ! Error at 1 layer. 
      parameter (nf_ehdferr = -101)
      integer nf_ecantread      ! Can't read. 
      parameter (nf_ecantread = -102)
      integer nf_ecantwrite     ! Can't write. 
      parameter (nf_ecantwrite = -103)
      integer nf_ecantcreate    ! Can't create. 
      parameter (nf_ecantcreate = -104)
      integer nf_efilemeta      ! Problem with file metadata. 
      parameter (nf_efilemeta = -105)
      integer nf_edimmeta       ! Problem with dimension metadata. 
      parameter (nf_edimmeta = -106)
      integer nf_eattmeta       ! Problem with attribute metadata. 
      parameter (nf_eattmeta = -107)
      integer nf_evarmeta       ! Problem with variable metadata. 
      parameter (nf_evarmeta = -108)
      integer nf_enocompound    ! Not a compound type. 
      parameter (nf_enocompound = -109)
      integer nf_eattexists     ! Attribute already exists. 
      parameter (nf_eattexists = -110)
      integer nf_enotnc4        ! Attempting netcdf-4 operation on netcdf-3 file.   
      parameter (nf_enotnc4 = -111)
      integer nf_estrictnc3     ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
      parameter (nf_estrictnc3 = -112)
      integer nf_enotnc3        ! Attempting netcdf-3 operation on netcdf-4 file.   
      parameter (nf_enotnc3 = -113)
      integer nf_enopar         ! Parallel operation on file opened for non-parallel access.   
      parameter (nf_enopar = -114)
      integer nf_eparinit       ! Error initializing for parallel access.   
      parameter (nf_eparinit = -115)
      integer nf_ebadgrpid      ! Bad group ID.   
      parameter (nf_ebadgrpid = -116)
      integer nf_ebadtypid      ! Bad type ID.   
      parameter (nf_ebadtypid = -117)
      integer nf_etypdefined    ! Type has already been defined and may not be edited. 
      parameter (nf_etypdefined = -118)
      integer nf_ebadfield      ! Bad field ID.   
      parameter (nf_ebadfield = -119)
      integer nf_ebadclass      ! Bad class.   
      parameter (nf_ebadclass = -120)
      integer nf_emaptype       ! Mapped access for atomic types only.   
      parameter (nf_emaptype = -121)
      integer nf_elatefill      ! Attempt to define fill value when data already exists. 
      parameter (nf_elatefill = -122)
      integer nf_elatedef       ! Attempt to define var properties, like deflate, after enddef. 
      parameter (nf_elatedef = -123)
      integer nf_edimscale      ! Probem with 1 dimscales. 
      parameter (nf_edimscale = -124)
      integer nf_enogrp       ! No group found.
      parameter (nf_enogrp = -125)


!     New functions.

!     Parallel I/O.
      integer nf_create_par
      external nf_create_par

      integer nf_open_par
      external nf_open_par

      integer nf_var_par_access
      external nf_var_par_access

!     Functions to handle groups.
      integer nf_inq_ncid
      external nf_inq_ncid

      integer nf_inq_grps
      external nf_inq_grps

      integer nf_inq_grpname
      external nf_inq_grpname

      integer nf_inq_grpname_full
      external nf_inq_grpname_full

      integer nf_inq_grpname_len
      external nf_inq_grpname_len

      integer nf_inq_grp_parent
      external nf_inq_grp_parent

      integer nf_inq_grp_ncid
      external nf_inq_grp_ncid

      integer nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid

      integer nf_inq_varids
      external nf_inq_varids

      integer nf_inq_dimids
      external nf_inq_dimids

      integer nf_def_grp
      external nf_def_grp

!     New rename grp function

      integer nf_rename_grp
      external nf_rename_grp

!     New options for netCDF variables.
      integer nf_def_var_deflate
      external nf_def_var_deflate

      integer nf_inq_var_deflate
      external nf_inq_var_deflate

      integer nf_def_var_szip
      external nf_def_var_szip

      integer nf_inq_var_szip
      external nf_inq_var_szip

      integer nf_def_var_fletcher32
      external nf_def_var_fletcher32

      integer nf_inq_var_fletcher32
      external nf_inq_var_fletcher32

      integer nf_def_var_chunking
      external nf_def_var_chunking

      integer nf_inq_var_chunking
      external nf_inq_var_chunking

      integer nf_def_var_fill
      external nf_def_var_fill

      integer nf_inq_var_fill
      external nf_inq_var_fill

      integer nf_def_var_endian
      external nf_def_var_endian

      integer nf_inq_var_endian
      external nf_inq_var_endian

      integer nf_def_var_filter
      external nf_def_var_filter

      integer nf_inq_var_filter
      external nf_inq_var_filter

!     User defined types.
      integer nf_inq_typeids
      external nf_inq_typeids

      integer nf_inq_typeid
      external nf_inq_typeid

      integer nf_inq_type
      external nf_inq_type

      integer nf_inq_user_type
      external nf_inq_user_type

!     User defined types - compound types.
      integer nf_def_compound
      external nf_def_compound

      integer nf_insert_compound
      external nf_insert_compound

      integer nf_insert_array_compound
      external nf_insert_array_compound

      integer nf_inq_compound
      external nf_inq_compound

      integer nf_inq_compound_name
      external nf_inq_compound_name

      integer nf_inq_compound_size
      external nf_inq_compound_size

      integer nf_inq_compound_nfields
      external nf_inq_compound_nfields

      integer nf_inq_compound_field
      external nf_inq_compound_field

      integer nf_inq_compound_fieldname
      external nf_inq_compound_fieldname

      integer nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex

      integer nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset

      integer nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype

      integer nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims

      integer nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes

!     User defined types - variable length arrays.
      integer nf_def_vlen
      external nf_def_vlen

      integer nf_inq_vlen
      external nf_inq_vlen

      integer nf_free_vlen
      external nf_free_vlen

!     User defined types - enums.
      integer nf_def_enum
      external nf_def_enum

      integer nf_insert_enum
      external nf_insert_enum

      integer nf_inq_enum
      external nf_inq_enum

      integer nf_inq_enum_member
      external nf_inq_enum_member

      integer nf_inq_enum_ident
      external nf_inq_enum_ident

!     User defined types - opaque.
      integer nf_def_opaque
      external nf_def_opaque

      integer nf_inq_opaque
      external nf_inq_opaque

!     Write and read attributes of any type, including user defined
!     types.
      integer nf_put_att
      external nf_put_att
      integer nf_get_att
      external nf_get_att

!     Write and read variables of any type, including user defined
!     types.
      integer nf_put_var
      external nf_put_var
      integer nf_put_var1
      external nf_put_var1
      integer nf_put_vara
      external nf_put_vara
      integer nf_put_vars
      external nf_put_vars
      integer nf_get_var
      external nf_get_var
      integer nf_get_var1
      external nf_get_var1
      integer nf_get_vara
      external nf_get_vara
      integer nf_get_vars
      external nf_get_vars

!     For helping F77 users with VLENs.
      integer nf_get_vlen_element
      external nf_get_vlen_element
      integer nf_put_vlen_element
      external nf_put_vlen_element

!     For dealing with file level chunk cache.
      integer nf_set_chunk_cache
      external nf_set_chunk_cache
      integer nf_get_chunk_cache
      external nf_get_chunk_cache

!     For dealing with per variable chunk cache.
      integer nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer nf_get_var_chunk_cache
      external nf_get_var_chunk_cache

!     NetCDF-2.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!      
! functions in the fortran interface
!
      integer nccre
      integer ncopn
      integer ncddef
      integer ncdid
      integer ncvdef
      integer ncvid
      integer nctlen
      integer ncsfil

      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil


      integer ncrdwr
      integer nccreat
      integer ncexcl
      integer ncindef
      integer ncnsync
      integer nchsync
      integer ncndirty
      integer nchdirty
      integer nclink
      integer ncnowrit
      integer ncwrite
      integer ncclob
      integer ncnoclob
      integer ncglobal
      integer ncfill
      integer ncnofill
      integer maxncop
      integer maxncdim
      integer maxncatt
      integer maxncvar
      integer maxncnam
      integer maxvdims
      integer ncnoerr
      integer ncebadid
      integer ncenfile
      integer nceexist
      integer nceinval
      integer nceperm
      integer ncenotin
      integer nceindef
      integer ncecoord
      integer ncemaxds
      integer ncename
      integer ncenoatt
      integer ncemaxat
      integer ncebadty
      integer ncebadd
      integer ncests
      integer nceunlim
      integer ncemaxvs
      integer ncenotvr
      integer nceglob
      integer ncenotnc
      integer ncfoobar
      integer ncsyserr
      integer ncfatal
      integer ncverbos
      integer ncentool


!
! netcdf data types:
!
      integer ncbyte
      integer ncchar
      integer ncshort
      integer nclong
      integer ncfloat
      integer ncdouble

      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)

!     
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!     

!     read/write, 0 => readonly 
      parameter(ncrdwr = 1)
!     in create phase, cleared by ncendef 
      parameter(nccreat = 2)
!     on create destroy existing file 
      parameter(ncexcl = 4)
!     in define mode, cleared by ncendef 
      parameter(ncindef = 8)
!     synchronise numrecs on change (x'10')
      parameter(ncnsync = 16)
!     synchronise whole header on change (x'20')
      parameter(nchsync = 32)
!     numrecs has changed (x'40')
      parameter(ncndirty = 64)  
!     header info has changed (x'80')
      parameter(nchdirty = 128)
!     prefill vars on endef and increase of record, the default behavior
      parameter(ncfill = 0)
!     do not fill vars on endef and increase of record (x'100')
      parameter(ncnofill = 256)
!     isa link (x'8000')
      parameter(nclink = 32768)

!     
!     'mode' arguments for nccreate and ncopen
!     
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)

!     
!     'size' argument to ncdimdef for an unlimited dimension
!     
      integer ncunlim
      parameter(ncunlim = 0)

!     
!     attribute id to put/get a global attribute
!     
      parameter(ncglobal  = 0)

!     
!     advisory maximums:
!     
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
!     not enforced 
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)

!     
!     global netcdf error status variable
!     initialized in error.c
!     

!     no error 
      parameter(ncnoerr = nf_noerr)
!     not a netcdf id 
      parameter(ncebadid = nf_ebadid)
!     too many netcdfs open 
      parameter(ncenfile = -31)   ! nc_syserr
!     netcdf file exists && ncnoclob
      parameter(nceexist = nf_eexist)
!     invalid argument 
      parameter(nceinval = nf_einval)
!     write to read only 
      parameter(nceperm = nf_eperm)
!     operation not allowed in data mode 
      parameter(ncenotin = nf_enotindefine )   
!     operation not allowed in define mode 
      parameter(nceindef = nf_eindefine)   
!     coordinates out of domain 
      parameter(ncecoord = nf_einvalcoords)
!     maxncdims exceeded 
      parameter(ncemaxds = nf_emaxdims)
!     string match to name in use 
      parameter(ncename = nf_enameinuse)   
!     attribute not found 
      parameter(ncenoatt = nf_enotatt)
!     maxncattrs exceeded 
      parameter(ncemaxat = nf_emaxatts)
!     not a netcdf data type 
      parameter(ncebadty = nf_ebadtype)
!     invalid dimension id 
      parameter(ncebadd = nf_ebaddim)
!     ncunlimited in the wrong index 
      parameter(nceunlim = nf_eunlimpos)
!     maxncvars exceeded 
      parameter(ncemaxvs = nf_emaxvars)
!     variable not found 
      parameter(ncenotvr = nf_enotvar)
!     action prohibited on ncglobal varid 
      parameter(nceglob = nf_eglobal)
!     not a netcdf file 
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname) 
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)

!     
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!     
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)

!
!     default fill values.  these must be the same as in the c interface.
!
      integer filbyte
      integer filchar
      integer filshort
      integer fillong
      real filfloat
      doubleprecision fildoub

      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690e+36)
      parameter (fildoub = 9.9692099683868690e+36)

   integer, parameter    :: max_num_vars = 50         
   integer, parameter    :: max_num_dims = 20         
   integer, parameter    :: unit = 100                

   character (len=filename_len)   :: input_file                
   character (len=filename_len)   :: output_file               
   character (len=3)     :: ce                        

   integer               :: num_members               
   integer               :: nv                        
   integer               :: num_obs                   
   integer               :: naccumt1                  
   integer               :: naccumt2                  
   integer               :: nstartaccum1              
   integer               :: nstartaccum2              
   integer               :: nout                      
   integer               :: length                    
   integer               :: rcode                     
   integer               :: cdfid                     
   integer               :: member, v, o, i, j, k, ijkv 
   integer               :: ivtype                    
   integer               :: natts                     

   integer               :: index                     
   integer               :: nijkv                     
   integer               :: iend                      
   real                  :: num_members_inv           
   real                  :: tainflatinput             
   real                  :: rhoinput                  
   real                  :: ds                        

   character(len=10)     :: cv(1:max_num_vars)        
   integer               :: id_var(1:max_num_vars)    
   integer               :: ndims(1:max_num_vars)     
   integer               :: istart(1:max_num_vars)    
   integer               :: dimids(1:max_num_dims)    
   integer               :: one(1:max_num_dims)       
   integer               :: dims(1:max_num_vars,1:max_num_dims)      
   integer               :: dim_prod(1:max_num_dims)  
   character (len=NF_MAX_NAME) :: var_name            

   real (kind=4), allocatable     :: data_r_3d(:,:,:)               
   real (kind=4), allocatable     :: data_r_4d(:,:,:,:)             

   real, pointer         :: xf(:,:)                   
   real, pointer         :: xf_mean(:)                
   real, pointer         :: xf_vari(:)                
   real, pointer         :: y(:,:)                    
   real, pointer         :: sigma_o2(:)               
   real, pointer         :: yo(:)                     
   real, pointer         :: ens_mean(:)               
   real, pointer         :: ens_stdv_pert_prior(:)    
   real, pointer         :: ens_stdv_pert_poste(:)    

   character(len=150)    :: infl_let_file             
   character(len=150)    :: infl_fac_file             
   character(len=150)    :: eigen_val_file            
   character(len=150)    :: inno2_val_file            
   character(len=150)    :: proj2_val_file            
   real                  :: num_members_m1_inv        
   
   real, pointer                 :: apm_yo(:,:)
   real, pointer                 :: apm_y(:,:)
   real, pointer                 :: apm_sigma_o2(:,:)
   character(len=10), pointer    :: apm_obs_type(:,:)
   real, pointer                 :: apm_type(:,:)
   real, pointer                 :: apm_subtype(:,:)
   real, pointer                 :: apm_obslat(:,:)
   real, pointer                 :: apm_obslon(:,:)
   real, pointer                 :: apm_obsprs(:,:)
   real, pointer                 :: apm_obstim(:,:)
   integer, pointer              :: apm_irec(:,:)
   character(len=10), pointer    :: yo_obs_typ(:)
   real, pointer                 :: yo_typ(:)
   real, pointer                 :: yo_subtyp(:)
   real, pointer                 :: yo_lat(:)
   real, pointer                 :: yo_lon(:)
   real, pointer                 :: yo_prs(:)
   real, pointer                 :: yo_tim(:)
   real, pointer                 :: xf_lon(:)
   real, pointer                 :: xf_lat(:)
   real (kind=4), allocatable    :: ens_lon_T(:,:),ens_lon_u(:,:),ens_lon_v(:,:)
   real (kind=4), allocatable    :: ens_lat_T(:,:),ens_lat_u(:,:),ens_lat_v(:,:)
   real (kind=4), allocatable    :: ens_prs_T(:,:,:),ens_pbs_T(:,:,:)
   integer                       :: ityp,idim
   integer                       :: id_lat_T,id_lat_u,id_lat_v
   integer                       :: id_lon_T,id_lon_u,id_lon_v
   integer                       :: id_prs_T,id_pbs_T
   integer                       :: ndim_lat_T,ndim_lat_u,ndim_lat_v
   integer                       :: ndim_lon_T,ndim_lon_u,ndim_lon_v
   integer                       :: ndim_prs_T,ndim_pbs_T
   integer, pointer              :: dim_lat_T(:),dim_lat_u(:),dim_lat_v(:)
   integer, pointer              :: dim_lon_T(:),dim_lon_u(:),dim_lon_v(:)
   integer, pointer              :: dim_prs_T(:),dim_pbs_T(:)

   integer                       :: apm_num_obs
   integer                       :: apm_icnt
   integer                       :: apm_reject
   integer                       :: apm_imem
   logical                       :: infl_fac_TRNK,infl_fac_WG03,infl_fac_WG07,letkf_flg
   logical                       :: infl_fac_BOWL, rand_filt 
   integer, allocatable          :: ijk_idx(:,:,:,:)
   integer, allocatable          :: apl_idx(:,:)
   integer, allocatable          :: apl_ndim(:)
   integer                       :: napl1, napl2, nxs, nys, nzs, ixy, idx, apm_unit
   character(len=150),dimension(:),allocatable    :: filt_ob_etkf              
   integer                       :: rnd_seed, rnd_nobs, nobs_flt, nseed

   logical                       :: etkf_erro_flg, etkf_inno_flg, etkf_wrfda  
   real                          :: etkf_erro_max, etkf_erro_min
   real                          :: etkf_inno_max, etkf_inno_min

   namelist / gen_be_etkf_nl / num_members, nv, cv, &
                               naccumt1, naccumt2, nstartaccum1, nstartaccum2, &
                               nout, tainflatinput, rhoinput, infl_fac_file, &
                               eigen_val_file, inno2_val_file, proj2_val_file, &
                               infl_fac_TRNK, infl_fac_WG03, infl_fac_WG07, &
                               infl_fac_BOWL, letkf_flg, rand_filt, rnd_seed, &
                               rnd_nobs, infl_let_file, etkf_erro_max, etkf_erro_min, &
                               etkf_inno_max, etkf_inno_min, etkf_erro_flg, etkf_inno_flg, &
                               etkf_wrfda
   

   write(stdout,'(/a)')' [1] Initialize information.'


   num_members = 56
   nv = 1
   cv = "U"
   naccumt1 = 0
   naccumt2 = 0
   nstartaccum1 = 0
   nstartaccum2 = 0
   nout = 1 
   tainflatinput = 0.0
   rhoinput = 0.0
   infl_fac_file  = "inflation_factor.dat"
   infl_let_file  = "inflation_letkf.dat"
   eigen_val_file = "eigen_value.dat"
   inno2_val_file = "innovation_value.dat"
   proj2_val_file = "projection_value.dat"

   open(unit=unit, file='gen_be_etkf_nl.nl', &
        form='formatted', status='old', action='read')
   read(unit, gen_be_etkf_nl)
   close(unit)

   write(stdout,'(a,i4)')'   Number of ensemble members = ', num_members
   write(stdout,'(a,i4)')'   Number of prognostic variables = ', nv
   write(stdout,'(50a)')'    List of prognostic variables = ', cv(1:nv)
   write(stdout,'(a,i4)')'   naccumt1 = ', naccumt1
   write(stdout,'(a,i4)')'   naccumt2 = ', naccumt2
   write(stdout,'(a,i4)')'   nstartaccum1 = ', nstartaccum1
   write(stdout,'(a,i4)')'   nstartaccum2 = ', nstartaccum2
   write(stdout,'(a,i4)')'   nout = ', nout
   write(stdout,'(a,f15.5)')'   tainflatinput = ', tainflatinput
   write(stdout,'(a,f15.5)')'   rhoinput = ', rhoinput
   write(stdout,'(150a)')'   infl_fac_file = ',infl_fac_file
   write(stdout,'(150a)')'   infl_let_file = ',infl_let_file
   write(stdout,'(150a)')'   eigen_val_file = ',eigen_val_file
   write(stdout,'(150a)')'   inno2_val_file = ',inno2_val_file
   write(stdout,'(150a)')'   proj2_val_file = ',proj2_val_file
   write(stdout,'(a,l10)')'   infl_fac_TRNK = ',infl_fac_TRNK
   write(stdout,'(a,l10)')'   infl_fac_WG03 = ',infl_fac_WG03
   write(stdout,'(a,l10)')'   infl_fac_WG07 = ',infl_fac_WG07
   write(stdout,'(a,l10)')'   infl_fac_BOWL = ',infl_fac_BOWL
   write(stdout,'(a,l10)')'   letkf_flag = ',letkf_flg
   write(stdout,'(a,i20)')'   rnd_seed = ',rnd_seed
   write(stdout,'(a,i20)')'   rnd_nobs = ',rnd_nobs
   write(stdout,'(a,f15.5)')' etkf_erro_max = ',etkf_erro_max
   write(stdout,'(a,f15.5)')' etkf_erro_min = ',etkf_erro_min
   write(stdout,'(a,f15.5)')' etkf_inno_max = ',etkf_inno_max
   write(stdout,'(a,f15.5)')' etkf_inno_min = ',etkf_inno_min
   write(stdout,'(a,l10)')'   etkf_erro_flg = ',etkf_erro_flg
   write(stdout,'(a,l10)')'   etkf_inno_flg = ',etkf_inno_flg
   write(stdout,'(a,l10)')'   etkf_wrfda = ',etkf_wrfda

   num_members_inv = 1.0 / real(num_members)
   num_members_m1_inv = 1.0 / real(num_members)

   allocate( ens_mean(1:nv) )
   allocate( ens_stdv_pert_prior(1:nv) )
   allocate( ens_stdv_pert_poste(1:nv) )


   write(stdout,'(/a)')' [2] Read observation information.'



   do member = 1, num_members
      write(unit=ce,FMT='(i3.3)') member
      input_file = 'ob.etkf.e'//ce  
      open(unit, file = input_file, status='old')
      read(unit,*) apm_num_obs
      print *, ' Number of unfiltered observations = ',apm_num_obs,member
      if ( member == 1 ) then
         allocate( apm_yo(1:apm_num_obs,1:num_members) )
         allocate( apm_y(1:apm_num_obs,1:num_members) )
         allocate( apm_sigma_o2(1:apm_num_obs,1:num_members) )
         allocate( apm_obs_type(1:apm_num_obs,1:num_members) )
         allocate( apm_type(1:apm_num_obs,1:num_members) )
         allocate( apm_subtype(1:apm_num_obs,1:num_members) )
         allocate( apm_obslat(1:apm_num_obs,1:num_members) )
         allocate( apm_obslon(1:apm_num_obs,1:num_members) )
         allocate( apm_obsprs(1:apm_num_obs,1:num_members) )
         allocate( apm_obstim(1:apm_num_obs,1:num_members) )
         allocate( apm_irec(1:apm_num_obs,1:num_members) )
      end if


      if(etkf_wrfda) then
         do o = 1, apm_num_obs
            read(unit,'(3f17.7)') apm_yo(o,member),apm_y(o,member),apm_sigma_o2(o,member)
         end do
         close(unit)
      else
         do o = 1, apm_num_obs
            read(unit,'(3f17.7,2x,a10,2x,2(f6.0,2x),4(f8.2,2x),i10)') &
            apm_yo(o,member),apm_y(o,member),apm_sigma_o2(o,member),apm_obs_type(o,member), &
            apm_type(o,member),apm_subtype(o,member),apm_obslat(o,member), &
            apm_obslon(o,member),apm_obsprs(o,member),apm_obstim(o,member),apm_irec(o,member)
         end do
      end if
   end do    

   apm_icnt=0



   if(etkf_wrfda) then
      write(stdout,'(a)') 'etkf_wrfda is true: NO ETKF OBS FILTERING '
      apm_icnt = apm_num_obs
   else
      do o = 1, apm_num_obs
         if((apm_type(o,1).ne.220.) .and. (apm_type(o,1).ne.120.)) then
            apm_reject=1
            do member=1,num_members
               apm_irec(o,member)=-999        
            end do 
         else
            apm_reject=0
            do member=1, num_members


               if((etkf_erro_flg .and. &
                  (apm_sigma_o2(o,member)/apm_yo(o,member).le.etkf_erro_min) .or. &
                  (apm_sigma_o2(o,member)/apm_yo(o,member).ge.etkf_erro_max)) .or. &


                  (etkf_inno_flg .and. & 
                  (abs(apm_y(o,member))/apm_yo(o,member).le.etkf_inno_min) .or. &
                  (abs(apm_y(o,member))/apm_yo(o,member).ge.etkf_inno_max))) then  
                  apm_reject=1
               endif            
            end do
            if(apm_reject.eq.1) then 
               do member=1,num_members
                  apm_irec(o,member)=-999        
               end do 
            else         
               apm_icnt=apm_icnt+1
            end if
         end if
      end do
   end if     


   write(stdout,'(a,i10)')'   Number of filtered observations = ', apm_icnt
   nobs_flt = apm_icnt


   if(etkf_wrfda) then
      num_obs=nobs_flt
      allocate( y(1:num_obs,1:num_members) )
      allocate( sigma_o2(1:num_obs) )
      allocate( yo(1:num_obs) )
      apm_icnt=0
      do o=1,num_obs
         yo(o)=apm_yo(o,1)
         sigma_o2(o)=apm_sigma_o2(o,1)
         do member=1,num_members
            y(o,member)=apm_y(o,member)
         end do
      end do 
   else
      if(rand_filt) then
         nseed=rnd_seed
         num_obs=rnd_nobs
         allocate( y(1:num_obs,1:num_members) )
         allocate( sigma_o2(1:num_obs) )
         allocate( yo(1:num_obs) )
         allocate( yo_obs_typ(1:num_obs) )
         allocate( yo_typ(1:num_obs) )
         allocate( yo_subtyp(1:num_obs) )
         allocate( yo_lat(1:num_obs) )
         allocate( yo_lon(1:num_obs) )
         allocate( yo_prs(1:num_obs) )
         allocate( yo_tim(1:num_obs) )
         call rand_filter(apm_y,apm_sigma_o2,apm_yo,apm_obs_type,apm_type,apm_subtype, &
                         apm_obslat,apm_obslon,apm_obsprs,apm_obstim,y,sigma_o2,yo,yo_obs_typ, &
                         yo_typ,yo_subtyp,yo_lat,yo_lon,yo_prs,yo_tim,apm_irec,apm_num_obs, &
                         num_obs,nobs_flt,num_members,nseed)
      else   
         num_obs=nobs_flt
         allocate( y(1:num_obs,1:num_members) )
         allocate( sigma_o2(1:num_obs) )
         allocate( yo(1:num_obs) )
         allocate( yo_obs_typ(1:num_obs) )
         allocate( yo_typ(1:num_obs) )
         allocate( yo_subtyp(1:num_obs) )
         allocate( yo_lat(1:num_obs) )
         allocate( yo_lon(1:num_obs) )
         allocate( yo_prs(1:num_obs) )
         allocate( yo_tim(1:num_obs) )
         apm_icnt=0
         do o=1,apm_num_obs
            if(apm_irec(o,1).ne.-999) then
               apm_icnt=apm_icnt+1
               yo(apm_icnt)=apm_yo(o,1)
               yo_obs_typ(apm_icnt)=apm_obs_type(o,1)
               yo_typ(apm_icnt)=apm_type(o,1)
               yo_subtyp(apm_icnt)=apm_subtype(o,1)
               yo_lat(apm_icnt)=apm_obslat(o,1)
               yo_lon(apm_icnt)=apm_obslon(o,1)
               yo_prs(apm_icnt)=apm_obsprs(o,1)
               yo_tim(apm_icnt)=apm_obstim(o,1)
               sigma_o2(apm_icnt)=apm_sigma_o2(o,1)
               do member=1,num_members
                  y(apm_icnt,member)=apm_y(o,member)
               end do



            endif
         enddo
      endif


      print *, "begin setup of ob.etkf.filt "
      allocate(filt_ob_etkf(num_members))
      do member=1,num_members
         write(unit=ce,FMT='(i3.3)') member
         filt_ob_etkf(member) = 'ob.etkf.filt.e'//ce 
         apm_unit=200+member 
         open(apm_unit, file = filt_ob_etkf(member), status='unknown')
         write(apm_unit,'(i17)') num_obs
      end do
      print *, "complete setup of ob.etkf.filt "


      do member=1,num_members
         apm_icnt=0
         do o=1, num_obs
               apm_icnt=apm_icnt+1
               apm_unit=200+member
               write(apm_unit,'(3f17.7,2x,a10,2x,2(f6.0,2x),4(f8.2,2x),i10)') &
               yo(o),y(o,member),sigma_o2(o),yo_obs_typ(o),yo_typ(o),yo_subtyp(o), &
               yo_lat(o),yo_lon(o),yo_prs(o),yo_tim(o),apm_icnt
         end do
      end do
      do member=1,num_members
         apm_unit=200+member
         close(apm_unit)
      enddo
      deallocate(filt_ob_etkf)
      deallocate(apm_yo,apm_y,apm_sigma_o2,apm_obs_type,apm_type,apm_subtype, &
      apm_obslat,apm_obslon,apm_obsprs,apm_obstim,apm_irec)  
   end if      


   do o = 1, num_obs


      sigma_o2(o) = sigma_o2(o) * sigma_o2(o)


      do member = 1, num_members
         y(o,member) = yo(o) - y(o,member)
      end do
   end do


   write(stdout,'(/a)')' [3] Set up arrays using input ensemble mean forecast'



      input_file = 'etkf_input'
      length = len_trim(input_file)
      rcode = nf_open( input_file(1:length), NF_NOWRITE, cdfid )


      do v = 1, nv


          rcode = nf_inq_varid ( cdfid, cv(v), id_var(v) )
         if ( rcode /= 0 ) then
            write(UNIT=message(1),FMT='(A,A)') &
            cv(v), ' variable is not in input file'
            call da_error("gen_be_etkf.b",429,message(1:1)) 
         end if 


         dimids = 0
         rcode = nf_inq_var( cdfid, id_var(v), cv(v), ivtype, ndims(v), dimids, natts )
         if ( ivtype /= 5 ) then
            write(UNIT=message(1),FMT='(A,A)') cv(v), ' variable is not real type'
            call da_error("gen_be_etkf.b",437,message(1:1))
         end if


         dims(v,:) = 0
         do i = 1, ndims(v)
            rcode = nf_inq_dimlen( cdfid, dimids(i), dims(v,i) )
         end do
      end do


      one = 1
      istart(1) = 1
      do v = 2, nv
         istart(v) = istart(v-1) + product(dims(v-1,1:ndims(v-1)-1))
      end do


      nijkv = istart(nv) + product(dims(nv,1:ndims(nv)-1)) - 1
      allocate( xf_mean(1:nijkv) )


      do v = 1, nv
         index = istart(v)
         if(ndims(v) == 3) then
            allocate( data_r_3d(dims(v,1),dims(v,2),dims(v,3)))
            rcode = nf_get_vara_real( cdfid, id_var(v), one, dims(v,:), data_r_3d)
            do k = 1, dims(v,3)
               do j = 1, dims(v,2)
                  do i = 1, dims(v,1)
                     xf_mean(index) = data_r_3d(i,j,k)
                     index = index + 1
                  end do
               end do
            end do
            deallocate( data_r_3d )
         endif
         if(ndims(v) == 4) then
            allocate( data_r_4d(dims(v,1),dims(v,2),dims(v,3),dims(v,4)))
            rcode = nf_get_vara_real( cdfid, id_var(v), one, dims(v,:), data_r_4d)
            do k = 1, dims(v,3)
               do j = 1, dims(v,2)
                  do i = 1, dims(v,1)
                     xf_mean(index) = data_r_4d(i,j,k,1)
                     index = index + 1
                  end do
               end do
            end do
            deallocate( data_r_4d )
         endif
      enddo


      rcode=nf_inq_varid(cdfid,'XLONG',id_lon_T)
      rcode=nf_inq_varid(cdfid,'XLONG_U',id_lon_u)
      rcode=nf_inq_varid(cdfid,'XLONG_V',id_lon_v)
      rcode=nf_inq_varid(cdfid,'XLAT',id_lat_T)
      rcode=nf_inq_varid(cdfid,'XLAT_U',id_lat_u)
      rcode=nf_inq_varid(cdfid,'XLAT_V',id_lat_v)
      rcode=nf_inq_varid(cdfid,'P',id_prs_T)
      rcode=nf_inq_varid(cdfid,'PB',id_pbs_T)










      dimids(:)=0
      var_name='XLONG'
      rcode=nf_inq_var(cdfid,id_lon_T,var_name,ityp,ndim_lon_T,dimids,natts)
      allocate(dim_lon_T(ndim_lon_T))
      do idim=1,ndim_lon_T
         rcode=nf_inq_dimlen(cdfid,dimids(idim),dim_lon_T(idim))
      enddo

      dimids(:)=0
      var_name='XLONG_U'
      rcode=nf_inq_var(cdfid,id_lon_u,var_name,ityp,ndim_lon_u,dimids,natts)
      allocate(dim_lon_u(ndim_lon_u))
      do idim=1,ndim_lon_u
         rcode=nf_inq_dimlen(cdfid,dimids(idim),dim_lon_u(idim))
      enddo

      dimids(:)=0
      var_name='XLONG_V'
      rcode=nf_inq_var(cdfid,id_lon_v,var_name,ityp,ndim_lon_v,dimids,natts)
      allocate(dim_lon_v(ndim_lon_v))
      do idim=1,ndim_lon_v
         rcode=nf_inq_dimlen(cdfid,dimids(idim),dim_lon_v(idim))
      enddo

      dimids(:)=0
      var_name='XLAT'
      rcode=nf_inq_var(cdfid,id_lat_T,var_name,ityp,ndim_lat_T,dimids,natts)
      allocate(dim_lat_T(ndim_lat_T))
      do idim=1,ndim_lat_T
         rcode=nf_inq_dimlen(cdfid,dimids(idim),dim_lat_T(idim))
      enddo

      dimids(:)=0
      var_name='XLAT_U'
      rcode=nf_inq_var(cdfid,id_lat_u,var_name,ityp,ndim_lat_u,dimids,natts)
      allocate(dim_lat_u(ndim_lat_u))
      do idim=1,ndim_lat_u
         rcode=nf_inq_dimlen(cdfid,dimids(idim),dim_lat_u(idim))
      enddo

      dimids(:)=0
      var_name='XLAT_V'
      rcode=nf_inq_var(cdfid,id_lat_v,var_name,ityp,ndim_lat_v,dimids,natts)
      allocate(dim_lat_v(ndim_lat_v))
      do idim=1,ndim_lat_v
         rcode=nf_inq_dimlen(cdfid,dimids(idim),dim_lat_v(idim))
      enddo

      dimids(:)=0
      var_name='P'
      rcode=nf_inq_var(cdfid,id_prs_T,var_name,ityp,ndim_prs_T,dimids,natts)
      allocate(dim_prs_T(ndim_prs_T))
      do idim=1,ndim_prs_T
         rcode=nf_inq_dimlen(cdfid,dimids(idim),dim_prs_T(idim))
      enddo

      dimids(:)=0
      var_name='PB'
      rcode=nf_inq_var(cdfid,id_pbs_T,var_name,ityp,ndim_pbs_T,dimids,natts)
      allocate(dim_pbs_T(ndim_pbs_T))
      do idim=1,ndim_pbs_T
         rcode=nf_inq_dimlen(cdfid,dimids(idim),dim_pbs_T(idim))
      enddo



      allocate(ens_lon_T(dim_lon_T(1),dim_lon_T(2)))
      allocate(ens_lon_u(dim_lon_u(1),dim_lon_u(2)))
      allocate(ens_lon_v(dim_lon_v(1),dim_lon_v(2)))
      allocate(ens_lat_T(dim_lat_T(1),dim_lat_T(2)))
      allocate(ens_lat_u(dim_lat_u(1),dim_lat_u(2)))
      allocate(ens_lat_v(dim_lat_v(1),dim_lat_v(2)))
      allocate(ens_prs_T(dim_prs_T(1),dim_prs_T(2),dim_prs_T(3)))
      allocate(ens_pbs_T(dim_pbs_T(1),dim_pbs_T(2),dim_pbs_T(3)))


      one(:)=1
      rcode=nf_get_vara_real(cdfid,id_lon_T,one,dim_lon_T,ens_lon_T)

      rcode=nf_get_vara_real(cdfid,id_lon_u,one,dim_lon_u,ens_lon_u)

      rcode=nf_get_vara_real(cdfid,id_lon_v,one,dim_lon_v,ens_lon_v)

      rcode=nf_get_vara_real(cdfid,id_lat_T,one,dim_lat_T,ens_lat_T)

      rcode=nf_get_vara_real(cdfid,id_lat_u,one,dim_lat_u,ens_lat_u)

      rcode=nf_get_vara_real(cdfid,id_lat_v,one,dim_lat_v,ens_lat_v)

      rcode=nf_get_vara_real(cdfid,id_prs_T,one,dim_prs_T,ens_prs_T)

      rcode=nf_get_vara_real(cdfid,id_pbs_T,one,dim_pbs_T,ens_pbs_T)

      rcode = nf_close( cdfid )


   write(stdout,'(/a)')' [4] Extract necessary fields from WRF ensemble forecasts.'



      allocate( xf(1:nijkv,1:num_members) )
      allocate( xf_lon(1:nijkv),xf_lat(1:nijkv))
      allocate( xf_vari(1:nijkv) )


      do member = 1, num_members

         write(UNIT=ce,FMT='(i3.3)')member
         input_file = 'etkf_input.e'//ce
         length = len_trim(input_file)
         rcode = nf_open( input_file(1:length), NF_NOWRITE, cdfid )
         do v = 1, nv
            index = istart(v)
            if(ndims(v) == 3) then
               allocate( data_r_3d(dims(v,1),dims(v,2),dims(v,3)))
               rcode = nf_get_vara_real( cdfid, id_var(v), one, dims(v,:), data_r_3d)
               do k = 1, dims(v,3)
                  do j = 1, dims(v,2)
                     do i = 1, dims(v,1)
                        xf(index,member) = data_r_3d(i,j,k)
                        index = index + 1
                     end do
                  end do
               end do
               deallocate( data_r_3d )
            endif
            if(ndims(v) == 4) then
               allocate( data_r_4d(dims(v,1),dims(v,2),dims(v,3),dims(v,4)))
               rcode = nf_get_vara_real( cdfid, id_var(v), one, dims(v,:), data_r_4d)
               do k = 1, dims(v,3)
                  do j = 1, dims(v,2)
                     do i = 1, dims(v,1)
                        xf(index,member) = data_r_4d(i,j,k,1)
                        index = index + 1
                     end do
                  end do
               end do
               deallocate( data_r_4d )
            endif
         end do
         rcode = nf_close( cdfid )
      end do





      nxs=dims(1,1)
      nys=dims(2,2)
      nzs=dims(3,3)

      napl1=(nxs-1)*(nys-1)
      napl2=2*nv*nzs

      allocate (ijk_idx(1:nv,1:nxs,1:nys,1:nzs))
      allocate (apl_idx(1:napl1,1:napl2))
      allocate (apl_ndim(1:napl1))
      
      ijk_idx=0
      apl_idx=0
      apl_ndim=0

      do v = 1, nv
         index = istart(v)
         if(cv(v) .eq. 'U') then
            do k = 1, dims(v,3)
               do j = 1, dims(v,2)
                  do i = 1, dims(v,1)
                     xf_lat(index) = ens_lat_u(i,j)
                     xf_lon(index) = ens_lon_u(i,j)
                     ijk_idx(v,i,j,k)=index
                     index = index + 1
                  end do
               end do
            end do
         end if
         if(cv(v) .eq. 'V') then
            do k = 1, dims(v,3)
               do j = 1, dims(v,2)
                  do i = 1, dims(v,1)
                     xf_lat(index) = ens_lat_v(i,j)
                     xf_lon(index) = ens_lon_v(i,j)
                     ijk_idx(v,i,j,k)=index
                     index = index + 1
                  end do
               end do
            end do
         end if
         if(cv(v).eq.'W' .or. cv(v).eq.'T' .or. cv(v).eq.'W' .or. &
         cv(v).eq.'PH' .or. cv(v).eq.'QVAPOR') then
            do k = 1, dims(v,3)
               do j = 1, dims(v,2)
                  do i = 1, dims(v,1)
                     xf_lat(index) = ens_lat_T(i,j)
                     xf_lon(index) = ens_lon_T(i,j)
                     ijk_idx(v,i,j,k)=index
                     index = index + 1
                  end do
               end do
            end do
         end if
         if(cv(v) .eq. 'MU') then
            do j = 1, dims(v,2)
               do i = 1, dims(v,1)
                  xf_lat(index) = ens_lat_v(i,j)
                  xf_lon(index) = ens_lon_v(i,j)
                  ijk_idx(v,i,j,1)=index
                  index = index + 1
               end do
            end do
         end if
      end do


      ixy=0
      do j=1,nys-1
         do i=1,nxs-1
            ixy=ixy+1

            idx=0
            do k=1,dims(1,3)
               idx=idx+1
               apl_idx(ixy,idx)=ijk_idx(1,i,j,k)
               if(i.eq.nxs-1) then
                 idx=idx+1
                 apl_idx(ixy,idx)=ijk_idx(1,nxs,j,k)
               endif
            enddo

            do k=1,dims(2,3)
               idx=idx+1
               apl_idx(ixy,idx)=ijk_idx(2,i,j,k)
               if(j.eq.nys-1) then
                 idx=idx+1
                 apl_idx(ixy,idx)=ijk_idx(2,i,nys,k)
               endif
            enddo

            do k=1,dims(3,3)
               idx=idx+1
               apl_idx(ixy,idx)=ijk_idx(3,i,j,k)
            enddo

            do k=1,dims(4,3)
               idx=idx+1
               apl_idx(ixy,idx)=ijk_idx(4,i,j,k)
            enddo

            do k=1,dims(5,3)
               idx=idx+1
               apl_idx(ixy,idx)=ijk_idx(5,i,j,k)
            enddo

            do k=1,dims(6,3)
               idx=idx+1
               apl_idx(ixy,idx)=ijk_idx(6,i,j,k)
            enddo

            idx=idx+1
            apl_idx(ixy,idx)=ijk_idx(7,i,j,1)
            apl_ndim(ixy)=idx
         enddo  
      enddo


      do ijkv = 1, nijkv
         xf(ijkv,1:num_members) = xf(ijkv,1:num_members) - xf_mean(ijkv)
         xf_vari(ijkv) = sum(xf(ijkv,1:num_members)**2) * num_members_m1_inv
      end do



     do v = 1, nv
        iend = istart(v) + product(dims(v,1:ndims(v)-1)) - 1
        ens_mean(v) = sum(xf_mean(istart(v):iend)) / real(iend - istart(v) + 1)
        ens_stdv_pert_prior(v) =sqrt( sum(xf_vari(istart(v):iend)) / &
                                    real(iend - istart(v) + 1) )
     end do


   write(stdout,'(/a)')' [5] Call ETKF:'



     if (letkf_flg) then
        call da_solve_letkf( nijkv, num_members, num_obs, xf, y, sigma_o2, yo, nout, &
                       naccumt1, naccumt2, nstartaccum1, nstartaccum2, tainflatinput, &
                       rhoinput, infl_fac_file, eigen_val_file, inno2_val_file, &
                       proj2_val_file, infl_fac_TRNK, infl_fac_WG03, infl_fac_WG07, &
                       infl_fac_BOWL, yo_lat, yo_lon, yo_prs, xf_lat, xf_lon, ijk_idx, &
                       apl_idx, apl_ndim, napl1, napl2, nxs, nys, nzs, nv, infl_let_file)
     else
        call da_solve_etkf( nijkv, num_members, num_obs, xf, y, sigma_o2, yo, nout, &
                       naccumt1, naccumt2, nstartaccum1, nstartaccum2, tainflatinput, &
                       rhoinput, infl_fac_file, eigen_val_file, inno2_val_file, &
                       proj2_val_file, infl_fac_TRNK, infl_fac_WG03, infl_fac_WG07, &
                       infl_fac_BOWL)
     endif


     do ijkv = 1, nijkv
        xf_vari(ijkv) = sum(xf(ijkv,1:num_members)**2) * num_members_m1_inv
     end do
     write(stdout,'(5a)')'   v', ' Variable  ', '    Ensemble Mean', &
                       '  Prior Pert StDv', ' Post. Pert. StDv'
     do v = 1, nv
        iend = istart(v) + product(dims(v,1:ndims(v)-1)) - 1
        ens_stdv_pert_poste(v) =sqrt( sum(xf_vari(istart(v):iend)) / &
                                    real(iend - istart(v) + 1) )
        write(stdout,'(i4,1x,a10,3f17.7)')v, cv(v), &
        ens_mean(v), ens_stdv_pert_prior(v), ens_stdv_pert_poste(v)
     end do


   write(stdout,'(/a)')' [6] Output ETKF analysis ensemble:'

     do member = 1, num_members
        write(UNIT=ce,FMT='(i3.3)')member
        input_file = 'etkf_output.e'//ce
        rcode = nf_open( trim(input_file), NF_WRITE, cdfid )
        if ( rcode /= 0 ) then
           print *, 'MISSING FILE ',trim(input_file)
           print *, 'ERROR OPENING OUTPUT FILE: ', nf_strerror(rcode)
           stop
        end if
        do v = 1, nv
           index = istart(v)
           if(ndims(v) == 3) then
              allocate( data_r_3d(dims(v,1),dims(v,2),dims(v,3)))
              rcode = nf_get_vara_real( cdfid, id_var(v), one, dims(v,:), data_r_3d)
              do k = 1, dims(v,3)
                 do j = 1, dims(v,2)
                    do i = 1, dims(v,1)



                       data_r_3d(i,j,k) = xf(index,member)
                       index = index + 1
                    end do
                 end do
              end do
              call ncvpt( cdfid, id_var(v), one, dims(v,:), data_r_3d, rcode)
              deallocate( data_r_3d )
           endif
           if(ndims(v) == 4) then
              allocate( data_r_4d(dims(v,1),dims(v,2),dims(v,3),dims(v,4)))
              rcode = nf_get_vara_real( cdfid, id_var(v), one, dims(v,:), data_r_4d)
              do k = 1, dims(v,3)
                 do j = 1, dims(v,2)
                    do i = 1, dims(v,1)



                       data_r_4d(i,j,k,1) = xf(index,member)
                       index = index + 1
                    end do
                 end do
              end do
              call ncvpt( cdfid, id_var(v), one, dims(v,:), data_r_4d, rcode)
              deallocate( data_r_4d )
           endif
        enddo
        rcode = nf_close( cdfid )
     end do
     deallocate( ens_mean )
     deallocate( ens_stdv_pert_prior )
     deallocate( ens_stdv_pert_poste )

end program gen_be_etkf
