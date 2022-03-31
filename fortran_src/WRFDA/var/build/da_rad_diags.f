












program da_rad_diags

















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



   namelist /record1/ nproc, instid, file_prefix, start_date, end_date, cycle_period
           
           
           
           
           
   integer, parameter                     :: maxnum = 20, maxlvl = 100
   integer                                :: nml_unit = 87
   integer                                :: nproc, nlev, ilev, ich
   integer                                :: cycle_period, nlev_rtm, nlev_mdl
   character(len=20), dimension(maxnum)   :: instid
   character(len=6)                       :: file_prefix
   character(len=10)                      :: start_date, end_date



   character(len=200)                     :: ncname
   integer                                :: ncid, dimid, varid
   integer, dimension(3)                  :: ishape, istart, icount

   logical                                :: amsr2
   logical                                :: isfile, prf_found, jac_found
   integer, parameter                     :: datelen1 = 10
   integer, parameter                     :: datelen2 = 19
   real*4, parameter                      :: missing_r = -888888.00
   character(len=20)                      :: rtm_option    
   character(len=250)                     :: buf, inst
   character(len=7)                       :: numbuf
   character(len=datelen1)                :: valid_date
   integer                                :: ninst, iinst, iproc, ipixel, ifirst
   integer                                :: ios, i, n, ips, ipe, nerr, itime, itmp
   integer                                :: ntime, nchan, total_npixel
   integer, dimension(:), allocatable     :: ichan, npixel, iunit, scanpos, isflg
   integer, dimension(:), allocatable     :: landsea_mask, soiltyp, vegtyp
   real*4,  dimension(:), allocatable     :: lat, lon, elv, elev
   real*4,  dimension(:), allocatable     :: ret_clw
   real*4,  dimension(:), allocatable     :: satzen, satazi, t2m, mr2m, u10, v10, ps, ts
   real*4,  dimension(:), allocatable     :: smois, tslb, snowh, vegfra, clwp
   integer, dimension(:,:), allocatable   :: tb_qc
   real*4,  dimension(:,:), allocatable   :: tb_obs, tb_bak, tb_inv, tb_oma, tb_err, ems, ems_jac
   real*4,  dimension(:,:), allocatable   :: prf_pfull, prf_phalf, prf_t, prf_q, prf_water
   real*4,  dimension(:,:), allocatable   :: prf_ice, prf_rain, prf_snow, prf_grau, prf_hail
   real*4,  dimension(:,:), allocatable   :: prf_water_reff, prf_ice_reff, prf_rain_reff
   real*4,  dimension(:,:), allocatable   :: prf_snow_reff, prf_grau_reff, prf_hail_reff
   real*4,  dimension(:,:), allocatable   :: rtm_prf_p, rtm_prf_t, rtm_prf_q
   real*4,  dimension(:,:), allocatable   :: mdl_prf_p, mdl_prf_t, mdl_prf_q, mdl_prf_qcw, mdl_prf_qrn
   real*4,  dimension(:,:,:), allocatable :: prf_t_jac, prf_q_jac, prf_der_trans, prf_trans_jac, prf_trans, prf_lod_jac, prf_lod
   real*4,  dimension(:,:,:), allocatable :: prf_water_jac, prf_ice_jac, prf_rain_jac, prf_snow_jac, prf_grau_jac, prf_hail_jac
   real*4,  dimension(:,:,:), allocatable :: prf_water_reff_jac, prf_ice_reff_jac, prf_rain_reff_jac
   real*4,  dimension(:,:,:), allocatable :: prf_snow_reff_jac, prf_grau_reff_jac, prf_hail_reff_jac
   character(len=200), dimension(:), allocatable      :: inpname
   character(len=datelen1), dimension(:), allocatable :: datestr1          
   character(len=datelen2), dimension(:), allocatable :: datestr2          



   instid   = ''
   ninst    = 0
   file_prefix = 'inv'
   ilev = 0
   nlev = 0
   nlev_rtm = 0
   nlev_mdl = 0
   prf_found = .false.
   jac_found = .false.
   rtm_option = 'unknown'



   open(unit=nml_unit, file='namelist.da_rad_diags', status='old', form='formatted')
   read(unit=nml_unit, nml=record1, iostat=ios)
   write(0,nml=record1)
   if ( ios /= 0 ) then
      write(0,*) 'Error reading namelist record1'
      stop
   end if



   do i = 1, maxnum
      if ( len_trim(instid(i)) /= 0 ) then
         ninst = ninst + 1
      end if
   end do



   ntime = 0
   valid_date = start_date
   do while ( valid_date <= end_date )
      ntime = ntime + 1
      call advance_cymdh(valid_date, cycle_period, valid_date)
   end do
   write(0,*) 'ntime = ', ntime

   allocate ( datestr1(ntime) )
   valid_date = start_date
   datestr1(1) = start_date
   do i = 2, ntime
      call advance_cymdh(datestr1(i-1), cycle_period, datestr1(i))
   end do

ntime_loop: do itime = 1, ntime

   write(0,*) '=================='
   write(0,*) trim(datestr1(itime))

   ninst_loop: do iinst = 1, ninst

      write(0,*) '------------------'
      write(0,*) trim(instid(iinst))

      amsr2 = index(instid(iinst),'amsr2') > 0

      nerr = 0
      total_npixel = 0
      ips = 0
      ipe = 0

      allocate ( npixel(0:nproc-1) )
      allocate ( iunit(0:nproc-1) )
      allocate ( inpname(0:nproc-1) )

      nproc_loop_1: do iproc = 0, nproc - 1   

         write(unit=inpname(iproc), fmt='(a,i4.4)')  &
            trim(datestr1(itime))//'/'//trim(adjustl(file_prefix))//'_'//trim(instid(iinst))//'.', iproc
         iunit(iproc) = 101 + iproc
         inquire(file=trim(inpname(iproc)), exist=isfile)
         if ( .not. isfile ) Then
            write(0,*) 'Error opening innovation radiance file ', trim(inpname(iproc))
            nerr = nerr + 1
            if ( nerr == nproc ) then
               write(0,*) 'found no vaild files for ', trim(instid(iinst))
               deallocate ( npixel )
               deallocate ( iunit )
               deallocate ( inpname )
               cycle ninst_loop
            end if
            cycle nproc_loop_1
         end if
 
         open(unit=iunit(iproc),file=trim(inpname(iproc)),form='formatted',iostat=ios)
         read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf   
 
         inst = buf(1:(index(buf,'number-of-pixels')-2))   
         
         
         
         numbuf = buf((index(buf,'channel-number')-8):(index(buf,'channel-number')-2))
         read(numbuf,'(i7)') npixel(iproc)

         total_npixel = total_npixel + npixel(iproc)

         itmp = 0
         do while ( ios == 0 .and. itmp < 2 )
            read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf
            if ( index(buf,'INFO :') > 0 ) then
               itmp = itmp + 1
            else
               if ( index(buf,'EMS_JACOBIAN') > 0 ) then
                  jac_found = .true.
               end if
               if ( index(buf,'level fullp(mb)') > 0 ) then
                  prf_found = .true.
                  rtm_option = 'CRTM'
               end if
               if ( index(buf,'RTM_level pres(mb)') > 0 ) then
                  prf_found = .true.
                  rtm_option = 'RTTOV'
               end if
            end if
            if ( rtm_option /= 'unknown' ) exit
         end do

      end do nproc_loop_1

      write(0,*) 'total_npixel = ', total_npixel

      ifirst = 1
      nproc_loop_2: do iproc = 0, nproc - 1

         inquire(file=trim(inpname(iproc)), exist=isfile)
         if ( .not. isfile ) cycle nproc_loop_2
         rewind(iunit(iproc))
         read(unit=iunit(iproc),fmt='(a)') buf
         
         
         
         numbuf = buf((index(buf,'index-of-channels')-6):(index(buf,'index-of-channels')-2))
         read(numbuf,'(i5)') nchan

         if ( .not. allocated(ichan) ) allocate (   ichan(1:nchan) )

         read(unit=iunit(iproc),fmt='(10i5)',iostat=ios) ichan
         read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf      
         read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf      

         if ( ifirst == 1 ) then

            allocate (     datestr2(1:total_npixel) )
            allocate (      scanpos(1:total_npixel) )
            allocate ( landsea_mask(1:total_npixel) )
            allocate (          elv(1:total_npixel) )
            allocate (          lat(1:total_npixel) )
            allocate (          lon(1:total_npixel) )
            allocate (       satzen(1:total_npixel) )
            allocate (       satazi(1:total_npixel) )
            allocate (      ret_clw(1:total_npixel) ) 
            allocate (          t2m(1:total_npixel) )
            allocate (         mr2m(1:total_npixel) )
            allocate (          u10(1:total_npixel) )
            allocate (          v10(1:total_npixel) )
            allocate (           ps(1:total_npixel) )
            allocate (           ts(1:total_npixel) )
            allocate (        smois(1:total_npixel) )
            allocate (         tslb(1:total_npixel) )
            allocate (        snowh(1:total_npixel) )
            allocate (        isflg(1:total_npixel) )
            allocate (      soiltyp(1:total_npixel) )
            allocate (       vegtyp(1:total_npixel) )
            allocate (       vegfra(1:total_npixel) )
            allocate (         elev(1:total_npixel) )
            allocate (         clwp(1:total_npixel) ) 
            allocate ( tb_obs(1:nchan,1:total_npixel) )
            allocate ( tb_bak(1:nchan,1:total_npixel) )
            allocate ( tb_inv(1:nchan,1:total_npixel) )
            allocate ( tb_oma(1:nchan,1:total_npixel) )
            allocate ( tb_err(1:nchan,1:total_npixel) )
            allocate (  tb_qc(1:nchan,1:total_npixel) )
            allocate (    ems(1:nchan,1:total_npixel) )
            if ( jac_found ) then
               allocate ( ems_jac(1:nchan,1:total_npixel) )
            end if
            if ( prf_found .and. (rtm_option == 'CRTM') ) then
               allocate ( prf_pfull(1:maxlvl,1:total_npixel) )
               allocate ( prf_phalf(1:maxlvl,1:total_npixel) )
               allocate ( prf_t(1:maxlvl,1:total_npixel) )
               allocate ( prf_q(1:maxlvl,1:total_npixel) )
               allocate ( prf_water(1:maxlvl,1:total_npixel) )
               allocate ( prf_ice(1:maxlvl,1:total_npixel) )
               allocate ( prf_rain(1:maxlvl,1:total_npixel) )
               allocate ( prf_snow(1:maxlvl,1:total_npixel) )
               allocate ( prf_grau(1:maxlvl,1:total_npixel) )
               allocate ( prf_hail(1:maxlvl,1:total_npixel) )
               allocate ( prf_water_reff(1:maxlvl,1:total_npixel) )
               allocate ( prf_ice_reff(1:maxlvl,1:total_npixel) )
               allocate ( prf_rain_reff(1:maxlvl,1:total_npixel) )
               allocate ( prf_snow_reff(1:maxlvl,1:total_npixel) )
               allocate ( prf_grau_reff(1:maxlvl,1:total_npixel) )
               allocate ( prf_hail_reff(1:maxlvl,1:total_npixel) )
               if ( jac_found ) then
                  allocate ( prf_t_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_q_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_der_trans(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_trans_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_trans(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_lod_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_lod(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_water_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_ice_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_rain_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_snow_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_grau_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_hail_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_water_reff_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_ice_reff_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_rain_reff_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_snow_reff_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_grau_reff_jac(1:maxlvl,1:nchan,1:total_npixel) )
                  allocate ( prf_hail_reff_jac(1:maxlvl,1:nchan,1:total_npixel) )
               end if
            end if
            if ( prf_found .and. (rtm_option == 'RTTOV') ) then
               allocate ( rtm_prf_p(1:maxlvl,1:total_npixel) )
               allocate ( rtm_prf_t(1:maxlvl,1:total_npixel) )
               allocate ( rtm_prf_q(1:maxlvl,1:total_npixel) )  
               allocate ( mdl_prf_p(1:maxlvl,1:total_npixel) )
               allocate ( mdl_prf_t(1:maxlvl,1:total_npixel) )
               allocate ( mdl_prf_q(1:maxlvl,1:total_npixel) )  
               allocate ( mdl_prf_qcw(1:maxlvl,1:total_npixel) )
               allocate ( mdl_prf_qrn(1:maxlvl,1:total_npixel) )
               
               rtm_prf_p = missing_r
               rtm_prf_t = missing_r
               rtm_prf_q = missing_r
               mdl_prf_p = missing_r
               mdl_prf_t = missing_r
               mdl_prf_q = missing_r
               mdl_prf_qcw = missing_r
               mdl_prf_qrn = missing_r
            end if
            
            tb_obs = missing_r
            tb_bak = missing_r
            tb_inv = missing_r
            tb_oma = missing_r
            tb_err = missing_r
            ncname = 'diags_'//trim(instid(iinst))//"_"//datestr1(itime)//'.nc'
            ios = NF_CREATE(trim(ncname), NF_CLOBBER, ncid)  
                                                             
                                                             
            if ( ios /= 0 ) then
               write(0,*) 'Error creating netcdf file: ', ncname
               stop
            end if
             
            ifirst = 0

         end if   



         ips = ipe + 1
         ipe = ipe + npixel(iproc)
         write(0,*) 'Processing pixels ', ips, ' to ', ipe

         npixel_loop: do ipixel = ips, ipe

            read(unit=iunit(iproc),fmt='(7x,i7,2x,a19,i3,i3,f6.0,4f8.2,f8.3)',iostat=ios)  &
               n, datestr2(ipixel), scanpos(ipixel), landsea_mask(ipixel), elv(ipixel), &
               lat(ipixel), lon(ipixel), satzen(ipixel), satazi(ipixel), ret_clw(ipixel)
            if ( scanpos(ipixel) > 360 .or. scanpos(ipixel) == 0 ) then
               backspace(iunit(iproc))
               read(unit=iunit(iproc),fmt='(7x,i7,2x,a19,i6,i3,f6.0,4f8.2,f8.3)',iostat=ios) &
               n, datestr2(ipixel), scanpos(ipixel), landsea_mask(ipixel), elv(ipixel), &
               lat(ipixel), lon(ipixel), satzen(ipixel), satazi(ipixel), ret_clw(ipixel)
            end if
            read(unit=iunit(iproc),fmt='(14x,9f10.2,3i3,3f10.2)',iostat=ios)  &
               t2m(ipixel), mr2m(ipixel), u10(ipixel), v10(ipixel), ps(ipixel), ts(ipixel), &
               smois(ipixel), tslb(ipixel), snowh(ipixel), isflg(ipixel), soiltyp(ipixel),  &
               vegtyp(ipixel), vegfra(ipixel), elev(ipixel), clwp(ipixel)
            read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf           
            read(unit=iunit(iproc),fmt='(10f11.2)',iostat=ios) tb_obs(:,ipixel)
            read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf           
            read(unit=iunit(iproc),fmt='(10f11.2)',iostat=ios) tb_bak(:,ipixel)
            read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf           
            read(unit=iunit(iproc),fmt='(10f11.2)',iostat=ios) tb_inv(:,ipixel)
            read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf           
            if ( buf(1:3) == "OMA" ) then
               read(unit=iunit(iproc),fmt='(10f11.2)',iostat=ios) tb_oma(:,ipixel)
               read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf        
            end if
            read(unit=iunit(iproc),fmt='(10f11.2)',iostat=ios) ems(:,ipixel)
            read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf           
            if ( buf(1:12) == "EMS_JACOBIAN" ) then
               read(unit=iunit(iproc),fmt='(10f10.3)',iostat=ios) ems_jac(:,ipixel)
               read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf        
            end if
            read(unit=iunit(iproc),fmt='(10f11.2)',iostat=ios) tb_err(:,ipixel)
            read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf           
            read(unit=iunit(iproc),fmt='(10i11)',iostat=ios  ) tb_qc(:,ipixel)
            read(unit=iunit(iproc),fmt='(a)',iostat=ios) buf
            if ( buf(1:4) == "INFO" ) then
               backspace(iunit(iproc))
               cycle npixel_loop
            else
               if ( index(buf(1:6),'level') > 0 ) then    
                  ilev = 0
                  do while ( ios == 0 )
                     ilev = ilev + 1
                     read(unit=iunit(iproc),fmt='(3x,2f10.2,f8.2,13f8.3)',iostat=ios)        &
                        prf_pfull(ilev,ipixel), prf_phalf(ilev,ipixel), prf_t(ilev,ipixel),  &
                        prf_q(ilev,ipixel), prf_water(ilev,ipixel), prf_ice(ilev,ipixel),    &
                        prf_rain(ilev,ipixel), prf_snow(ilev,ipixel), prf_grau(ilev,ipixel), &
                        prf_hail(ilev,ipixel),prf_water_reff(ilev,ipixel),                   &
                        prf_ice_reff(ilev,ipixel), prf_rain_reff(ilev,ipixel),               &
                        prf_snow_reff(ilev,ipixel), prf_grau_reff(ilev,ipixel),              &
                        prf_hail_reff(ilev,ipixel)
                  end do
                  nlev = ilev - 1
                  
               else if ( index(buf, 'RTM_level') > 0 ) then  
                  
                  ilev = 0
                  do while ( ios == 0 )
                     ilev = ilev + 1
                     read(unit=iunit(iproc),fmt='(3x,f10.2,f8.2,e11.4)',iostat=ios)        &
                        rtm_prf_p(ilev,ipixel), rtm_prf_t(ilev,ipixel),                    &
                        rtm_prf_q(ilev,ipixel)  
                  end do
                  nlev_rtm = ilev - 1
                  
                  ios  = 0
                  ilev = 0
                  do while ( ios == 0 )
                     ilev = ilev + 1
                     read(unit=iunit(iproc),fmt='(3x,f10.2,f8.2,3e11.4)',iostat=ios)        &
                        mdl_prf_p(ilev,ipixel), mdl_prf_t(ilev,ipixel),                     &
                        mdl_prf_q(ilev,ipixel), mdl_prf_qcw(ilev,ipixel),                   &
                        mdl_prf_qrn(ilev,ipixel)
                  end do
                  nlev_mdl = ilev - 1
               end if
    
               if ( jac_found ) then
                  ios = 0
                  do while ( ios == 0 )
                     
                     read(unit=iunit(iproc),fmt='(i5,i3,10x,19f14.7)',iostat=ios)               &
                        ich, ilev, prf_t_jac(ilev,ich,ipixel), prf_q_jac(ilev,ich,ipixel),     & 
                        prf_der_trans(ilev,ich,ipixel),                                        & 
                        prf_trans_jac(ilev,ich,ipixel), prf_trans(ilev,ich,ipixel),       & 
                        prf_lod_jac(ilev,ich,ipixel), prf_lod(ilev,ich,ipixel),           & 
                        prf_water_jac(ilev,ich,ipixel), prf_ice_jac(ilev,ich,ipixel),          &
                        prf_rain_jac(ilev,ich,ipixel), prf_snow_jac(ilev,ich,ipixel),          &
                        prf_grau_jac(ilev,ich,ipixel), prf_hail_jac(ilev,ich,ipixel),          &
                        prf_water_reff_jac(ilev,ich,ipixel), prf_ice_reff_jac(ilev,ich,ipixel), &
                        prf_rain_reff_jac(ilev,ich,ipixel), prf_snow_reff_jac(ilev,ich,ipixel), &
                        prf_grau_reff_jac(ilev,ich,ipixel), prf_hail_reff_jac(ilev,ich,ipixel)
                     nlev = max(nlev,ilev)
                  end do
                  backspace(iunit(iproc))
               else
                  backspace(iunit(iproc))
               end if
            end if 

         end do npixel_loop

         close(iunit(iproc))

      end do nproc_loop_2

      write(0,*) 'Writing out data in netCDF format...'
      
      
      
      ios = NF_DEF_DIM(ncid, 'nchan', nchan, dimid)
      ios = NF_DEF_DIM(ncid, 'npixel', total_npixel, dimid)
      ios = NF_DEF_DIM(ncid, 'DateStrLen', datelen2, dimid)
      if ( prf_found .or. jac_found ) then
         if ( rtm_option == 'RTTOV' ) then
            nlev = max(nlev_rtm,nlev_mdl)
         end if
         ios = NF_DEF_DIM(ncid, 'nlev', nlev, dimid)
      end if
      
      
      
      if ( trim(rtm_option) == 'unknown' ) then
         if ( maxval(mr2m) < -8000.0 ) then
            rtm_option = 'CRTM'
         else
            rtm_option = 'RTTOV'
         end if
      end if
      ios = NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'rtm_option', 20, rtm_option)
      
      
      
      
      
      ios = NF_INQ_DIMID(ncid, 'DateStrLen', dimid)
      ishape(1) = dimid
      ios = NF_INQ_DIMID(ncid, 'npixel', dimid)
      ishape(2) = dimid
      ios = NF_DEF_VAR(ncid, 'date',   NF_CHAR,  2, ishape(1:2), varid)
      
      
      
      ios = NF_INQ_DIMID(ncid, 'nchan', dimid)
      ishape(1) = dimid
      ios = NF_INQ_DIMID(ncid, 'npixel', dimid)
      ishape(2) = dimid
      ios = NF_DEF_VAR(ncid, 'tb_obs', NF_FLOAT, 2, ishape(1:2), varid)
      ios = NF_PUT_ATT_REAL(ncid, varid, 'missing_value', NF_FLOAT, 1, missing_r)
      ios = NF_DEF_VAR(ncid, 'tb_bak', NF_FLOAT, 2, ishape(1:2), varid)
      ios = NF_PUT_ATT_REAL(ncid, varid, 'missing_value', NF_FLOAT, 1, missing_r)
      ios = NF_DEF_VAR(ncid, 'tb_inv', NF_FLOAT, 2, ishape(1:2), varid)
      ios = NF_PUT_ATT_REAL(ncid, varid, 'missing_value', NF_FLOAT, 1, missing_r)
      ios = NF_DEF_VAR(ncid, 'tb_oma', NF_FLOAT, 2, ishape(1:2), varid)
      ios = NF_PUT_ATT_REAL(ncid, varid, 'missing_value', NF_FLOAT, 1, missing_r)
      ios = NF_DEF_VAR(ncid, 'ems',    NF_FLOAT, 2, ishape(1:2), varid)
      if ( jac_found ) then
         ios = NF_DEF_VAR(ncid, 'ems_jac',NF_FLOAT, 2, ishape(1:2), varid)
      end if
      ios = NF_DEF_VAR(ncid, 'tb_err', NF_FLOAT, 2, ishape(1:2), varid)
      ios = NF_DEF_VAR(ncid, 'tb_qc',  NF_INT,   2, ishape(1:2), varid)
      
      
      
      if ( prf_found ) then
         ios = NF_INQ_DIMID(ncid, 'nlev', dimid)
         ishape(1) = dimid
         ios = NF_INQ_DIMID(ncid, 'npixel', dimid)
         ishape(2) = dimid
         if ( rtm_option == 'CRTM' ) then
            ios = NF_DEF_VAR(ncid, 'prf_pfull', NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_phalf', NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_t',     NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_q',     NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_water', NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_ice',   NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_rain',  NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_snow',  NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_grau',  NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_hail',  NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_water_reff', NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_ice_reff',   NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_rain_reff',  NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_snow_reff',  NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_grau_reff',  NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'prf_hail_reff',  NF_FLOAT, 2, ishape(1:2), varid)
         else if ( rtm_option == 'RTTOV' ) then
            ios = NF_DEF_VAR(ncid, 'rtm_prf_p', NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'rtm_prf_t', NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'rtm_prf_q', NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'mdl_prf_p', NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'mdl_prf_t', NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'mdl_prf_q', NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'mdl_prf_qcw', NF_FLOAT, 2, ishape(1:2), varid)
            ios = NF_DEF_VAR(ncid, 'mdl_prf_qrn', NF_FLOAT, 2, ishape(1:2), varid)
         end if
      end if
      
      
      
      if ( jac_found ) then
         ios = NF_INQ_DIMID(ncid, 'nlev', dimid)
         ishape(1) = dimid
         ios = NF_INQ_DIMID(ncid, 'nchan', dimid)
         ishape(2) = dimid
         ios = NF_INQ_DIMID(ncid, 'npixel', dimid)
         ishape(3) = dimid
         ios = NF_DEF_VAR(ncid, 'prf_t_jac',     NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_q_jac',     NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_der_trans',     NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_trans_jac', NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_trans', NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_lod_jac', NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_lod', NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_water_jac', NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_ice_jac',   NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_rain_jac',  NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_snow_jac',  NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_grau_jac',  NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_hail_jac',  NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_water_reff_jac', NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_ice_reff_jac',   NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_rain_reff_jac',  NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_snow_reff_jac',  NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_grau_reff_jac',  NF_FLOAT, 3, ishape, varid)
         ios = NF_DEF_VAR(ncid, 'prf_hail_reff_jac',  NF_FLOAT, 3, ishape, varid)
      end if
      
      
      
      ios = NF_INQ_DIMID(ncid, 'nchan', dimid)
      ishape(1) = dimid
      ios = NF_DEF_VAR(ncid, 'ichan',  NF_INT,   1, ishape(1), varid)   
      
      
      
      ios = NF_INQ_DIMID(ncid, 'npixel', dimid)
      ishape(1) = dimid
      ios = NF_DEF_VAR(ncid, 'scanpos',      NF_INT,   1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'landsea_mask', NF_INT,   1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'elv',          NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'lat',          NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'lon',          NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'satzen',       NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'satazi',       NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 't2m',          NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'mr2m',         NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'u10',          NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'v10',          NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'ps',           NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'ts',           NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'tslb',         NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'smois',        NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'snowh',        NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'isflg',        NF_INT,   1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'soiltyp',      NF_INT,   1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'vegtyp',       NF_INT,   1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'vegfra',       NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'elev',         NF_FLOAT, 1, ishape(1), varid)
      ios = NF_DEF_VAR(ncid, 'clwp',         NF_FLOAT, 1, ishape(1), varid)
      if ( amsr2 ) then
         ios = NF_DEF_VAR(ncid, 'ret_clw',   NF_FLOAT, 1, ishape(1), varid)
      end if

      ios = NF_ENDDEF(ncid)
      
      
      istart(1) = 1
      istart(2) = 1
      icount(1) = datelen2
      icount(2) = total_npixel
      ios = NF_INQ_VARID (ncid, 'date', varid)
      ios = NF_PUT_VARA_TEXT(ncid, varid, istart, icount, datestr2)
      
      
      
      istart(1) = 1
      istart(2) = 1
      icount(1) = nchan
      icount(2) = total_npixel
      ios = NF_INQ_VARID (ncid, 'tb_obs', varid)
      ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), tb_obs)
      ios = NF_INQ_VARID (ncid, 'tb_bak', varid)
      ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), tb_bak)
      ios = NF_INQ_VARID (ncid, 'tb_inv', varid)
      ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), tb_inv)
      ios = NF_INQ_VARID (ncid, 'tb_oma', varid)
      ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), tb_oma)
      ios = NF_INQ_VARID (ncid, 'ems', varid)
      ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), ems)
      if ( jac_found ) then
         ios = NF_INQ_VARID (ncid, 'ems_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), ems_jac)
      end if
      ios = NF_INQ_VARID (ncid, 'tb_err', varid)
      ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), tb_err)
      ios = NF_INQ_VARID (ncid, 'tb_qc', varid)
      ios = NF_PUT_VARA_INT(ncid,  varid, istart(1:2), icount(1:2), tb_qc)
      
      
      
      if ( prf_found ) then
         istart(1) = 1
         istart(2) = 1
         icount(1) = nlev
         icount(2) = total_npixel
         if ( rtm_option == 'CRTM' ) then
            ios = NF_INQ_VARID (ncid, 'prf_pfull', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_pfull(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_phalf', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_phalf(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_t', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_t(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_q', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_q(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_water', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_water(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_ice', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_ice(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_rain', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_rain(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_snow', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_snow(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_grau', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_grau(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_hail', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_hail(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_water_reff', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_water_reff(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_ice_reff', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_ice_reff(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_rain_reff', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_rain_reff(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_snow_reff', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_snow_reff(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_grau_reff', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_grau_reff(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'prf_hail_reff', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), prf_hail_reff(1:nlev,:))
         else if ( rtm_option == 'RTTOV' ) then
            ios = NF_INQ_VARID (ncid, 'rtm_prf_p', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), rtm_prf_p(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'rtm_prf_t', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), rtm_prf_t(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'rtm_prf_q', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), rtm_prf_q(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'mdl_prf_p', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), mdl_prf_p(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'mdl_prf_t', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), mdl_prf_t(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'mdl_prf_q', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), mdl_prf_q(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'mdl_prf_qcw', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), mdl_prf_qcw(1:nlev,:))
            ios = NF_INQ_VARID (ncid, 'mdl_prf_qrn', varid)
            ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:2), icount(1:2), mdl_prf_qrn(1:nlev,:))
         end if
      end if
      
      
      
      if ( jac_found ) then
         istart(1) = 1
         istart(2) = 1
         istart(3) = 1
         icount(1) = nlev
         icount(2) = nchan
         icount(3) = total_npixel
         ios = NF_INQ_VARID (ncid, 'prf_t_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_t_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_q_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_q_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_der_trans', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_der_trans(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_trans_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_trans_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_trans', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_trans(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_lod_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_lod_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_lod', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_lod(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_water_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_water_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_ice_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_ice_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_rain_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_rain_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_snow_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_snow_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_grau_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_grau_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_hail_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_hail_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_water_reff_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_water_reff_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_ice_reff_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_ice_reff_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_rain_reff_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_rain_reff_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_snow_reff_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_snow_reff_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_grau_reff_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_grau_reff_jac(1:nlev,:,:))
         ios = NF_INQ_VARID (ncid, 'prf_hail_reff_jac', varid)
         ios = NF_PUT_VARA_REAL(ncid, varid, istart(1:3), icount(1:3), prf_hail_reff_jac(1:nlev,:,:))
      end if
      
      
      
      istart(1) = 1
      icount(1) = nchan
      ios = NF_INQ_VARID (ncid, 'ichan', varid)
      ios = NF_PUT_VARA_INT(ncid,  varid, istart(1), icount(1), ichan)
      
      
      
      istart(2) = 1
      icount(2) = total_npixel
      ios = NF_INQ_VARID (ncid, 'scanpos', varid)
      ios = NF_PUT_VARA_INT(ncid,  varid, istart(2), icount(2), scanpos)
      ios = NF_INQ_VARID (ncid, 'landsea_mask', varid)
      ios = NF_PUT_VARA_INT(ncid,  varid, istart(2), icount(2), landsea_mask)
      ios = NF_INQ_VARID (ncid, 'elv', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), elv)
      ios = NF_INQ_VARID (ncid, 'lat', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), lat)
      ios = NF_INQ_VARID (ncid, 'lon', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), lon)
      ios = NF_INQ_VARID (ncid, 'satzen', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), satzen)
      ios = NF_INQ_VARID (ncid, 'satazi', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), satazi)
      ios = NF_INQ_VARID (ncid, 't2m', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), t2m)
      ios = NF_INQ_VARID (ncid, 'mr2m', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), mr2m)
      ios = NF_INQ_VARID (ncid, 'u10', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), u10)
      ios = NF_INQ_VARID (ncid, 'v10', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), v10)
      ios = NF_INQ_VARID (ncid, 'ps', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), ps)
      ios = NF_INQ_VARID (ncid, 'ts', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), ts)
      ios = NF_INQ_VARID (ncid, 'smois', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), smois)
      ios = NF_INQ_VARID (ncid, 'tslb', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), tslb)
      ios = NF_INQ_VARID (ncid, 'snowh', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), snowh)
      ios = NF_INQ_VARID (ncid, 'isflg', varid)
      ios = NF_PUT_VARA_INT(ncid,  varid, istart(2), icount(2), isflg)
      ios = NF_INQ_VARID (ncid, 'soiltyp', varid)
      ios = NF_PUT_VARA_INT(ncid,  varid, istart(2), icount(2), soiltyp)
      ios = NF_INQ_VARID (ncid, 'vegtyp', varid)
      ios = NF_PUT_VARA_INT(ncid,  varid, istart(2), icount(2), vegtyp)
      ios = NF_INQ_VARID (ncid, 'vegfra', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), vegfra)
      ios = NF_INQ_VARID (ncid, 'elev', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), elev)
      ios = NF_INQ_VARID (ncid, 'clwp', varid)
      ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), clwp)
      if ( amsr2 ) then
         ios = NF_INQ_VARID (ncid, 'ret_clw', varid)
         ios = NF_PUT_VARA_REAL(ncid,  varid, istart(2), icount(2), ret_clw)
      end if
      
      ios = NF_CLOSE(ncid)

      deallocate ( npixel )
      deallocate ( iunit )
      deallocate ( inpname )

      deallocate ( ichan )
      deallocate ( datestr2 )
      deallocate ( scanpos )
      deallocate ( landsea_mask )
      deallocate ( elv )
      deallocate ( lat )
      deallocate ( lon )
      deallocate ( satzen )
      deallocate ( satazi )
      deallocate ( t2m )
      deallocate ( mr2m )
      deallocate ( u10 )
      deallocate ( v10 )
      deallocate ( ps )
      deallocate ( ts )
      deallocate ( smois )
      deallocate ( tslb )
      deallocate ( snowh )
      deallocate ( isflg )
      deallocate ( soiltyp )
      deallocate ( vegtyp )
      deallocate ( vegfra )
      deallocate ( elev )
      deallocate ( clwp )
      deallocate ( ret_clw )
      deallocate ( tb_obs )
      deallocate ( tb_bak )
      deallocate ( tb_inv )
      deallocate ( tb_oma )
      deallocate ( ems )
      if ( jac_found ) deallocate ( ems_jac )
      deallocate ( tb_err )
      deallocate ( tb_qc )
      if ( prf_found .and. (rtm_option == 'CRTM') ) then
         deallocate ( prf_pfull )
         deallocate ( prf_phalf )
         deallocate ( prf_t )
         deallocate ( prf_q )
         deallocate ( prf_water )
         deallocate ( prf_ice )
         deallocate ( prf_rain )
         deallocate ( prf_snow )
         deallocate ( prf_grau )
         deallocate ( prf_hail )
         deallocate ( prf_water_reff )
         deallocate ( prf_ice_reff )
         deallocate ( prf_rain_reff )
         deallocate ( prf_snow_reff )
         deallocate ( prf_grau_reff )
         deallocate ( prf_hail_reff )
         if ( jac_found ) then
            deallocate ( prf_t_jac )
            deallocate ( prf_q_jac )
            deallocate ( prf_der_trans )
            deallocate ( prf_trans_jac )
            deallocate ( prf_trans )
            deallocate ( prf_lod_jac )
            deallocate ( prf_lod )
            deallocate ( prf_water_jac )
            deallocate ( prf_ice_jac )
            deallocate ( prf_rain_jac )
            deallocate ( prf_snow_jac )
            deallocate ( prf_grau_jac )
            deallocate ( prf_hail_jac )
            deallocate ( prf_water_reff_jac )
            deallocate ( prf_ice_reff_jac )
            deallocate ( prf_rain_reff_jac )
            deallocate ( prf_snow_reff_jac )
            deallocate ( prf_grau_reff_jac )
            deallocate ( prf_hail_reff_jac )
         end if
      end if
      if ( prf_found .and. (rtm_option == 'RTTOV') ) then
         deallocate ( rtm_prf_p )
         deallocate ( rtm_prf_t )
         deallocate ( rtm_prf_q )
         deallocate ( mdl_prf_p )
         deallocate ( mdl_prf_t )
         deallocate ( mdl_prf_q )
         deallocate ( mdl_prf_qcw )
         deallocate ( mdl_prf_qrn )
      end if

   end do ninst_loop

end do ntime_loop

deallocate ( datestr1 )

end program da_rad_diags

subroutine advance_cymdh(currentdate,dh,newdate)
   
   implicit none       
      
   character(len=10), intent(in)  :: currentdate
   integer,           intent(in)  :: dh
   character(len=10), intent(out) :: newdate
   
   integer :: ccyy, mm, dd, hh

   read(currentdate(1:10), fmt='(i4, 3i2)')  ccyy, mm, dd, hh
   hh = hh + dh
   do while (hh < 0)
      hh = hh + 24
      call change_date ( ccyy, mm, dd, -1 )
   end do  
   do while (hh > 23)                     
      hh = hh - 24                        
      call change_date ( ccyy, mm, dd, 1 )
   end do
   write(newdate,'(i4.4,3(i2.2))') ccyy, mm, dd, hh
   
end subroutine advance_cymdh

subroutine change_date( ccyy, mm, dd, delta )

      implicit none

      integer, intent(inout) :: ccyy, mm, dd
      integer, intent(in)    :: delta
      integer, dimension(12) :: mmday

      mmday = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      mmday(2) = 28
      if (mod(ccyy,4) == 0) then
         mmday(2) = 29
         if ( mod(ccyy,100) == 0) then
            mmday(2) = 28
         endif
         if(mod(ccyy,400) == 0) then
            mmday(2) = 29
         end if
      endif
      dd = dd + delta
      if(dd == 0) then
         mm = mm - 1
         if(mm == 0) then
            mm = 12
            ccyy = ccyy - 1
         endif
         dd = mmday(mm)
      elseif ( dd .gt. mmday(mm) ) then
         dd = 1
         mm = mm + 1
         if(mm > 12 ) then
            mm = 1
            ccyy = ccyy + 1
         end if
      end if

end subroutine change_date
