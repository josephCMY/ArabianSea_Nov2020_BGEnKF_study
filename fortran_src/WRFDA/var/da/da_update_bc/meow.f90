!  program to update BC file from 3dvar output.
!  current version reads only wrf-netcdf file format

program update_wrf_bc

   USE module_netcdf_interface
   USE module_couple_uv

   implicit none

   integer, parameter :: max_3d_variables = 20, &
                         max_2d_variables = 20
 
   character(len=512) :: wrf_3dvar_output_file, &
                         wrf_bdy_file, &
                         wrf_bdy_file_real, &
                         wrf_input_from_si, &
                         wrf_input_from_si_randmean, &
                         wrf_3dvar_random_draw
 
   character(len=20) :: var_pref, var_name, vbt_name

   character(len=20) :: var3d(max_3d_variables), &
                        var2d(max_2d_variables)

   character(len=10), dimension(4) :: bdyname, tenname

   integer           :: ids, ide, jds, jde, kds, kde
   integer           :: num3d, num2d, ndims
   integer           :: time_level
   integer           :: i,j,k,l,m,n,n_1

   integer, dimension(4) :: dims
 
   real, allocatable, dimension(:,:,:) :: tend3d, scnd3d, frst3d, full3d
   real, allocatable, dimension(:,:,:) :: scnd3d_p, frst3d_p, frst3d_r, frst3d_re, &
                                          full3d_r, full3d_0, full3d_1

   real, allocatable, dimension(:,:,:) :: u, v, u_0, v_0, u_1, v_1, u_r, v_r

   real, allocatable, dimension(:,  :) :: mu, mub, msfu, msfv, msfm, &
                                          tend2d, scnd2d, frst2d, full2d

   real, allocatable, dimension(:,  :) :: mu_0,mu_1,mu_r

   real, allocatable, dimension(:,  :) :: frst2d_p, scnd2d_p, frst2d_r, frst2d_re

   real, allocatable, dimension(:,  :) :: tsk, tsk_3dvar

   integer, allocatable, dimension(:,:) :: ivgtyp

   character(len=80), allocatable, dimension(:) :: times, time_var, &
                                                   thisbdytime, nextbdytime
 
   integer :: east_end, north_end, io_status

   logical :: cycling, debug, low_bdy_only, perturb_bdy

   real :: bdyfrq, a

   integer, parameter :: namelist_unit = 7, &
                         ori_unit = 11, &
                         new_unit = 12

   namelist /control_param/ wrf_3dvar_output_file, &
                            wrf_bdy_file, &
                            wrf_bdy_file_real, &
                            wrf_input_from_si, &
                            wrf_input_from_si_randmean, &
                            wrf_3dvar_random_draw, &
                            cycling, debug, low_bdy_only, perturb_bdy, &
                            n_1

!---------------------------------------------------------------------
   wrf_3dvar_output_file = 'wrfinput_d01_update'
   wrf_bdy_file          = 'wrfbdy_d01_update'
   wrf_bdy_file_real     = 'wrfbdy_d01_pri'
   wrf_input_from_si     = 'wrfinput_d01'
   wrf_3dvar_random_draw = 'random_draw'
   wrf_input_from_si_randmean = 'wrfinput_d01_randmean'

   cycling = .false.
   debug   = .false. 
   low_bdy_only = .false.
   perturb_bdy = .false.
   n_1=1
!---------------------------------------------------------------------
!  Read namelist
!----------------------------------------------------------------------------
   io_status = 0

   open(unit = namelist_unit, file = 'parame.in', &
          status = 'old' , access = 'sequential', &
          form   = 'formatted', action = 'read', &
          iostat = io_status )

   if(io_status /= 0) then
      print *, 'Error to open namelist file: parame.in.'
      print *, 'Will work for updating lateral boundary only.'
   else
      read(unit=namelist_unit, nml = control_param , iostat = io_status)

      if(io_status /= 0) then
         print *, 'Error to read control_param. Stopped.'
         stop
      endif

      WRITE(unit=*, fmt='(2a)') &
           'wrf_3dvar_output_file = ', trim(wrf_3dvar_output_file), &
           'wrf_bdy_file          = ', trim(wrf_bdy_file), &
           'wrf_bdy_file_real     = ', trim(wrf_bdy_file_real), &
           'wrf_input_from_si     = ', trim(wrf_input_from_si), &
           'wrf_3dvar_random_draw = ', trim(wrf_3dvar_random_draw)

      WRITE(unit=*, fmt='(a, L10)') &
           'cycling = ', cycling

      close(unit=namelist_unit)
   endif
!----------------------------------------------------------------------------

!--3D need update
   num3d=6
   var3d(1)='U'
   var3d(2)='V'
   var3d(3)='W'
   var3d(4)='T'
   var3d(5)='PH'
   var3d(6)='QVAPOR'

!--2D need update
   num2d=9
   var2d(1)='MUB'
   var2d(2)='MU'
   var2d(3)='MAPFAC_U'
   var2d(4)='MAPFAC_V'
   var2d(5)='MAPFAC_M'
   var2d(6)='TMN'
   var2d(7)='SST'
   var2d(8)='TSK'
   var2d(9)='VEGFRA'
!   var2d(10)='ALBBCK'

!---------------------------------------------------------------------

!--First, the boundary times
   call get_dims_cdf(wrf_bdy_file, 'Times', dims, ndims, debug)

   if(debug) then
      write(unit=*, fmt='(a,i2,2x,a,4i6)') &
           'Times: ndims=', ndims, 'dims=', (dims(i), i=1,ndims)
   endif

   time_level = dims(2)

   if(time_level < 1) then
      write(unit=*, fmt='(a,i2/a)') &
           'time_level = ', time_level, &
           'We need at least one time-level BDY.'
      stop 'Wrong BDY file.'
   endif

   allocate(times(dims(2)))
   allocate(thisbdytime(dims(2)))
   allocate(nextbdytime(dims(2)))

   call get_times_cdf(wrf_bdy_file, times, dims(2), dims(2), debug )

   call get_bdytimestr_cdf(wrf_bdy_file, 'thisbdytime', thisbdytime, dims(2), debug )
   call get_bdytimestr_cdf(wrf_bdy_file, 'nextbdytime', nextbdytime, dims(2), debug )

   call get_bdyfrq(thisbdytime(1), nextbdytime(1), bdyfrq, debug)

   if(debug) then
      do n=1, dims(2)
         write(unit=*, fmt='(3(a, i2, 2a,2x))') &
           '       times(', n, ')=', trim(times(n)), &
           'thisbdytime (', n, ')=', trim(thisbdytime(n)), &
           'nextbdytime (', n, ')=', trim(nextbdytime(n))
      enddo
   endif

   !zmeng
   print*, 'processing   ',trim(times(n_1))

!---------------------------------------------------------------------
   east_end=0
   north_end=0
!---------------------------------------------------------------------
!--For 2D variables
!--Get mu, mub, msfu, and msfv

   do n=1,num2d
      call get_dims_cdf( wrf_3dvar_output_file, trim(var2d(n)), dims, ndims, debug )

      select case(trim(var2d(n)))
         case ('MU') ;
            if(low_bdy_only) cycle

            allocate(mu(dims(1), dims(2)))
            call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), mu, &
                                      dims(1), dims(2), 1, debug )

            allocate(mu_0(dims(1), dims(2)))
            call get_var_2d_real_cdf( wrf_input_from_si_randmean, trim(var2d(n)), mu_0, &
                                      dims(1), dims(2), 1, debug )

            allocate(mu_r(dims(1), dims(2)))
            call get_var_2d_real_cdf( wrf_3dvar_random_draw, trim(var2d(n)), mu_r, &
                                      dims(1), dims(2), 1, debug )

            east_end=dims(1)+1
            north_end=dims(2)+1
         case ('MUB') ;
            if(low_bdy_only) cycle

            allocate(mub(dims(1), dims(2)))

            call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), mub, &
                                      dims(1), dims(2), 1, debug )
         case ('MAPFAC_U') ;
            if(low_bdy_only) cycle

            allocate(msfu(dims(1), dims(2)))

            call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), msfu, &
                                      dims(1), dims(2), 1, debug )
         case ('MAPFAC_V') ;
            if(low_bdy_only) cycle

            allocate(msfv(dims(1), dims(2)))

            call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), msfv, &
                                      dims(1), dims(2), 1, debug )
         case ('MAPFAC_M') ;
            if(low_bdy_only) cycle

            allocate(msfm(dims(1), dims(2)))

            call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), msfm, &
                                      dims(1), dims(2), 1, debug )
         case ('TSK') ;
            if(.not. cycling) cycle

            allocate(tsk(dims(1), dims(2)))
            allocate(tsk_3dvar(dims(1), dims(2)))
            allocate(ivgtyp(dims(1), dims(2)))

            call get_var_2d_real_cdf( wrf_input_from_si, trim(var2d(n)), tsk, &
                                      dims(1), dims(2), 1, debug )
            call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), tsk_3dvar, &
                                      dims(1), dims(2), 1, debug )
            call get_var_2d_int_cdf( wrf_3dvar_output_file, 'IVGTYP', ivgtyp, &
                                      dims(1), dims(2), 1, debug )

!-----------update TSK.
            do j=1,dims(2)
            do i=1,dims(1)
               if(ivgtyp(i,j) /= 16) &
                  tsk(i,j)=tsk_3dvar(i,j)
            enddo
            enddo

            call put_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), tsk, &
                                      dims(1), dims(2), 1, debug )
            deallocate(tsk)
            deallocate(ivgtyp)
            deallocate(tsk_3dvar)
         case ('TMN', 'SST', 'VEGFRA', 'ALBBCK') ;
            if(.not. cycling) cycle

            allocate(full2d(dims(1), dims(2)))

            call get_var_2d_real_cdf( wrf_input_from_si, trim(var2d(n)), full2d, &
                                      dims(1), dims(2), 1, debug )

            call put_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), full2d, &
                                      dims(1), dims(2), 1, debug )
            deallocate(full2d)
         case default ;
            print *, 'It is impossible here. var2d(n)=', trim(var2d(n))
      end select
   enddo

!---------------------------------------------------------------------
  
   if(low_bdy_only) then
      print *, 'Only low boundary updated.'
      stop 
   endif

   if(east_end < 1 .or. north_end < 1) then
      write(unit=*, fmt='(a)') 'Wrong data for Boundary.'
      stop
   endif

!---------------------------------------------------------------------
!--boundary variables
   bdyname(1)='_BXS'
   bdyname(2)='_BXE'
   bdyname(3)='_BYS'
   bdyname(4)='_BYE'

!--boundary tendancy variables
   tenname(1)='_BTXS'
   tenname(2)='_BTXE'
   tenname(3)='_BTYS'
   tenname(4)='_BTYE'

!---------------------------------------------------------------------

   do m=1,4
      var_name='MU' // trim(bdyname(m))
      vbt_name='MU' // trim(tenname(m))

      call get_dims_cdf( wrf_bdy_file, trim(var_name), dims, ndims, debug )

      allocate(frst2d(dims(1), dims(2)))
      allocate(frst2d_r(dims(1), dims(2)))
      allocate(frst2d_re(dims(1), dims(2)))
      allocate(scnd2d(dims(1), dims(2)))
      allocate(tend2d(dims(1), dims(2)))

!-----Get variable at second time level
      if(time_level > 1 .and. n_1 .lt. time_level ) then
         call get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), scnd2d, &
                                   dims(1), dims(2), n_1+1, debug )
         if ( n_1 > 1 ) &
            call get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), frst2d, &
                                   dims(1), dims(2), n_1, debug )
      else
         call get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), frst2d, &
                                   dims(1), dims(2), n_1, debug )
         call get_var_2d_real_cdf( wrf_bdy_file, trim(vbt_name), tend2d, &
                                   dims(1), dims(2), n_1, debug )
         if (time_level > 1 .and. n_1 .eq. time_level ) &
         call get_var_2d_real_cdf( wrf_bdy_file_real, trim(var_name), frst2d_re, &
                                   dims(1), dims(2), n_1, debug )
      endif

      if(debug) then
         write(unit=ori_unit, fmt='(a,i2,2x,2a/a,i2,2x,a,4i6)') &
              'No.', m, 'Variable: ', trim(vbt_name), &
              'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)

         call get_var_2d_real_cdf( wrf_bdy_file, trim(vbt_name), tend2d, &
                                   dims(1), dims(2), n_1, debug )

         write(unit=ori_unit, fmt='(a, 10i12)') &
              ' old ', (i, i=1,dims(2))
         do j=1,dims(1)
            write(unit=ori_unit, fmt='(i4, 1x, 5e20.7)') &
                  j, (tend2d(j,i), i=1,dims(2))
         enddo
      endif

!-----calculate variable at first time level
      select case(m)
         case (1) ;		! West boundary
            do l=1,dims(2)
            do j=1,dims(1)
               if(time_level < 2 ) &
               scnd2d(j,l)=frst2d(j,l)+tend2d(j,l)*bdyfrq
               if(time_level > 1 .and. n_1 .eq. time_level) &
               scnd2d(j,l)=frst2d_re(j,l)+tend2d(j,l)*bdyfrq
               if (n_1 < 2) frst2d(j,l)=mu(l,j)
               frst2d_r(j,l)=mu_r(l,j)-mu_0(l,j)
            enddo
            enddo
         case (2) ;		! East boundary
            do l=1,dims(2)
            do j=1,dims(1)
               if(time_level < 2 ) &
               scnd2d(j,l)=frst2d(j,l)+tend2d(j,l)*bdyfrq
               if(time_level > 1 .and. n_1 .eq. time_level) &
               scnd2d(j,l)=frst2d_re(j,l)+tend2d(j,l)*bdyfrq
               if (n_1 < 2) frst2d(j,l)=mu(east_end-l,j)
               frst2d_r(j,l)=mu_r(east_end-l,j)-mu_0(east_end-l,j)
            enddo
            enddo
         case (3) ;		! South boundary
            do l=1,dims(2)
            do i=1,dims(1)
               if(time_level < 2 ) &
               scnd2d(i,l)=frst2d(i,l)+tend2d(i,l)*bdyfrq
               if(time_level > 1 .and. n_1 .eq. time_level) &
               scnd2d(i,l)=frst2d_re(i,l)+tend2d(i,l)*bdyfrq
               if (n_1 < 2) frst2d(i,l)=mu(i,l)
               frst2d_r(i,l)=mu_r(i,l)-mu_0(i,l)
            enddo
            enddo
         case (4) ;		! North boundary
            do l=1,dims(2)
            do i=1,dims(1)
               if(time_level < 2 ) &
               scnd2d(i,l)=frst2d(i,l)+tend2d(i,l)*bdyfrq
               if(time_level > 1 .and. n_1 .eq. time_level) &
               scnd2d(i,l)=frst2d_re(i,l)+tend2d(i,l)*bdyfrq
               if (n_1 < 2) frst2d(i,l)=mu(i,north_end-l)
               frst2d_r(i,l)=mu_r(i,north_end-l)-mu_0(i,north_end-l)
            enddo
            enddo
         case default ;
            print *, 'It is impossible here. mu, m=', m
      end select

!-----calculate perturbation of mu  at second time level-zmeng
      if ( perturb_bdy ) then

      do l=1,dims(2)
      do i=1,dims(1)
!         scnd2d_p(i,l)=a*frst2d_p(i,l)+sqrt(1-a**2)*frst2d_r(i,l)
         scnd2d(i,l)=scnd2d(i,l)+frst2d_r(i,l)
      enddo
      enddo
      endif

!-----calculate new tendancy 
      do l=1,dims(2)
      do i=1,dims(1)
         tend2d(i,l)=(scnd2d(i,l)-frst2d(i,l))/bdyfrq
      enddo
      enddo

      if(debug) then
         write(unit=new_unit, fmt='(a,i2,2x,2a/a,i2,2x,a,4i6)') &
              'No.', m, 'Variable: ', trim(vbt_name), &
              'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)

         write(unit=new_unit, fmt='(a, 10i12)') &
              ' new ', (i, i=1,dims(2))

         do j=1,dims(1)
            write(unit=new_unit, fmt='(i4, 1x, 5e20.7)') &
                  j, (tend2d(j,i), i=1,dims(2))
         enddo
      endif

!-----output new variable at first time level
      if ( n_1 < 2 ) &
      call put_var_2d_real_cdf( wrf_bdy_file, trim(var_name), frst2d, &
                                dims(1), dims(2), n_1, debug )

!-----output new variable at second time level
      if ( time_level > 1  .and. n_1 .lt. time_level ) then
      call put_var_2d_real_cdf( wrf_bdy_file, trim(var_name), scnd2d, &
                                dims(1), dims(2), n_1+1, debug )
      endif 
!-----output new tendancy 
      call put_var_2d_real_cdf( wrf_bdy_file, trim(vbt_name), tend2d, &
                                dims(1), dims(2), n_1, debug )

      if (debug) then
        print*,'m, frst2d(2,2)',m, frst2d(2,2)
        print*,'m, frst2d_r(2,2)',m, frst2d_r(2,2)
        print*,'m, scnd2d(2,2)',m, scnd2d(2,2)
      endif

      deallocate(frst2d)
      deallocate(frst2d_r)
      deallocate(frst2d_re)
      deallocate(scnd2d)
      deallocate(tend2d)
   enddo

!---------------------------------------------------------------------
!--For 3D variables

!--Get U
   call get_dims_cdf( wrf_3dvar_output_file, 'U', dims, ndims, debug )

!  call get_att_cdf( wrf_3dvar_output_file, 'U', debug )

   allocate(u(dims(1), dims(2), dims(3)))
   allocate(u_0(dims(1), dims(2), dims(3)))
   allocate(u_r(dims(1), dims(2), dims(3)))

   ids=1
   ide=dims(1)-1
   jds=1
   jde=dims(2)
   kds=1
   kde=dims(3)

   call get_var_3d_real_cdf( wrf_3dvar_output_file, 'U', u, &
                             dims(1), dims(2), dims(3), 1, debug )

   call get_var_3d_real_cdf( wrf_input_from_si_randmean, 'U', u_0, &
                             dims(1), dims(2), dims(3), 1, debug )

   call get_var_3d_real_cdf( wrf_3dvar_random_draw, 'U', u_r, &
                             dims(1), dims(2), dims(3), 1, debug )

!  do j=1,dims(2)
!     write(unit=*, fmt='(2(a,i5), a, f12.8)') &
!          'u(', dims(1), ',', j, ',1)=', u(dims(1),j,1)
!  enddo

!--Get V
   call get_dims_cdf( wrf_3dvar_output_file, 'V', dims, ndims, debug )

!  call get_att_cdf( wrf_3dvar_output_file, 'V', debug )

   allocate(v(dims(1), dims(2), dims(3)))
   allocate(v_0(dims(1), dims(2), dims(3)))
   allocate(v_r(dims(1), dims(2), dims(3)))

   call get_var_3d_real_cdf( wrf_3dvar_output_file, 'V', v, &
                             dims(1), dims(2), dims(3), 1, debug )
   call get_var_3d_real_cdf( wrf_input_from_si_randmean, 'V', v_0, &
                             dims(1), dims(2), dims(3), 1, debug )
   call get_var_3d_real_cdf( wrf_3dvar_random_draw, 'V', v_r, &
                             dims(1), dims(2), dims(3), 1, debug )

!  do i=1,dims(1)
!     write(unit=*, fmt='(2(a,i5), a, f12.8)') &
!          'v(', i, ',', dims(2), ',1)=', v(i,dims(2),1)
!  enddo

   if(debug) then
      write(unit=*, fmt='(a,e20.12,4x)') &
           'Before couple Sample u=', u(dims(1)/2,dims(2)/2,dims(3)/2), &
           'Before couple Sample v=', v(dims(1)/2,dims(2)/2,dims(3)/2)
   endif

!---------------------------------------------------------------------
!--Couple u, v.
   call couple_uv ( u, v, mu, mub, msfu, msfv, ids, ide, jds, jde, kds, kde )
   call couple_uv ( u_0, v_0, mu_0, mub, msfu, msfv, ids, ide, jds, jde, kds, kde )
   call couple_uv ( u_r, v_r, mu_r, mub, msfu, msfv, ids, ide, jds, jde, kds, kde )

   if(debug) then
      write(unit=*, fmt='(a,e20.12,4x)') &
           'After  couple Sample u=', u(dims(1)/2,dims(2)/2,dims(3)/2), &
           'After  couple Sample v=', v(dims(1)/2,dims(2)/2,dims(3)/2)
   endif

!---------------------------------------------------------------------
!--For 3D variables

   do n=1,num3d
      write(unit=*, fmt='(a, i3, 2a)') 'Processing: var3d(', n, ')=', trim(var3d(n))

      call get_dims_cdf( wrf_3dvar_output_file, trim(var3d(n)), dims, ndims, debug )

      allocate(full3d(dims(1), dims(2), dims(3)))
      allocate(full3d_0(dims(1), dims(2), dims(3)))
      allocate(full3d_r(dims(1), dims(2), dims(3)))

      east_end=dims(1)+1
      north_end=dims(2)+1

      select case(trim(var3d(n)))
         case ('U') ;		! U
!           var_pref='R' // trim(var3d(n))
            var_pref=trim(var3d(n))
            full3d(:,:,:)=u(:,:,:)
            full3d_0(:,:,:)=u_0(:,:,:)
            full3d_r(:,:,:)=u_r(:,:,:)
         case ('V') ;		! V
!           var_pref='R' // trim(var3d(n))
            var_pref=trim(var3d(n))
            full3d(:,:,:)=v(:,:,:)
            full3d_0(:,:,:)=v_0(:,:,:)
            full3d_r(:,:,:)=v_r(:,:,:)
         case ('W') ;
!           var_pref = 'R' // trim(var3d(n))
            var_pref = trim(var3d(n))

            call get_var_3d_real_cdf( wrf_3dvar_output_file, trim(var3d(n)), full3d, &
                                      dims(1), dims(2), dims(3), 1, debug )
            call get_var_3d_real_cdf( wrf_input_from_si_randmean, trim(var3d(n)), full3d_0, &
                                      dims(1), dims(2), dims(3), 1, debug )
            call get_var_3d_real_cdf( wrf_3dvar_random_draw, trim(var3d(n)), full3d_r, &
                                      dims(1), dims(2), dims(3), 1, debug )

            if(debug) then
               write(unit=*, fmt='(3a,e20.12,4x)') &
                    'Before couple Sample ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            endif

            do k=1,dims(3)
            do j=1,dims(2)
            do i=1,dims(1)
               full3d(i,j,k)=full3d(i,j,k)*(mu(i,j)+mub(i,j))/msfm(i,j)
               full3d_0(i,j,k)=full3d_0(i,j,k)*(mu_0(i,j)+mub(i,j))/msfm(i,j)
               full3d_r(i,j,k)=full3d_r(i,j,k)*(mu_r(i,j)+mub(i,j))/msfm(i,j)
            enddo
            enddo
            enddo

            if(debug) then
               write(unit=*, fmt='(3a,e20.12,4x)') &
                    'After  couple Sample ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            endif
         case ('T', 'PH') ;
            var_pref=trim(var3d(n))
 
            call get_var_3d_real_cdf( wrf_3dvar_output_file, trim(var3d(n)), full3d, &
                                      dims(1), dims(2), dims(3), 1, debug )
            call get_var_3d_real_cdf( wrf_input_from_si_randmean, trim(var3d(n)), full3d_0, &
                                      dims(1), dims(2), dims(3), 1, debug )
            call get_var_3d_real_cdf( wrf_3dvar_random_draw, trim(var3d(n)), full3d_r, &
                                      dims(1), dims(2), dims(3), 1, debug )

            if(debug) then
               write(unit=*, fmt='(3a,e20.12,4x)') &
                    'Before couple Sample ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            endif

            do k=1,dims(3)
            do j=1,dims(2)
            do i=1,dims(1)
               full3d(i,j,k)=full3d(i,j,k)*(mu(i,j)+mub(i,j))
               full3d_0(i,j,k)=full3d_0(i,j,k)*(mu_0(i,j)+mub(i,j))
               full3d_r(i,j,k)=full3d_r(i,j,k)*(mu_r(i,j)+mub(i,j))
            enddo
            enddo
            enddo

            if(debug) then
               write(unit=*, fmt='(3a,e20.12,4x)') &
                    'After  couple Sample ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            endif
         case ('QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP') ;
!           var_pref='R' // var3d(n)(1:2)
!           var_pref=var3d(n)(1:2)
            var_pref=var3d(n)
 
            call get_var_3d_real_cdf( wrf_3dvar_output_file, trim(var3d(n)), full3d, &
                                      dims(1), dims(2), dims(3), 1, debug )
            call get_var_3d_real_cdf( wrf_input_from_si_randmean, trim(var3d(n)), full3d_0, &
                                      dims(1), dims(2), dims(3), 1, debug )
            call get_var_3d_real_cdf( wrf_3dvar_random_draw, trim(var3d(n)), full3d_r, &
                                      dims(1), dims(2), dims(3), 1, debug )

            if(debug) then
               write(unit=*, fmt='(3a,e20.12,4x)') &
                    'Before couple Sample ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            endif

            do k=1,dims(3)
            do j=1,dims(2)
            do i=1,dims(1)
               full3d(i,j,k)=full3d(i,j,k)*(mu(i,j)+mub(i,j))
               full3d_0(i,j,k)=full3d_0(i,j,k)*(mu_0(i,j)+mub(i,j))
               full3d_r(i,j,k)=full3d_r(i,j,k)*(mu_r(i,j)+mub(i,j))
            enddo
            enddo
            enddo

            if(debug) then
               write(unit=*, fmt='(3a,e20.12,4x)') &
                    'After  couple Sample ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            endif
         case default ;
            print *, 'It is impossible here. var3d(', n, ')=', trim(var3d(n))
      end select

      do m=1,4
         var_name=trim(var_pref) // trim(bdyname(m))
         vbt_name=trim(var_pref) // trim(tenname(m))

         write(unit=*, fmt='(a, i3, 2a)') 'Processing: bdyname(', m, ')=', trim(var_name)

         call get_dims_cdf( wrf_bdy_file, trim(var_name), dims, ndims, debug )

         allocate(frst3d(dims(1), dims(2), dims(3)))
         allocate(frst3d_r(dims(1), dims(2), dims(3)))
         allocate(frst3d_re(dims(1), dims(2), dims(3)))
         allocate(scnd3d(dims(1), dims(2), dims(3)))
         allocate(tend3d(dims(1), dims(2), dims(3)))

!--------Get variable at second time level
         if(time_level > 1  .and. n_1 .lt. time_level ) then
            call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), scnd3d, &
                                      dims(1), dims(2), dims(3), n_1+1, debug )
            if( n_1 > 1 ) &
            call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                                      dims(1), dims(2), dims(3), n_1, debug )
         else
            call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                                      dims(1), dims(2), dims(3), n_1, debug )
            call get_var_3d_real_cdf( wrf_bdy_file, trim(vbt_name), tend3d, &
                                      dims(1), dims(2), dims(3), n_1, debug )
            if (time_level > 1 .and. n_1 .eq. time_level ) & 
            call get_var_3d_real_cdf( wrf_bdy_file_real, trim(var_name), frst3d_re, &
                                      dims(1), dims(2), dims(3), n_1, debug )
         endif

         if(debug) then
            write(unit=ori_unit, fmt='(a,i2,2x,2a/a,i2,2x,a,4i6)') &
                 'No.', m, 'Variable: ', trim(vbt_name), &
                 'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)

            call get_var_3d_real_cdf( wrf_bdy_file, trim(vbt_name), tend3d, &
                                      dims(1), dims(2), dims(3), n_1, debug )

            write(unit=ori_unit, fmt='(a, 10i12)') &
                 ' old ', (i, i=1,dims(3))
            do j=1,dims(1)
               write(unit=ori_unit, fmt='(i4, 1x, 5e20.7)') &
                     j, (tend3d(j,dims(2)/2,i), i=1,dims(3))
            enddo
         endif
   
         select case(trim(bdyname(m)))
            case ('_BXS') ;		! West boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do j=1,dims(1)
                  if(time_level < 2 ) &
                  scnd3d(j,k,l)=frst3d(j,k,l)+tend3d(j,k,l)*bdyfrq
                  if (time_level > 1 .and. n_1 .eq. time_level ) &
                  scnd3d(j,k,l)=frst3d_re(j,k,l)+tend3d(j,k,l)*bdyfrq
                  if (n_1 < 2) frst3d(j,k,l)=full3d(l,j,k)
                  frst3d_r(j,k,l)=full3d_r(l,j,k)-full3d_0(l,j,k)
               enddo
               enddo
               enddo
            case ('_BXE') ;		! East boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do j=1,dims(1)
                  if(time_level < 2 ) &
                  scnd3d(j,k,l)=frst3d(j,k,l)+tend3d(j,k,l)*bdyfrq
                  if (time_level > 1 .and. n_1 .eq. time_level ) &
                  scnd3d(j,k,l)=frst3d_re(j,k,l)+tend3d(j,k,l)*bdyfrq
                  if (n_1 < 2) frst3d(j,k,l)=full3d(east_end-l,j,k)
                  frst3d_r(j,k,l)=full3d_r(east_end-l,j,k)-full3d_0(east_end-l,j,k)
               enddo
               enddo
               enddo
            case ('_BYS') ;		! South boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do i=1,dims(1)
                  if(time_level < 2 ) &
                  scnd3d(i,k,l)=frst3d(i,k,l)+tend3d(i,k,l)*bdyfrq
                  if (time_level > 1 .and. n_1 .eq. time_level ) &
                  scnd3d(i,k,l)=frst3d_re(i,k,l)+tend3d(i,k,l)*bdyfrq
                  if (n_1 < 2) frst3d(i,k,l)=full3d(i,l,k)
                  frst3d_r(i,k,l)=full3d_r(i,l,k)-full3d_0(i,l,k)
               enddo
               enddo
               enddo
            case ('_BYE') ;		! North boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do i=1,dims(1)
                  if(time_level < 2 ) &
                  scnd3d(i,k,l)=frst3d(i,k,l)+tend3d(i,k,l)*bdyfrq
                  if (time_level > 1 .and. n_1 .eq. time_level ) &
                  scnd3d(i,k,l)=frst3d_re(i,k,l)+tend3d(i,k,l)*bdyfrq
                  if (n_1 < 2) frst3d(i,k,l)=full3d(i,north_end-l,k)
                  frst3d_r(i,k,l)=full3d_r(i,north_end-l,k)-full3d_0(i,north_end-l,k)
               enddo
               enddo
               enddo
            case default ;
               print *, 'It is impossible here.'
               print *, 'bdyname(', m, ')=', trim(bdyname(m))
               stop
         end select

         write(unit=*, fmt='(a, i3, 2a)') 'cal. tend: bdyname(', m, ')=', trim(vbt_name)

!--------calculate new tendancy
         if  ( perturb_bdy ) then

         do l=1,dims(3)
         do k=1,dims(2)
         do i=1,dims(1)
!           scnd3d_p(i,k,l)=a*frst3d_p(i,k,l)+sqrt(1-a**2)*frst3d_r(i,k,l)
           scnd3d(i,k,l)=scnd3d(i,k,l)+frst3d_r(i,k,l)
         enddo
         enddo
         enddo
         
         endif

!--------calculate new tendancy 
         do l=1,dims(3)
         do k=1,dims(2)
         do i=1,dims(1)
            tend3d(i,k,l)=(scnd3d(i,k,l)-frst3d(i,k,l))/bdyfrq
         enddo
         enddo
         enddo

         if(debug) then
            write(unit=new_unit, fmt='(a,i2,2x,2a/a,i2,2x,a,4i6)') &
                 'No.', m, 'Variable: ', trim(vbt_name), &
                 'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)

            write(unit=new_unit, fmt='(a, 10i12)') &
                 ' new ', (i, i=1,dims(3))

            do j=1,dims(1)
               write(unit=new_unit, fmt='(i4, 1x, 5e20.7)') &
                     j, (tend3d(j,dims(2)/2,i), i=1,dims(3))
            enddo
         endif

!--------output new variable at first time level

         if ( n_1 < 2 ) &
         call put_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                                dims(1), dims(2), dims(3), n_1, debug )
         if (time_level > 1  .and. n_1 .lt. time_level ) then
         call put_var_3d_real_cdf( wrf_bdy_file, trim(var_name), scnd3d, &
                                dims(1), dims(2), dims(3), n_1+1, debug )
         endif
         call put_var_3d_real_cdf( wrf_bdy_file, trim(vbt_name), tend3d, &
                                   dims(1), dims(2), dims(3), n_1, debug )

         deallocate(frst3d)
         deallocate(frst3d_r)
         deallocate(frst3d_re)
         deallocate(scnd3d)
         deallocate(tend3d)
      enddo
      
      deallocate(full3d)
      deallocate(full3d_r)
      deallocate(full3d_0)
   enddo

   deallocate(mu)
   deallocate(mu_r)
   deallocate(mu_0)
   deallocate(mub)
   deallocate(msfu)
   deallocate(msfv)
   deallocate(u)
   deallocate(u_0)
   deallocate(u_r)
   deallocate(v)
   deallocate(v_0)
   deallocate(v_r)

   if(io_status /= 0) then
      print *, 'Only lateral boundary updated.'
   else
      if(low_bdy_only) then
         if(cycling) then
            print *, 'Both low boudary and lateral boundary updated.'
         else
            print *, 'Only low boudary updated.'
         endif
      else
         print *, 'Only lateral boundary updated.'
      endif
   endif
   print *, '*** Update_bc completed successfully ***'
end program update_wrf_bc


