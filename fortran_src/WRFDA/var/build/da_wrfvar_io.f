












module da_wrfvar_io

   use module_configure, only : grid_config_rec_type, model_config_rec, &
      model_to_grid_config_rec
   use module_date_time, only : get_julgmt, geth_julgmt, current_date, start_date
   use module_domain, only : domain, get_ijk_from_grid
   use module_io_domain, only : open_r_dataset,close_dataset, &
      input_input, open_w_dataset,output_input, &
      input_boundary, output_boundary, output_auxhist4, &
      input_auxhist6, input_auxhist4
   use module_state_description, only : p_qv

   use da_control, only : trace_use, ierr, var4d, var4d_lbc, num_fgat_time, rootproc
   use da_reporting, only : da_error, message, da_message
   use da_tracing, only : da_trace_entry, da_trace_exit, da_trace

   use, intrinsic :: iso_c_binding,                       &
                     ONLY: c_int32_t, C_CHAR, C_NULL_CHAR



contains

subroutine da_med_initialdata_input (grid, config_flags, filename, in_date)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)                :: grid
   type (grid_config_rec_type), intent(inout) :: config_flags

   character(*),  intent (in)                 :: filename
   character(*),  intent (in),  optional      :: in_date

   integer                 :: fid , status, n, nsave

   integer :: julyr, julday
   real    :: gmt

   if (trace_use) call da_trace_entry("da_med_initialdata_input")
   if (trace_use) call da_trace("da_med_initialdata_input", &
      Message="Reading "//trim(filename))

   ! Initialize the mother domain.

   grid%input_from_file = .true.
   ! Initialize the array grid%pb (YRG, 08/26/2010) for use of wrftoxb:
   grid%pb = 0.0
   
   call ext_ncd_open_for_read(trim(filename), 0, 0, "", fid, ierr)

   if (ierr /= 0) then
      write(unit=message(1), fmt='(2a)') &
         'Netcdf error opening file:', trim(filename)
      call da_error("da_med_initialdata_input.inc",35,message(1:1))
   end if

   call ext_ncd_get_next_time(fid, current_date, Status)

   if (present(in_date)) then
      ! Set start_date to current_date.
      read(in_date(1:19), fmt='(i4, 5(1x, i2))') &
         grid%start_year,   &
         grid%start_month,  &
         grid%start_day,    &
         grid%start_hour,   &
         grid%start_minute, &
         grid%start_second 

      nsave = -1
      do n=1, 1000
         if (current_date(1:19) == in_date(1:19)) then
            nsave = n - 1
            exit
         end if
         call ext_ncd_get_next_time(fid, current_date, Status)
      end do

      if (nsave < 0) then
         call da_error("da_med_initialdata_input.inc",60,(/"Cannot find the needed time"/))
      end if
   else
      ! Set start_date to current_date.
      read(current_date(1:19), fmt='(i4, 5(1x, i2))') &
           grid%start_year,  &
           grid%start_month, &
           grid%start_day,   &
           grid%start_hour,  &
           grid%start_minute,&
           grid%start_second
   end if

   call geth_julgmt(julyr, julday, gmt)
   call nl_set_gmt (grid%id, gmt)
   call nl_set_julyr (grid%id, julyr)
   call nl_set_julday (grid%id, julday)

   call nl_set_iswater (grid%id, grid%iswater)
   call nl_set_cen_lat (grid%id , grid%cen_lat)
   call nl_set_cen_lon (grid%id , grid%cen_lon)
   call nl_set_truelat1 (grid%id , grid%truelat1)
   call nl_set_truelat2 (grid%id , grid%truelat2)
   call nl_set_moad_cen_lat (grid%id , grid%moad_cen_lat)
   call nl_set_stand_lon (grid%id , grid%stand_lon)
   call nl_set_pole_lat (grid%id , grid%pole_lat)
   call nl_set_map_proj (grid%id , grid%map_proj)
   start_date=current_date

   call geth_julgmt(julyr, julday, gmt)
   config_flags%gmt = gmt
   config_flags%julyr = julyr
   config_flags%julday = julday

   call ext_ncd_ioclose(fid, ierr)

   call da_trace("da_med_initialdata_input", &
       message="open_r_dataset for "//trim(filename))
   call open_r_dataset (fid, trim(filename), grid , config_flags , &
      "DATASET=INPUT", ierr)

   if (ierr .NE. 0) then
      write(unit=message(1),fmt='(A,A,A,I5)') 'Error opening ', &
        trim(filename),' for reading ierr=',ierr
      call da_error("da_med_initialdata_input.inc",104,message(1:1))
   end if

   if (present(in_date)) then
      do n=1, nsave
         call da_message((/"current_date="//current_date// &
            ', in_date='//in_date/))
         call ext_ncd_get_next_time(fid, current_date, Status)
      end do
   end if

   call input_input (fid ,   grid , config_flags , ierr)

   call nl_get_mminlu (grid%id , grid%mminlu)

   call close_dataset (fid , config_flags , "DATASET=INPUT")

   if (trace_use) call da_trace_exit("da_med_initialdata_input")

end subroutine da_med_initialdata_input


subroutine da_med_initialdata_output (grid , config_flags, out_filename)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(domain), intent(in)                   :: grid
   type (grid_config_rec_type) , intent(inout) :: config_flags       
   character(*),  intent (in),  optional       :: out_filename

   integer                :: fid
   character (len=80)     :: file_name

   integer :: julyr, julday
   real    :: gmt

   if (trace_use) call da_trace_entry("da_med_initialdata_output")

   if (present(out_filename)) then
      file_name = trim(out_filename)
   else
      file_name = 'wrfvar_output'
   end if

   call da_trace ("da_med_initialdata_ouput",Message="Writing wrfvar output")
   call open_w_dataset (fid, trim(file_name), grid , config_flags , &
                         output_input , "DATASET=INPUT,REAL_OUTPUT_SIZE=4", ierr)

   if (ierr /= 0) then
      write(unit=message(1),fmt=*) 'Error opening ', &
         trim(file_name),' for writing ierr=',ierr
      call da_error("da_med_initialdata_output.inc",34,message(1:1))
   end if

   start_date=current_date

   call geth_julgmt(julyr, julday, gmt)
   config_flags%gmt = gmt
   config_flags%julyr = julyr
   config_flags%julday = julday

   call output_input (fid, grid , config_flags , ierr)

   call close_dataset (fid , config_flags, "DATASET=INPUT,REAL_OUTPUT_SIZE=4")

   if (trace_use) call da_trace_exit("da_med_initialdata_output")

end subroutine da_med_initialdata_output


subroutine da_med_initialdata_output_lbc (grid , config_flags, out_filename)

   !-----------------------------------------------------------------------
   ! Purpose: Write out LBC condition at the t=0 for 4DVAR with LBC control
   !          We only use the perturbation of LBC
   ! Xin Zhang: 10/8/2010
   !-----------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)                 :: grid
   type (grid_config_rec_type) , intent(inout) :: config_flags       
   character(*),  intent (in),  optional       :: out_filename

end subroutine da_med_initialdata_output_lbc


subroutine da_update_firstguess(grid, out_filename) 

!---------------------------------------------------------------------------
!  Purpose: Only replce the fields touched by WRFDA, other than re-generate 
!           whole wrfvar_output from the scratch.
!
!  Update  :   jliu@ucar.edu Apr 1, 2013
!              Minor change for Purpose
!              Added copy file mods instead of copying file outside 
!              Requires Fortran 2003 standard compiler
!  Creator :   Junmei Ban, Mar 14, 2013
!         
!  The following WRF fields are modified:
!    grid%u_2
!    grid%v_2
!    grid%w_2
!    grid%mu_2
!    grid%ph_2
!    grid%t_2
!    grid%moist
!    grid%p
!    grid%psfc
!    grid%t2, grid%q2, grid%u10, grid%v10, grid%th2
!
!---------------------------------------------------------------------------

  use module_domain, only : get_ijk_from_grid, program_name
  use da_control,    only : use_radarobs, use_rad, crtm_cloud, &
                            use_radar_rhv, use_radar_rqv    
  use module_state_description, only : p_qv, p_qc, p_qr, p_qi, &
                                       p_qs, p_qg

  implicit none

  INTERFACE
    integer(c_int32_t) function copyfile(ifile, ofile) bind(c)
      import :: c_int32_t, C_CHAR
      CHARACTER(KIND=C_CHAR), DIMENSION(*), intent(in) :: ifile, ofile
    END function copyfile
  END INTERFACE

  include 'netcdf.inc'

  type(domain), intent(in)             :: grid
  character(*), intent(in), optional   :: out_filename

! Declare local parameters
  character(len=120) :: file_name 
  character(len=19)  :: DateStr1
  character(len=4)   :: staggering=' N/A'
  character(len=3)   :: ordering  
  character(len=80), dimension(3)  ::  dimnames
  character(len=31)  :: rmse_var  
  integer            :: dh1
  integer            :: i,j,k
  integer            :: ndim1
  integer            :: WrfType
  integer            :: it, ierr, Status, Status_next_time
  integer            :: wrf_real
  integer            :: nlon_regional,nlat_regional,nsig_regional
  integer            :: ids, ide, jds, jde, kds, kde,    &
                        ims, ime, jms, jme, kms, kme,    &
                        ips, ipe, jps, jpe, kps, kpe    
  integer, dimension(4)           :: start_index, end_index1
  real, dimension(:), allocatable :: globbuf  
  real*4,allocatable :: field3(:,:,:),field2(:,:)
  real*4,allocatable :: field3u(:,:,:),field3v(:,:,:),field3ph(:,:,:)
  character(len=4) ::  fgname
  integer :: julyr, julday
  real    :: gmt

  wrf_real=104
  end_index1=0

  call get_ijk_from_grid (  grid ,                         &
                          ids, ide, jds, jde, kds, kde,    &
                          ims, ime, jms, jme, kms, kme,    &
                          ips, ipe, jps, jpe, kps, kpe    )

!
! update wrfvar_output file with analysis variables from 3dvar
!
  if (present(out_filename)) then
     file_name = trim(out_filename)
     if (file_name == 'ana02') then
        fgname = 'fg02'
     else
        fgname = 'fg'
     endif
  else
     file_name = 'wrfvar_output'
     fgname = 'fg'
  endif

  if (rootproc) then
     ierr = copyfile(trim(fgname)//C_NULL_CHAR, trim(file_name)//C_NULL_CHAR)
     if ( ierr /= 0 ) then
        write(unit=message(1),fmt='(a)') "Failed to create "//trim(file_name)//" from "//trim(fgname)
        call da_error("da_update_firstguess.inc",99,message(1:1))
     endif

     call ext_ncd_open_for_update( trim(file_name), 0, 0, "", dh1, Status)
     if ( Status /= 0 )then
        write(unit=message(1),fmt='(a)') "Failed to open "//trim(file_name)
        call da_error("da_update_firstguess.inc",105,message(1:1))
     endif

!-------------  get date info

     call ext_ncd_get_next_time(dh1, DateStr1, Status_next_time)
     if ( var4d .or. num_fgat_time == 1 )  then ! Don't do it for FGAT
        if ( DateStr1 /= start_date )then
           ! impossible scenario
           ! start_date is set to be equal to file date in da_med_initialdata_input.inc
           write(unit=message(1),fmt='(a)') 'date info mismatch '//trim(DateStr1)//" != "//trim(start_date)
           call da_error("da_update_firstguess.inc",116,message(1:1))
        endif
     endif

     ! update analysis time info in global attributes
     ! needs to be done before the ext_ncd_write_field calls
     if ( var4d .or. num_fgat_time == 1 ) then  ! For 4dvar or 3dvar
        call get_julgmt(start_date, julyr, julday, gmt)
        CALL ext_ncd_put_dom_ti_char (dh1, 'TITLE', ' OUTPUT FROM '//trim(program_name), ierr)
        CALL ext_ncd_put_dom_ti_char (dh1, 'START_DATE', trim(start_date), ierr)
        CALL ext_ncd_put_dom_ti_char (dh1, 'SIMULATION_START_DATE', trim(start_date), ierr)
     else                                  ! For FGAT
        call get_julgmt(DateStr1//'     ', julyr, julday, gmt)
        CALL ext_ncd_put_dom_ti_char (dh1, 'TITLE', ' OUTPUT FROM '//trim(program_name), ierr)
        CALL ext_ncd_put_dom_ti_char (dh1, 'START_DATE', trim(DateStr1), ierr)
        CALL ext_ncd_put_dom_ti_char (dh1, 'SIMULATION_START_DATE', trim(DateStr1), ierr)
     endif
     CALL ext_ncd_put_dom_ti_real    (dh1, 'GMT',    gmt,    1, ierr)
     CALL ext_ncd_put_dom_ti_integer (dh1, 'JULYR',  julyr,  1, ierr)
     CALL ext_ncd_put_dom_ti_integer (dh1, 'JULDAY', julday, 1, ierr)

!-------------  get grid info
     rmse_var='T'
     call ext_ncd_get_var_info (dh1,rmse_var,ndim1,ordering,staggering, &
                               start_index,end_index1, WrfType, ierr    )
     nlon_regional=end_index1(1)
     nlat_regional=end_index1(2)
     nsig_regional=end_index1(3)

!    write(6,*)' nlon,lat,sig_regional=',nlon_regional,nlat_regional,nsig_regional
     allocate(field2(nlon_regional,nlat_regional),field3(nlon_regional,nlat_regional,nsig_regional))
     allocate(field3u(nlon_regional+1,nlat_regional,nsig_regional))
     allocate(field3v(nlon_regional,nlat_regional+1,nsig_regional))
     allocate(field3ph(nlon_regional,nlat_regional,nsig_regional+1))

  end if ! end of rootproc 

!
! update MU
! 
  if (rootproc) then
     ALLOCATE(globbuf((ide-1-ids+3)*3*(jde-1-jds+3)))    
  else
     ALLOCATE(globbuf(1))
  end if
  globbuf=0.0
  call wrf_patch_to_global_double(grid%mu_2,globbuf,1,' ','xy', &
                                  ids, ide-1,          jds, jde-1,         1, 1,  &
                                  ims, ime,            jms, jme,           1, 1,  &
                                  ips, min(ipe,ide-1), jps, min(jpe,jde-1),1, 1)
  if (rootproc) then
     do j= 1,nlat_regional 
        do i= 1,nlon_regional
           field2(i,j)=globbuf(i+(j-1)*(nlon_regional+1))     
        end do 
     end do     
     rmse_var='MU'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),       &
                             field2,WRF_REAL,0,0,0,ordering,   &
                             staggering, dimnames ,              &
                             start_index,end_index1,             & !dom
                             start_index,end_index1,             & !mem
                             start_index,end_index1,             & !pat
                             ierr                               )
  end if  

  deallocate(globbuf) 

!
! update PSFC 
!
  if (rootproc) then
     ALLOCATE(globbuf((ide-1-ids+3)*3*(jde-1-jds+3)))    ! global_mu_2(ids:ide-1,jds:jde-1) )
  else
     ALLOCATE(globbuf(1))
  end if
  globbuf=0.0
  call wrf_patch_to_global_double(grid%psfc,globbuf,1,' ','xy', &
                                  ids, ide-1,          jds, jde-1,         1, 1,  &
                                  ims, ime,            jms, jme,           1, 1,  &
                                  ips, min(ipe,ide-1), jps, min(jpe,jde-1),1, 1)
  if (rootproc) then
     do j= 1,nlat_regional
        do i= 1,nlon_regional
           field2(i,j)=globbuf(i+(j-1)*(nlon_regional+1)) 
        end do
     end do
     rmse_var='PSFC'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
                               start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),        &
                             field2,WRF_REAL,0,0,0,ordering,    &
                             staggering, dimnames ,               &
                             start_index,end_index1,              & !dom
                             start_index,end_index1,              & !mem
                             start_index,end_index1,              & !pat
                             ierr                   )
  end if

  deallocate(globbuf) 

!
! update T2 
!
  if (rootproc) then
     ALLOCATE(globbuf((ide-1-ids+3)*3*(jde-1-jds+3)))    
  else
     ALLOCATE(globbuf(1))
  end if
  globbuf=0.0
  call wrf_patch_to_global_double(grid%t2,globbuf,1,' ','xy', &
                                  ids, ide-1,          jds, jde-1,         1, 1,  &
                                  ims, ime,            jms, jme,           1, 1,  &
                                  ips, min(ipe,ide-1), jps, min(jpe,jde-1),1, 1)
  if (rootproc) then
     do j= 1,nlat_regional
        do i= 1,nlon_regional
           field2(i,j)=globbuf(i+(j-1)*(nlon_regional+1))
        end do
     end do
     rmse_var='T2'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                             field2,WRF_REAL,0,0,0,ordering,     &
                             staggering, dimnames ,                &
                             start_index,end_index1,               & !dom
                             start_index,end_index1,               & !mem
                             start_index,end_index1,               & !pat
                             ierr                                 )
  end if

  deallocate(globbuf) 
  
!
! update TH2
!
  if (rootproc) then
     ALLOCATE(globbuf((ide-1-ids+3)*3*(jde-1-jds+3)))    
  else
     ALLOCATE(globbuf(1))
  end if
  globbuf=0.0
  call wrf_patch_to_global_double(grid%th2,globbuf,1,' ','xy', &
                                  ids, ide-1,          jds, jde-1,         1, 1,  &
                                  ims, ime,            jms, jme,           1, 1,  &
                                  ips, min(ipe,ide-1), jps, min(jpe,jde-1),1, 1)
  if (rootproc) then
     do j= 1,nlat_regional
        do i= 1,nlon_regional
           field2(i,j)=globbuf(i+(j-1)*(nlon_regional+1)) 
        end do
     end do
     rmse_var='TH2'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                             field2,WRF_REAL,0,0,0,ordering,     &
                             staggering, dimnames ,                &
                             start_index,end_index1,               & !dom
                             start_index,end_index1,               & !mem
                             start_index,end_index1,               & !pat
                             ierr                                 )
  end if

  deallocate(globbuf) 

!
! update Q2
!
  if (rootproc) then
     ALLOCATE(globbuf((ide-1-ids+3)*3*(jde-1-jds+3)))   
  else
     ALLOCATE(globbuf(1))
  end if
  globbuf=0.0
  call wrf_patch_to_global_double(grid%q2,globbuf,1,' ','xy', &
                                  ids, ide-1,          jds, jde-1,         1, 1,  &
                                  ims, ime,            jms, jme,           1, 1,  &
                                  ips, min(ipe,ide-1), jps, min(jpe,jde-1),1, 1)
  if (rootproc) then
     do j= 1,nlat_regional
        do i= 1,nlon_regional
           field2(i,j)=globbuf(i+(j-1)*(nlon_regional+1))   
        end do
     end do
     rmse_var='Q2'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),        &
                             field2,WRF_REAL,0,0,0,ordering,    &
                             staggering, dimnames ,               &
                             start_index,end_index1,              & !dom
                             start_index,end_index1,              & !mem
                             start_index,end_index1,              & !pat
                             ierr                                 )
  end if

  deallocate(globbuf) 

!
! update U10
!
  if (rootproc) then
     ALLOCATE(globbuf((ide-1-ids+3)*3*(jde-1-jds+3)))    ! global_mu_2(ids:ide-1,jds:jde-1) )
  else
     ALLOCATE(globbuf(1))
  end if
  globbuf=0.0
  call wrf_patch_to_global_double(grid%u10,globbuf,1,' ','xy', &
                                  ids, ide-1,          jds, jde-1,         1, 1,  &
                                  ims, ime,            jms, jme,           1, 1,  &
                                  ips, min(ipe,ide-1), jps, min(jpe,jde-1),1, 1)
  if (rootproc) then
     do j= 1,nlat_regional
        do i= 1,nlon_regional
           field2(i,j)=globbuf(i+(j-1)*(nlon_regional+1))   
        end do
     end do
     rmse_var='U10'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                             field2,WRF_REAL,0,0,0,ordering,     &
                             staggering, dimnames ,                &
                             start_index,end_index1,               & !dom
                             start_index,end_index1,               & !mem
                             start_index,end_index1,               & !pat
                             ierr                                 )
  end if

  deallocate(globbuf) 
  
!
! update V10
!
  if (rootproc) then
     ALLOCATE(globbuf((ide-1-ids+3)*3*(jde-1-jds+3)))   
  else
     ALLOCATE(globbuf(1))
  end if
  globbuf=0.0
  call wrf_patch_to_global_double(grid%v10,globbuf,1,' ','xy', &
                                  ids, ide-1,          jds, jde-1,         1, 1,  &
                                  ims, ime,            jms, jme,           1, 1,  &
                                  ips, min(ipe,ide-1), jps, min(jpe,jde-1),1, 1)
  if (rootproc) then
     do j= 1,nlat_regional
        do i= 1,nlon_regional
           field2(i,j)=globbuf(i+(j-1)*(nlon_regional+1))   
        end do
     end do
     rmse_var='V10'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                             field2,WRF_REAL,0,0,0,ordering,     &
                             staggering, dimnames ,                &
                             start_index,end_index1,               & !dom
                             start_index,end_index1,               & !mem
                             start_index,end_index1,               & !pat
                             ierr                                 )
  end if
  
  deallocate(globbuf) 

!
!  update P
!
  if (rootproc) then
    allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
  else
    allocate( globbuf( 1 ) )
  endif
  globbuf=0.0

  CALL wrf_patch_to_global_double ( grid%p, globbuf, 1, '', 'xyz' ,             &
                                    ids, ide-1,          jds, jde-1,          kds, kde-1, &
                                    ims, ime,            jms, jme,            kms, kme,   &
                                    ips, min(ipe,ide-1), jps, min(jpe,jde-1), kps, kpe-1  )
  if (rootproc) then 
     do k= 1,nsig_regional 
        do j= 1,nlat_regional 
           do i= 1,nlon_regional
              field3(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
           end do 
        end do
     end do 
     rmse_var='P'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                             field3,WRF_REAL,0,0,0,ordering,     &
                             staggering, dimnames ,                &
                             start_index,end_index1,               & !dom
                             start_index,end_index1,               & !mem
                             start_index,end_index1,               & !pat
                             ierr                                 )
  end if  ! end of rootproc
  
  deallocate(globbuf) 

!
!  update T
!
  if (rootproc) then
    allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
  else
    allocate( globbuf( 1 ) )
  endif
  globbuf=0.0
  CALL wrf_patch_to_global_double ( grid%t_2, globbuf, 1, '', 'xyz' ,                       &
                                    ids, ide-1,          jds, jde-1,          kds, kde-1, &
                                    ims, ime,            jms, jme,            kms, kme,   &
                                    ips, min(ipe,ide-1), jps, min(jpe,jde-1), kps, kpe-1  )
  if (rootproc) then 
     do k= 1,nsig_regional 
        do j= 1,nlat_regional 
           do i= 1,nlon_regional
              field3(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
           end do 
        end do
     end do 
     rmse_var='T'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                             field3,WRF_REAL,0,0,0,ordering,       &
                             staggering, dimnames ,                &
                             start_index,end_index1,               & !dom
                             start_index,end_index1,               & !mem
                             start_index,end_index1,               & !pat
                             ierr                                 )
  end if  
  
  deallocate(globbuf) 

!
!  update QVAPOR
!
  if (rootproc) then
    allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
  else
    allocate( globbuf( 1 ) )
  endif
  globbuf=0.0
  CALL wrf_patch_to_global_double ( grid%moist(:,:,:,p_qv), globbuf, 1, '', 'xyz' ,         &
                                    ids, ide-1,          jds, jde-1,          kds, kde-1, &
                                    ims, ime,            jms, jme,            kms, kme,   &
                                    ips, min(ipe,ide-1), jps, min(jpe,jde-1), kps, kpe-1  )
  if (rootproc) then 
     do k= 1,nsig_regional 
        do j= 1,nlat_regional 
           do i= 1,nlon_regional
              field3(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
           end do 
        end do
     end do 
     rmse_var='QVAPOR'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                             field3,WRF_REAL,0,0,0,ordering,       &
                             staggering, dimnames ,                &
                             start_index,end_index1,               & !dom
                             start_index,end_index1,               & !mem
                             start_index,end_index1,               & !pat
                             ierr                                 )
  end if  ! end of rootproc

  deallocate(globbuf) 

!
!  update PH
!
  if (rootproc) then
    allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
  else
    allocate( globbuf( 1 ) )
  endif
  globbuf=0.0
  CALL wrf_patch_to_global_double ( grid%ph_2, globbuf, 1, 'Z', 'xyz' ,                   &
                                    ids, ide-1,          jds, jde-1,          kds, kde, &
                                    ims, ime,            jms, jme,            kms, kme, &
                                    ips, min(ipe,ide-1), jps, min(jpe,jde-1), kps, kpe  )
  if (rootproc) then 
     do k= 1,nsig_regional+1 
        do j= 1,nlat_regional 
           do i= 1,nlon_regional
              field3ph(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
           end do 
        end do
     end do 
     rmse_var='PH'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                             field3ph,WRF_REAL,0,0,0,ordering,     &
                             staggering, dimnames ,                &
                             start_index,end_index1,               & !dom
                             start_index,end_index1,               & !mem
                             start_index,end_index1,               & !pat
                             ierr                                 )
  end if  ! end of rootproc

  deallocate(globbuf) 

!
!  update U
!
  if (rootproc) then
    allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
  else
    allocate( globbuf( 1 ) )
  endif
  globbuf=0.0
  CALL wrf_patch_to_global_double ( grid%u_2, globbuf, 1, 'X', 'xyz' ,                    &
                                    ids, ide,          jds, jde-1,          kds, kde-1, &
                                    ims, ime,          jms, jme,            kms, kme,   &
                                    ips, min(ipe,ide), jps, min(jpe,jde-1), kps, kpe-1  )
  if (rootproc) then 
     do k= 1,nsig_regional 
        do j= 1,nlat_regional 
           do i= 1,nlon_regional+1
              field3u(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
           end do 
        end do
     end do 
     rmse_var='U'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                             field3u,WRF_REAL,0,0,0,ordering,      &
                             staggering, dimnames ,                &
                             start_index,end_index1,               & !dom
                             start_index,end_index1,               & !mem
                             start_index,end_index1,               & !pat
                             ierr                                 )
  end if  ! end of rootproc

  deallocate(globbuf) 

!
!  update V
!
  if (rootproc) then
    allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
  else
    allocate( globbuf( 1 ) )
  endif
  globbuf=0.0
  CALL wrf_patch_to_global_double ( grid%v_2, globbuf, 1, 'Y', 'xyz' ,                    &
                                    ids, ide-1,          jds, jde,          kds, kde-1, &
                                    ims, ime,            jms, jme,          kms, kme,   &
                                    ips, min(ipe,ide-1), jps, min(jpe,jde), kps, kpe-1  )
  if (rootproc) then 
     do k= 1,nsig_regional 
        do j= 1,nlat_regional+1 
           do i= 1,nlon_regional
              field3v(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
           end do 
        end do
     end do 
     rmse_var='V'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                             field3v,WRF_REAL,0,0,0,ordering,      &
                             staggering, dimnames ,                &
                             start_index,end_index1,               & !dom
                             start_index,end_index1,               & !mem
                             start_index,end_index1,               & !pat
          ierr                                 )
  end if  ! end of rootproc

  deallocate(globbuf) 

!
!  update W
!
  if (rootproc) then
    allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
  else
    allocate( globbuf( 1 ) )
  endif
  globbuf=0.0
  CALL wrf_patch_to_global_double ( grid%w_2, globbuf, 1, 'Z', 'xyz' ,                    &
                                    ids, ide-1,          jds, jde-1,          kds, kde, &
                                    ims, ime,            jms, jme,            kms, kme, &
                                    ips, min(ipe,ide-1), jps, min(jpe,jde-1), kps, kpe  )
  if (rootproc) then 
     do k= 1,nsig_regional+1 
        do j= 1,nlat_regional 
           do i= 1,nlon_regional
              field3ph(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
           end do 
        end do
     end do 
     rmse_var='W'
     call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
          start_index,end_index1, WrfType, ierr    )
!    write(6,*)' rmse_var=',trim(rmse_var)
!    write(6,*)' ordering=',ordering
!    write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!    write(6,*)' ndim1=',ndim1
!    write(6,*)' staggering=',staggering
!    write(6,*)' start_index=',start_index
!    write(6,*)' end_index1=',end_index1
     call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                             field3ph,WRF_REAL,0,0,0,ordering,     &
                             staggering, dimnames ,                &
                             start_index,end_index1,               & !dom
                             start_index,end_index1,               & !mem
                             start_index,end_index1,               & !pat
                             ierr                                  )
  end if    

  deallocate(globbuf) 

!-------------Update QCLOUD, QRAIN, QICE, QSNOW & QGROUP 
if ( (use_radarobs .and. use_radar_rhv) .or. (use_radarobs .and. use_radar_rqv) .or. (use_rad .and. crtm_cloud) ) then
   if (size(grid%moist,dim=4) >= 4) then         ! update QCLOUD & QRAIN
!
!  update QCLOUD
!
      if (rootproc) then
         allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
      else
         allocate( globbuf( 1 ) )
      endif
      globbuf=0.0
      CALL wrf_patch_to_global_double ( grid%moist(:,:,:,p_qc), globbuf, 1, '', 'xyz' ,     &
                                    ids, ide-1,          jds, jde-1,          kds, kde-1, &
                                    ims, ime,            jms, jme,            kms, kme,   &
                                    ips, min(ipe,ide-1), jps, min(jpe,jde-1), kps, kpe-1  )
      if (rootproc) then 
         do k= 1,nsig_regional 
            do j= 1,nlat_regional 
               do i= 1,nlon_regional
                  field3(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
               end do 
            end do
         end do 
         rmse_var='QCLOUD'
         call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
              start_index,end_index1, WrfType, ierr    )
!        write(6,*)' rmse_var=',trim(rmse_var)
!        write(6,*)' ordering=',ordering
!        write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!        write(6,*)' ndim1=',ndim1
!        write(6,*)' staggering=',staggering
!        write(6,*)' start_index=',start_index
!        write(6,*)' end_index1=',end_index1
         call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                                 field3,WRF_REAL,0,0,0,ordering,       &
                                 staggering, dimnames ,                &
                                 start_index,end_index1,               & !dom
                                 start_index,end_index1,               & !mem
                                 start_index,end_index1,               & !pat
                                 ierr                                  )
      end if  ! end of rootproc

      deallocate(globbuf) 

!
!  update QRAIN
!
      if (rootproc) then
         allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
      else
         allocate( globbuf( 1 ) )
      endif
      globbuf=0.0
      CALL wrf_patch_to_global_double ( grid%moist(:,:,:,p_qr), globbuf, 1, '', 'xyz' ,     &
                                    ids, ide-1,          jds, jde-1,          kds, kde-1, &
                                    ims, ime,            jms, jme,            kms, kme,   &
                                    ips, min(ipe,ide-1), jps, min(jpe,jde-1), kps, kpe-1  )
      if (rootproc) then 
         do k= 1,nsig_regional 
            do j= 1,nlat_regional 
               do i= 1,nlon_regional
                  field3(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
               end do 
            end do
         end do 
         rmse_var='QRAIN'
         call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
              start_index,end_index1, WrfType, ierr    )
!        write(6,*)' rmse_var=',trim(rmse_var)
!        write(6,*)' ordering=',ordering
!        write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!        write(6,*)' ndim1=',ndim1
!        write(6,*)' staggering=',staggering
!        write(6,*)' start_index=',start_index
!        write(6,*)' end_index1=',end_index1
         call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                                 field3,WRF_REAL,0,0,0,ordering,       &
                                 staggering, dimnames ,                &
                                 start_index,end_index1,               & !dom
                                 start_index,end_index1,               & !mem
                                 start_index,end_index1,               & !pat
                                 ierr                                  )
      end if  ! end of rootproc

      deallocate(globbuf) 

   end if    ! end of update QCLOUD & QRAIN

   if (size(grid%moist,dim=4) >= 6) then     ! update QICE & QSNOW
!
!  update QICE
!
      if (rootproc) then
         allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
      else
         allocate( globbuf( 1 ) )
      endif
      globbuf=0.0
      CALL wrf_patch_to_global_double ( grid%moist(:,:,:,p_qi), globbuf, 1, '', 'xyz' ,     &
                                    ids, ide-1,          jds, jde-1,          kds, kde-1, &
                                    ims, ime,            jms, jme,            kms, kme,   &
                                    ips, min(ipe,ide-1), jps, min(jpe,jde-1), kps, kpe-1  )
      if (rootproc) then 
         do k= 1,nsig_regional 
            do j= 1,nlat_regional 
               do i= 1,nlon_regional
                  field3(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
               end do 
            end do
         end do 
         rmse_var='QICE'
         call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
              start_index,end_index1, WrfType, ierr    )
!        write(6,*)' rmse_var=',trim(rmse_var)
!        write(6,*)' ordering=',ordering
!        write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!        write(6,*)' ndim1=',ndim1
!        write(6,*)' staggering=',staggering
!        write(6,*)' start_index=',start_index
!        write(6,*)' end_index1=',end_index1
         call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                                 field3,WRF_REAL,0,0,0,ordering,       &
                                 staggering, dimnames ,                &
                                 start_index,end_index1,               & !dom
                                 start_index,end_index1,               & !mem
                                 start_index,end_index1,               & !pat
                                 ierr                                 )
      end if  ! end of rootproc

      deallocate(globbuf) 

!
!  update QSNOW
!
      if (rootproc) then
         allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
      else
         allocate( globbuf( 1 ) )
      endif
      globbuf=0.0
      CALL wrf_patch_to_global_double ( grid%moist(:,:,:,p_qs), globbuf, 1, '', 'xyz' ,     &
                                    ids, ide-1,          jds, jde-1,          kds, kde-1, &
                                    ims, ime,            jms, jme,            kms, kme,   &
                                    ips, min(ipe,ide-1), jps, min(jpe,jde-1), kps, kpe-1  )
      if (rootproc) then 
         do k= 1,nsig_regional 
            do j= 1,nlat_regional 
               do i= 1,nlon_regional
                  field3(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
               end do 
            end do
         end do 
         rmse_var='QSNOW'
         call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
              start_index,end_index1, WrfType, ierr    )
!        write(6,*)' rmse_var=',trim(rmse_var)
!        write(6,*)' ordering=',ordering
!        write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!        write(6,*)' ndim1=',ndim1
!        write(6,*)' staggering=',staggering
!        write(6,*)' start_index=',start_index
!        write(6,*)' end_index1=',end_index1
         call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                                 field3,WRF_REAL,0,0,0,ordering,       &
                                 staggering, dimnames ,                &
                                 start_index,end_index1,               & !dom
                                 start_index,end_index1,               & !mem
                                 start_index,end_index1,               & !pat
                                 ierr                                 )
      end if  ! end of rootproc

      deallocate(globbuf) 

   end if    ! end of update QICE & QSNOW 

   if (size(grid%moist,dim=4) >= 7) then      ! update QGRAUP
!
!  update QGRAUP
!
      if (rootproc) then
         allocate( globbuf((ide-1-ids+3)*(jde-1-jds+3)*(kde-1-kds+3)))
      else
         allocate( globbuf( 1 ) )
      endif
      globbuf=0.0
      CALL wrf_patch_to_global_double ( grid%moist(:,:,:,p_qg), globbuf, 1, '', 'xyz' ,     &
                                    ids, ide-1,          jds, jde-1,          kds, kde-1, &
                                    ims, ime,            jms, jme,            kms, kme,   &
                                    ips, min(ipe,ide-1), jps, min(jpe,jde-1), kps, kpe-1  )
      if (rootproc) then 
         do k= 1,nsig_regional 
            do j= 1,nlat_regional 
               do i= 1,nlon_regional
                  field3(i,j,k)=globbuf(i+(j-1)*(nlon_regional+1)+(k-1)*(nlon_regional+1)*(nlat_regional+1))
               end do 
            end do
         end do 
         rmse_var='QGRAUP'
         call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
              start_index,end_index1, WrfType, ierr    )
!        write(6,*)' rmse_var=',trim(rmse_var)
!        write(6,*)' ordering=',ordering
!        write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
!        write(6,*)' ndim1=',ndim1
!        write(6,*)' staggering=',staggering
!        write(6,*)' start_index=',start_index
!        write(6,*)' end_index1=',end_index1
         call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),         &
                                 field3,WRF_REAL,0,0,0,ordering,       &
                                 staggering, dimnames ,                &
                                 start_index,end_index1,               & !dom
                                 start_index,end_index1,               & !mem
                                 start_index,end_index1,               & !pat
                                 ierr                                 )
      end if  ! end of rootproc

      deallocate(globbuf) 

   end if     ! end of update QGRAUP
end if        ! end of radar or radiance
!-------------End of update QCLOUD, QRAIN, QICE, QSNOW & QGROUP

  if (rootproc) then
     deallocate(field2,field3,field3u,field3v,field3ph)
     call ext_ncd_ioclose(dh1, Status)
  end if   ! end of rootproc
  
end subroutine da_update_firstguess 

end module da_wrfvar_io
