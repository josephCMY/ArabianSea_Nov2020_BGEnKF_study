!--------------------------------------------------------------------------------------------------------------
! Subroutine for BT-based inflation - by Masashi Minamide 2017.
! ----------------------------------------
subroutine cal_abei(wrf_file,ix,jx,kx,kxs,ni,nj,nk,nv,nm,g_comm,s_comm,gid,sid,iid,jid,xm,x,ind,proj,times)
use constants
use namelist_define
use obs_define
use mpi_module
use netcdf
use radar
use wrf_field_indices_define
use boundary_variables
implicit none
! variable description:
! numbers_en = number of ensemble members (n)
! nicpu, njcpu = number of cpus used in x and y direction. nicpu*njcpu = number of cpus used for a domain
! sid = subdomain id for a cpu. iid, jid= cpu id in decomposed subdomain (coordinate in x and y direction)
! gid = cpugroup id. A cpugroup is cpus with the same gid, but different sid's.
!                    the ensemble members are distributed to each cpugroups in the same way.
! nmcpu = number of cpugroups
! g_comm = MPI_COMM splitted with same sid (all cpus in the same subdomain)
! s_comm = MPI_COMM splitted with same gid (all cpus in the same cpugroup, i.e. covering all subdomains)
! ind : randomly permutated indices
integer, intent(in) :: sid,gid,iid,jid
integer, intent(in) :: ix,jx,kx,kxs,ni,nj,nk,nv,nm, g_comm,s_comm
integer, dimension(obs%num), intent(in) :: ind
character (len=10), intent(in) :: wrf_file 
type (proj_info)    :: proj
character (len=10)  :: obstype
character (len=80), intent(in) :: times
! fac = 1/(n-1) used for averaging in variance calculation
! corr_coef = localization factor
! ngx, ngz = roi in horizontall and vertical
! var,cov = error variance,covariance of something
! y_hxm = y-hxm (the innovation vector, mean), hxa is the perturbation (of ensemble members).
real      :: fac,d,alpha,var,cov,y_hxm_inf,corr_coef,d_ogn
real      :: var_a,var_b,m_d2,m_var_a,m_var_b
real,dimension(numbers_en) :: hxa
integer   :: ngx, ngz
integer   :: i, j, k, m, n, iob, iiob, ie, iunit,ounit,ii, jj, kk, is, it, ig, iv, i1,j1, itot
integer   :: ist,ied,jst,jed,kst,ked, istart,iend,jstart,jend, uist,uied,ujst,ujed
integer   :: sist,sied,sjst,sjed, sistart,siend,sjstart,sjend,sis,sie,sjs,sje
real      :: gaussdev, error, xb
integer, dimension(8)  :: values
character (len=10) :: filename, date, time, zone, varname
character (len=20) :: format1
double precision :: timer,t0,t1,t2,t3,t4,t5
integer :: grid_id, fid, rcode, num_update_var, assimilated_obs_num, update_flag
real                    :: p_top
real, dimension(ix, jx) :: xlong, xlat, xland, lu_index
real, dimension(kx+1)   :: znw
real, dimension(kx)     :: znu
real, dimension(3)      :: center_xb
! state vectors (x): xm=ensemble mean; x=perturbation from ensemble mean
! _f=first guess; _t=truth
! ya=Hxa,Hxm: the state vector translated to observation space. ya(number of obs, number of members + mean)
! km        : kalman gain in updating x
real, dimension(ni,nj,nk,nv,nm), intent(inout) :: x
real, dimension(ni,nj,nk,nv), intent(inout)    :: xm 
!real, dimension(ni,nj,nk,nv,nm)                :: xf
!real, dimension(ni,nj,nk,nv)                   :: x_t, std_x,std_xf
!real, dimension(3,3,nk,nv,nm,nicpu*njcpu) :: xobsend, xob
real, dimension(3,3,nk,nv) :: tmp, tmpsend
real, dimension(obs%num) :: ya, yasend, yf
integer, dimension(obs%num) :: kick_flag
real, allocatable, dimension(:,:,:)   :: km
! for satellite radiance
integer :: iob_radmin,iob_radmax
real, dimension(obs%num) :: yasend_nc, ya_tb, ya_nc, ya_ca
real, dimension(obs%num,1) :: yasend_tb
real, dimension(ix+1, jx+1,kx+1) :: inflate_bt
!real, dimension(ni,nj,nk)     :: xq_n,xq_p
!real, dimension(ni,nj,nk,nm)  :: xq_nsend,xq_psend
!real, dimension(ix+1,jx+1,kx+1)  :: inflate_bt
real  :: inflate_tmp
! inflation parameters
integer :: irec, iost
real    :: inflate_update, error_update, infvar_update
real, parameter :: min_inflate = 1.
real, parameter :: max_inflate = 1.4
logical, parameter :: use_adaptive_bt_inf = .false.

! Variables for re-decomposition of domain
! sm prefix stands for small
integer :: sm_iid, sm_jid, sm_nicpu, sm_njcpu, sm_sid, sm_gid
integer :: sm_ni, sm_nj, sm_flag
integer :: sm_ist, sm_ied, sm_jst, sm_jed
integer :: sm_ist2, sm_ied2, sm_jst2, sm_jed2,sm_sid2, sm_iid2, sm_jid2
real, allocatable, dimension(:,:,:,:,:)   :: sm_xm
real, allocatable, dimension(:,:,:)   :: dat3d

integer :: sm_comm ! Communicator for processes doing crtm

read(wrf_file(6:10),'(i5)')iunit
if ( my_proc_id==0 ) then
   write(*, *)'   '
   write(*, *)'---------------------------------------------------'
   write(*, *)'.... Begin BT-based inflation ....'
   write(*, *)'   '
endif

! Get domain info
call get_variable2d( wrf_file, 'XLONG     ', ix, jx, 1, xlong )
call get_variable2d( wrf_file, 'XLAT      ', ix, jx, 1, xlat )
call get_variable2d( wrf_file, 'LANDMASK  ', ix, jx, 1, xland )
call get_variable2d( wrf_file, 'LU_INDEX  ', ix, jx, 1, lu_index )
call get_variable1d( wrf_file, 'ZNW       ', kx+1, 1, znw )
call get_variable1d( wrf_file, 'ZNU       ', kx, 1, znu )
call get_variable0d( wrf_file, 'P_TOP     ', 1, p_top )
call open_file(wrf_file, nf_nowrite, fid)
rcode = nf_get_att_int(fid, nf_global, 'GRID_ID', grid_id)
num_update_var = 0
do m = 1, 20
   if ( len_trim(updatevar(m))>=1 ) num_update_var=num_update_var+1
enddo

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

assimilated_obs_num = 0
if( raw%radiance%num > 0 ) then !!=======================================================

! I. calculate ya=Hxa,Hxm
! 1. find obs location in z direction, return as obs%position(:,3)
!    for simulated obs, calculate obs value.
if(my_proc_id==0) write(*,*) 'Calculating Hx...'
timer = MPI_Wtime()


! 2. Parallellized computation of Tb
!===============================================================================
! Note: xb_to_radiance operates uses x or xm to generate Tb.
! Thus, for efficient parallelization, we are re-decomposing the slabs
! temporarily for running crtm on xm.
! 
! This re-decomposition only applies to running xb_to_radiance.
! The inflation procedure still uses the nicpu*njcpu decomposition seen in
! enkf.f
!==============================================================================
ya=0.
yasend=0.
yasend_tb=0.

! Generating parameters for the re-decomposed xm field
! Note that the method of decomposition I am using doesnt use all
! the available threads
if (nmcpu > 1) then
  sm_nicpu = nicpu * int( sqrt( real( nmcpu ) ) )
  sm_njcpu = njcpu * int( sqrt( real( nmcpu ) ) )
  sm_sid = mod(my_proc_id, sm_nicpu*sm_njcpu)
  sm_iid = mod(sm_sid,sm_nicpu)
  sm_jid = int(sm_sid/sm_nicpu)
  sm_ni = int( (ix+1)/sm_nicpu ) + 1
  sm_nj = int( (jx+1)/sm_njcpu ) + 1
  sm_flag = int( my_proc_id/(sm_nicpu*sm_njcpu) )
else
  sm_nicpu = nicpu
  sm_njcpu = njcpu
  sm_sid = sid
  sm_iid = iid
  sm_jid = jid
  sm_ni  = ni
  sm_nj  = nj
  sm_flag = 0
endif


! Generate communicator for processes doing CRTM
call MPI_Comm_split( comm, sm_flag, sm_sid, sm_comm, ierr )

write( filename, '(a5,i5.5)') wrf_file(1:5), iunit+numbers_en


! Only threads assigned to do CRTM will operate from here on.
sm_procs: if (sm_flag == 0) then

  ! Allocate memory for re-decomposed mean
  allocate( sm_xm(sm_ni, sm_nj, nk, nv, 1) )

  ! 2.1 Loading the ensemble mean field
  load_xm_loop: do m = 1, nv

    call wrf_var_dimension(filename,enkfvar(m),ix,jx,kx, kxs, ii,jj,kk)
    ! All processes determine their sub-domains
    sm_ist=sm_iid*sm_ni+1
    sm_ied=(sm_iid+1)*sm_ni
    sm_jst=sm_jid*sm_nj+1
    sm_jed=(sm_jid+1)*sm_nj
    if(sm_iid==(sm_nicpu-1)) sm_ied=ii
    if(sm_jid==(sm_njcpu-1)) sm_jed=jj
 
    ! Root process reads in mean field and sends out to all
    ! sm procs
    root_sm_proc: if (sm_sid == 0) then

       ! Read-in procedure
       allocate(dat3d(ii,jj,kk))

       if(kk>1) then
         call get_variable3d(filename,enkfvar(m),ii,jj,kk,1,dat3d)
       else if (kk==1) then
         call get_variable2d(filename,enkfvar(m),ii,jj,1,dat3d)
       endif
  
       ! Pass decomposed info to self (sid = 0 procedure)
       sm_xm(1:sm_ied-sm_ist+1,1:sm_jed-sm_jst+1,1:kk,m,1) = &
         dat3d(sm_ist:sm_ied,sm_jst:sm_jed,1:kk)

       ! Sending info from sid = 0 to sid > 0
       send_slab_loop: do sm_sid2 = 1, sm_nicpu*sm_njcpu-1
          
         ! Send information from root sid
         ! Determine sid>0's subdomain
         sm_iid2=mod(sm_sid2,sm_nicpu)
         sm_jid2=int(sm_sid2/sm_nicpu)
         sm_ist2=sm_iid2*sm_ni+1
         sm_ied2=(sm_iid2+1)*sm_ni
         sm_jst2=sm_jid2*sm_nj+1
         sm_jed2=(sm_jid2+1)*sm_nj
         if(sm_iid2==(sm_nicpu-1)) sm_ied2=ii
         if(sm_jid2==(sm_njcpu-1)) sm_jed2=jj

         call MPI_Send( dat3d(sm_ist2:sm_ied2,sm_jst2:sm_jed2,1:kk), &
                        (sm_ied2-sm_ist2+1)*(sm_jed2-sm_jst2+1)*kk, &
                        MPI_REAL, sm_sid2, sm_sid2, sm_comm, ierr )
       enddo send_slab_loop
 
    endif root_sm_proc

    ! Members read info from root sid
    sm_recv_proc: if (sm_sid>0) then
      ! Receive data
      call MPI_Recv( &
        sm_xm(1:sm_ied-sm_ist+1,1:sm_jed-sm_jst+1,1:kk,m,1), &
        (sm_ied-sm_ist+1)*(sm_jed-sm_jst+1)*kk, &
        MPI_REAL, 0, sm_sid, sm_comm, status, ierr )
    endif sm_recv_proc

    if (sm_sid==0) deallocate(dat3d)

    call MPI_Barrier( sm_comm, ierr)
  
  enddo load_xm_loop

  ! 2.2: Run CRTM on the ensemble mean
  ! Situation where clouds are present
  call xb_to_radiance( sm_xm, proj, ix, jx, kx, sm_ni, sm_nj, nk, &
                       nv, 1, xlong, xlat, xland, iob_radmin, &
                       iob_radmax, yasend_tb, sm_iid, sm_jid, &
                       sm_nicpu, sm_njcpu, .true.)
  yasend(iob_radmin:iob_radmax) = yasend_tb(iob_radmin:iob_radmax,1)

  call MPI_Barrier( sm_comm, ierr)
  ! Situation where clouds are removed
  sm_xm(:,:,:,ind_qcloud,:) = 0.
  sm_xm(:,:,:,ind_qrain ,:) = 0.
  sm_xm(:,:,:,ind_qice  ,:) = 0.
  sm_xm(:,:,:,ind_qsnow ,:) = 0.
  sm_xm(:,:,:,ind_qgraup,:) = 0.
  call xb_to_radiance( sm_xm, proj, ix, jx, kx, sm_ni, sm_nj, nk, &
                       nv, 1, xlong, xlat, xland, iob_radmin, &
                       iob_radmax, yasend_tb, sm_iid, sm_jid, &
                       sm_nicpu, sm_njcpu, .false.)
  yasend_nc(iob_radmin:iob_radmax) = yasend_tb(iob_radmin:iob_radmax,1)

  ! Deallocate after running CRTM
  deallocate( sm_xm )
endif sm_procs

! Gather ya and ya_nc together
call MPI_Allreduce( yasend, ya, obs%num, MPI_REAL, MPI_SUM, comm, ierr)
call MPI_Allreduce( yasend_nc, ya_nc, obs%num, MPI_REAL, MPI_SUM, comm, ierr)

!if (my_proc_id==0) write(*,*)'ya',ya
!if (my_proc_id==0) write(*,*)ya_nc

!! Symmetric cloud parameter
!ya_ca = (abs(obs%dat-ya_nc)+abs(ya(:,numbers_en+1)-ya_nc))/2.
!! Asymmetric cloud parameter
ya_ca = abs(obs%dat-ya_nc) - abs(ya-ya_nc)

!----------------------------------------------------------------------------------------
! II. calculate inflation factor with radiance obs
if(my_proc_id==0) write(*,*) 'Seqencing obs...'

timer=MPI_Wtime()
t1=0.; t2=0.; t3=0.; t4=0.;

inflate_bt = 1.
obs_assimilate_cycle : do it = 1,obs%num
! 0. get basic information of the obs being assimilated
   iob=ind(it)
   obstype = obs%type(iob)

   if (obstype=='Radiance  ') then
      call get_inflate_facotr_bt(ya_ca(iob),inflate_tmp,obs%sat(iob),obs%ch(iob))
      y_hxm_inf = inflate_tmp - inflate_bt(int(obs%position(iob,1)),int(obs%position(iob,2)),int(obs%position(iob,3)))
      if ( my_proc_id==0 ) write(*,'(a,i6,a,f6.2,a,f6.2,a,f6.2,a,f6.2)') &
'No.',iob,' '//obstype//'(inf) =',inflate_tmp,', ca=', ya_ca(iob),', inf_b=',inflate_tmp-y_hxm_inf,', delta_inf=',y_hxm_inf
      ngx = 10
      ngz = 999 !obs%roi(iob,2)
      ! start and end indices of the update zone of the obs
      ist = max( 1,    int(obs%position(iob,1))-ngx )
      ied = min( ix+1, int(obs%position(iob,1))+ngx )
      jst = max( 1,    int(obs%position(iob,2))-ngx )
      jed = min( jx+1, int(obs%position(iob,2))+ngx )
      kst = max( 1,    int(obs%position(iob,3))-ngz )
      ked = min( kx+1, int(obs%position(iob,3))+ngz )
      call wrf_var_dimension ( wrf_file, varname, ix, jx, kx, kxs, ii, jj, kk )

      !write(*,*)'RGN'
      !write(*,*)varname,ix,jx,kx
      if ( kk == 1 ) then
        kst  = 1
        ked = 1
      endif
      
      ! 1. calculate corr_coef=localization factor, 1 at obs position, reduce to 0 at roi
      !    2-D variables are treated as at surface k=1
      allocate ( km(ied-ist+1, jed-jst+1, ked-kst+1) )

      km=0.
      do k = kst,ked
      do j = jst,jed
      do i = ist,ied
        call corr(real(i-obs%position(iob,1)),real(j-obs%position(iob,2)),0.,ngx,ngz,corr_coef)
        km(i-ist+1,j-jst+1,k-kst+1) = corr_coef
      enddo
      enddo
      enddo

      ! 2. Update inflate_bt
      inflate_bt(ist:ied,jst:jed,kst:ked) = inflate_bt(ist:ied,jst:jed,kst:ked) + km * y_hxm_inf

      deallocate(km)

   endif
end do obs_assimilate_cycle

if ( my_proc_id==0 ) write(*,*)'inf_mean= ',sum(inflate_bt(:,:,1))/float((ix+1)*(jx+1)),'inf_max= ',maxval(inflate_bt)

! QC on inflation factor
where(inflate_bt < min_inflate) inflate_bt = min_inflate
where(inflate_bt > max_inflate) inflate_bt = max_inflate

if ( my_proc_id==0 ) write(*,*)'final inf_mean= ', sum(inflate_bt(:,:,1))/float((ix+1)*(jx+1))

if(my_proc_id==10) then
  open(10,file='inflation_factor.txt',form='formatted')
  do j = 1,jx
    do i = 1,ix
      write(10,*) i,j,inflate_bt(i,j,1)
    enddo
  enddo
  close (10)
endif

if(my_proc_id==11) then
  open(11,file='cloud_parameter.txt',form='formatted')
  do iob = 1,obs%num
    obstype = obs%type(iob)
    if (obstype=='Radiance  ') then
      write(11,*) int(obs%position(iob,1)),int(obs%position(iob,2)),ya_ca(iob)
    endif
  enddo
  close (11)
endif

!---------------------------------------------------------------------------------
! II. inflate ensemble members

do n = 1, nm
   ie = (n-1)*nmcpu+gid+1
   inflate_x_loop: do m=1,nv
     ! Check if variable is to be inflated
     varname=enkfvar(m)
     update_flag = 0
     do iv = 1, num_update_var
       if ( varname .eq. updatevar(iv) ) update_flag = 1
     enddo

     ! Inflate slab
     if ( ie<=numbers_en+1 .and. update_flag == 1 ) x(:,:,:,m,n)=inflate_bt(istart:iend,jstart:jend,:)*(x(:,:,:,m,n)-xm(:,:,:,m)) + xm(:,:,:,m)

     ! Inflate slab boundaries
     if ( ie <= numbers_en .and. update_flag ==1) then

       if ( istart > 1 ) &
         x_we_bdy(1,:,:,m,n) = inflate_bt(istart-1, jstart:jend, :)*( x_we_bdy(1,:,:,m,n) - xm_we_bdy(1,:,:,m) ) +  xm_we_bdy(1,:,:,m)

       if ( iend+1 <= ix+1 ) &
         x_we_bdy(2,:,:,m,n) = inflate_bt(  iend+1, jstart:jend, :)*( x_we_bdy(2,:,:,m,n) - xm_we_bdy(2,:,:,m) ) +  xm_we_bdy(2,:,:,m)

       if ( jstart > 1 ) &
         x_sn_bdy(1,:,:,m,n) = inflate_bt(istart:iend, jstart-1, :)*( x_sn_bdy(1,:,:,m,n) - xm_sn_bdy(1,:,:,m) ) +  xm_sn_bdy(1,:,:,m)

       if ( jend+1 <= jx+1 ) &
         x_sn_bdy(2,:,:,m,n) = inflate_bt(istart:iend,   jend+1, :)*( x_sn_bdy(2,:,:,m,n) - xm_sn_bdy(2,:,:,m) ) +  xm_sn_bdy(2,:,:,m)

       if ( istart > 1 .and. jstart > 1 ) &
         x_crnrs(1,1,:,m,n) = inflate_bt(istart-1,jstart-1,:)*(x_crnrs(1,1,:,m,n)-xm_crnrs(1,1,:,m)) + xm_crnrs(1,1,:,m)

       if ( istart > 1 .and. jend+1 <= jx+1 ) &
         x_crnrs(1,2,:,m,n) = inflate_bt(istart-1,jend+1,:)  *(x_crnrs(1,2,:,m,n)-xm_crnrs(1,2,:,m)) + xm_crnrs(1,2,:,m)

       if ( iend+1 <= ix+1 .and. jstart > 1 ) &
         x_crnrs(2,1,:,m,n) = inflate_bt(iend+1,jstart-1,:)  *(x_crnrs(2,1,:,m,n)-xm_crnrs(2,1,:,m)) + xm_crnrs(2,1,:,m)

       if ( iend+1 <= ix+1 .and. jend+1 <= jx+1 ) &
         x_crnrs(2,2,:,m,n) = inflate_bt(iend+1,jend+1,:)    *(x_crnrs(2,2,:,m,n)-xm_crnrs(2,2,:,m)) + xm_crnrs(2,2,:,m)

     endif

   enddo inflate_x_loop
enddo


else !!=======================================================================

   if ( my_proc_id==0 ) write(*,*)'No radaince data'

endif !!======================================================================

call MPI_Comm_free(sm_comm,ierr)


end subroutine cal_abei
!==============================================================================
subroutine get_inflate_facotr_bt(ca,inf,satid,ch)
! This is routine to give empirical inflation factor based on cloud parameter
implicit none
real, intent(in)    :: ca
character(*), intent(in) :: satid
integer, intent(in)    :: ch
real, intent(inout) :: inf
! input parameters computed offline
!real, parameter     :: a = 0.015
real :: a

if ( trim(adjustl(satid)) == 'abi_gr' &
     .or. trim(adjustl(satid)) == 'ahi_h8' &
     .or. trim(adjustl(satid)) == 'abi_g16' ) then
   if (ch == 8) then
      a = 0.015
   else if (ch == 9) then
      a = 0.012
   else if (ch == 10) then
      a = 0.009
   else if (ch == 14) then
      a = 0.009
   else
      write(*,*) 'ABEI ERROR!!! Please prepare ABEI parameter for channel ',ch
      stop
   endif
else if (trim(adjustl(satid)) == 'imgr_g13') then
   if (ch == 3) then
      a = 0.014
   else if (ch == 4) then
      a = 0.010
   else
      write(*,*) 'ABEI ERROR!!! Please prepare ABEI parameter for channel ',ch
      stop
   endif
else if ( trim(adjustl(satid)) == 'mviriNOM_m07'     &
          .or. trim(adjustl(satid)) == 'mviri_m07') then
  if ( ch == 3 ) then
    a = 0.015       ! Roughly same central wavelength as ch 8 of ABI GOES-R
  else if (ch == 2) then
    a = 0.009       ! Roughly same central wavelength as ch 14 of ABI GOES-R
  else
    write(*,*) 'ABEI ERROR!!! Please prepare ABEI parameter for channel ',ch
    stop
  endif
! Inflation factor for seviri_m08
else if (trim(adjustl(satid)) == 'seviri_m08') then
  if ( ch == 5 ) then
    a = 0.015       ! Roughly same central wavelength as ch 8 of ABI GOES-R

  else
    write(*,*) 'ABEI ERROR!!! Please prepare ABEI parameter for channel ',ch
    stop
  endif

else
      write(*,*) 'ABEI ERROR!!! Please prepare ABEI parameter for ',trim(adjustl(satid))
      stop
endif 

if (ca <= 0.) then
   inf = 1.
else
   inf = a*ca + 1.
endif

!if (inf < min_inflate) inf = min_inflate
!if (inf > max_inflate) inf = max_inflate

end subroutine get_inflate_facotr_bt

