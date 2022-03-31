!--------------------------------------------------------------------------------------------------------------
! Subroutine for enkf update
!
! Algorithm description:
! ---------------------
!     The most time consuming calculation is in H. So instead of calculating HBH' each time in the update step,
! we calculate yf=Hxf beforehand, and update both xf and yf during the update. This enables us to use as many
! cpus we want during the calculation of yf. And the update step is much faster.
!
! Original EnKF equations:
! -----------------------
!     Hereafter, we use x and y for shorthand of xf and yf=Hxf, and xa = x-<x> = perturbation, xm = <x> =
! mean for ensemble members. <x> denotes the expected value of x. x' denotes the transpose of x. 
!     Dimensions of variables: x(n) is the state vector, yo(p) is the obs, H(pn) maps x->y, B(nn) is covariance
! of x, R(pp) is covariance of yo.
!     n = state vector length = total number of model grid points
!     p = number of obs
!     m = number of ensemble members
!     The update equations:
!         xm = xm + BH'(R+HBH')^-1 (yom - Hxm)
!         xa = xa + BH'(R+HBH')^-1 (yoa - Hxa) 
!         B  = cov(xa,xa) = <xa xa'> approx. by sum over m (xa xa')/(m-1)
!         R  = cov(yoa,yoa)
!     For x, on the left hand side is the analysis field (posterior, x_a), and on the right hand side is 
! evaluated by the first guess fields (prior, xf). For xf, B=P_f, for x_a, B=P_a.
!
! Serial EnKF implementation:
! --------------------------
!     If R is diagonal, with components (r1,r2,...,rp), we can update x with 1 obs at a time (sequentially).
! So, instead of xf->x_a, P_f->P_a, we have xf=x_1->x_2->...->x_i->...->x_p+1=x_a, and P_f->...->B_i->...->P_a.
! The update equations:
!         for i=1,2,...,p
!         xm_i+1 = xm_i + B_i hi'(ri + hi B_i hi')^-1 (yom_i - hi xm_i)
!         xa_i+1 = xa_i + B_i hi'(ri + hi B_i hi')^-1 (yoa_i - hi xa_i)
!         B_i+1  = <xa_i+1 xa_i+1'>
! Proof of equivalence:
!     First, we simplify the original EnKF equation, given: (1) yo=H xo, (2) A+ABA=A(I+BA)=(I+AB)A =>
! A (I+BA)^-1 = (I+AB)^-1 A, (3) R is diagonal, so there is S satisfying R=S'S, (si^2=ri).
!         x=x+BH'(R+HBH')^-1 H(xo-x)
!          =x+BH'[S'(I+S'^-1 HBH' S^-1)S]^-1 H(xo-x)
!          =x+BH'S^-1 (I+S'^-1 HBH' S^-1)S'^-1 H(xo-x)    (let L=S'^-1 H)
!          =x+BL'(I+LBL')^-1 L(xo-x)
!          =x+BL'L(I+BL'L)^-1 (xo-x)                      (let Q=BL'L)
!          =x+Q(I+Q)^-1 (xo-x)
!          =Q(I+Q)^-1 xo + (I+Q)^-1 x
!         B=...=B-Q(I+Q)^-1 B=(I+Q)^-1 B
!     The serial EnKF:
!         x_i+1=...=x_i+B_i li'li (I+B_i li'li)^-1 (xo-x_i)
!         B_i+1=...=B_i-B_i li'li (I+B_i li'li)^-1 B_i 
!     where li=si^-1 hi; Q=BL'L= Bl1'l1+Bl2'l2+...+Blp'lp = Q1+Q2+...+Qp
!     Using mathematical induction: if p=1, equivalence is true. Assume equivalence holds until the (p-1)th step.
! We can prove the equivalence by showing that p-1 -> p gives the same results as original equations.
!     L=(l1 l2 ... lp)=(Lp-1 lp), Q=Qp-1+Qp.
!     Equivalent until i=p-1: 
!         x_p=Qp-1(I+Qp-1)^-1 xo + (I+Qp-1)^-1 x_1
!         B_p=(I+Qp-1)^-1 B_1
!     p-1 -> p:
!         x_p+1 = x_p+B_p lp'lp(I+B_p lp'lp)^-1 (xo-x_p) = ...   (substitute x_p and B_p by previous 2 eqns)
!               = (I+Qp-1)^-1 Qp [I+(I+Qp-1)^-1 Qp]^-1 xo +
!                 [I+(I+Qp-1)^-1 Qp]^-1 [Qp-1(I+Qp-1)^-1 xo + (I+Qp-1)^-1 x_1]
!               = Q(I+Q)^-1 xo + (I+Q)^-1 x_1         (Q, I are symmetric, can switch places in multiplication)
!         B_p+1 = (I+B_p lp'lp)^-1 B_p = ... = (I+Q)^-1 B_1
!     Q.E.D.
!
! Square-root modification of EnKF (EnSRF):
! ----------------------------------------
!     The obs (yo) we use in assimilation is the mean value with a prescribed error, yo=yom, yoa=0 and r=error^2.
! R is prescribed, not calculated by R=<yoa yoa'>. So in the update eqns, the calculation of B is erroneous:
!         B_f = <xa_f xa_f'>
!         B_a = <xa_a xa_a'> = <(xa_f-KH xa_f) (xa_f-KH xa_f)'> = (I-KH) B_f (I-KH)'
!     it should be:
!         B_a = <(xa_f - KH xa_f + K yoa) (xa_f - KH xa_f + K yoa)'>
!             = (I-KH) B_f (I-KH)' + KRK' = (I-KH) B_f
!     This can be reconciled by changing K to aK, so that (I-aKH) B_f (I-aKH)' = (I-KH) B_f.
!     In serial EnKF, this "a" (alpha) can be solved as a scalar:
!         alpha = {1+sqrt[r/(r+hBh')]}^-1    (Whitaker and Hamill 2002)
!
! Localization:
! ------------
!     Due to deficient sampling of model error (small num of members), B contains spurious correlations between
! variables at large distances. This is reconciled with a localization technique. Gaspari and Cohn (1999) 5th-order
! correlation function is used to calculate coefficient over radius of influence (roi), from 1 at the obs location
! to 0 at roi.
! ****NOTE: the use of localization breaks the equivalence between serial EnKF and original EnKF. Now the sequence
! of obs being assimilated matters for the final analysis.
!
! Final version of update equations for x and y in this algorithm:
! ---------------------------------------------------------------
!     Update x: for i=1,2,...,p
!         B_i = <xa_i xa_i'> = SUM( xa_i xa_i' )/(m-1)
!         xm_i+1 = xm_i + B_i hi'(ri + hi B_i hi')^-1 (yo_i - hi xm_i) * localization         ...(mean)
!         xa_i+1 = xa_i + B_i hi'(ri + hi B_i hi')^-1 (   0 - hi xa_i) * alpha * localization ...(perturbation)
!     Update y: for i=1,2,...,p
!         for j=1,2,...,p, we update yf_j with the (i+1)th obs. We left-multiply the eqns for x with hj:
!         ym_j_i+1 = hj xm_i+1 = hj xm_i + hj B_i hi' (ri + hi B_i hi')^-1 (yo_i - hi xm_i) * localization
!         yf_j_i+1 = hj xa_i+1 = hj xa_i + hj B_i hi' (ri + hi B_i hi')^-1 (   0 - hi xa_i) * alpha * localization
!     where hj B_i hi' = <hj xa_i xa_i' hi'> = <yf_j_i yf_i_i'>
!     yf=Hxa, ym=Hxm is calculated all before update steps, so (h? x?_?) here in the eqns are all from y?.
!
!--------------------------------------------------------------------------------------------------------------
subroutine enkf (wrf_file,g_comm,s_comm,gid,sid,iid,jid, proj)
use constants
use namelist_define
use wrf_field_indices_define
use obs_define
use mpi_module
use netcdf
use radar 
use common_variables
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
integer, intent(in) :: sid,gid,iid,jid,g_comm,s_comm
character (len=10), intent(in) :: wrf_file 
type (proj_info)    :: proj
character (len=10)  :: obstype
! fac = 1/(n-1) used for averaging in variance calculation
! loc_fac = localization factor
! ngx, ngz = roi in horizontall and vertical
! var,cov = error variance,covariance of something
! y_hxm = y-hxm (the innovation vector, mean), hxa is the perturbation (of ensemble members).
real      :: fac,d,alpha,var,cov,y_hxm,d_ogn, loc_fac
real      :: var_a,var_b,m_d2,m_var_a,m_var_b,m_dab
real,dimension(numbers_en) :: hxa
integer   :: ngx, ngz
integer   :: i, j, k, m, n, iob, iiob, nob, ie, iunit,ounit,ii, jj, kk, is, it, ig, iv, i1,j1, itot, iit
integer   :: ist,ied,jst,jed,kst,ked, istart,iend,jstart,jend, uist,uied,ujst,ujed
integer   :: sist,sied,sjst,sjed, sistart,siend,sjstart,sjend,sis,sie,sjs,sje
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
real, dimension(ni,nj,nk,nv)                   :: x_t, std_x,std_xf
real, dimension(3,3,nk,nv,nm,nicpu*njcpu) :: xobsend, xob
real, dimension(3,3,nk,nv) :: tmp, tmpsend
!real, dimension(obs%num,numbers_en+1) :: ya, yasend, yf

real  :: ya_mean 
real, allocatable, dimension(:,:,:)   :: km1 , km1send !, km1, km1send
real, dimension( ni, nj, nk ) :: km, kmsend
real :: dx
! for satellite radiance
integer :: iob_radmin,iob_radmax
real, dimension(obs%num) :: ya_ca
real, dimension(obs%num,nm) :: yasend_tb
real, dimension(ni,nj,nk)     :: xq_n,xq_p
real, dimension(ni,nj,nk,nm)  :: xq_nsend,xq_psend
! Empirical Localization Functions
integer :: n_r=101
integer :: ich, iob_satch,num_satch
! parameter estimation
integer :: iost,update_file
! Variables for slab tagging
integer :: flag_slab_update, flag_run_xvar
integer, dimension(ix, jx ) :: slab_mask, obs_mask

read(wrf_file(6:10),'(i5)')iunit
if ( my_proc_id==0 ) then
   write(*, *)'   '
   write(*, *)'---------------------------------------------------'
   write(*, *)'.... Begin EnKF Assimilation ....'
   write(*, *)'   '
endif


! Get domain info
call open_file(wrf_file, nf_nowrite, fid)
rcode = nf_get_att_int(fid, nf_global, 'GRID_ID', grid_id)
rcode = nf_get_att_real(fid, nf_global, 'DX', dx)
dx=dx/1000.
num_update_var = 0
do m = 1, 20
   if ( len_trim(updatevar(m))>=1 ) num_update_var=num_update_var+1
enddo
call close_file(fid)


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



assimilated_obs_num = 0


! Immediate termination of subroutine if no obs are present
if ( obs%num == 0 ) return



!----------------------------------------------------------------------------------------
! III. assimilate obs
if(my_proc_id==0) write(*,*) 'Assimilating obs...'

! Description of parallelization (by Man-Yau Chan, 2020)
! ======================================================
! At this point in the code, each process contains a copy of obs and ya. Ie,
! each process knows all of the observations that needs to be assimilated, and
! the corresponding ensemble of simulated observations. Furthermore, each
! process contains a slab of the domain.
! 
! The central idea of the parallelization is simply update each process' slab
! with all the observations, using the copies of obs and ya held in said
! process. Ie, each slab is updated independent of other slabs. Furthermore, we
! will be using the joint state-observation space update in each slab.
!
! For instance, suppose we have 3 ensemble members and 4 observations (obs A, B, 
! C and D). If we only use 2 processes and 2 slabs, then each process contains a
! slab that covers half of the WRF domain. In this setup half a domain's worth of 
!WRF fields for all three members and the mean state.
!
! In each of the 2 processes:
! 1) Select an unassimilated obs (lets's say A).
! 2) Compute ensemble covariance between the slab's WRF variables and simulated 
!    ensemble's version of selected obs.
! 3) Perform Kalman filter update on all three members and the mean state using
!    the covariance computed in step 2).
! 4) Compute ensemble covariance between ensemble simulated version of selected
!    obs and the ensemble simulated versions of all 4 obs.
! 5) Using covariance computed in step 4), perform Kalman filter update on 
!    the ensemble simulated observations.
! 6) Return to step 1), and assimilate a different obs until all obs are 
!    assimilated.
! 
! As can be seen in the steps above, there is no communication between
! processes holding to different slabs. Ie, each slab is updated independently 
! of other slabs. 
! 
!
! Optimal NMCPU, NICPU and NJCPU setup:
! --------------------------------------
! For optimal use of resources, always set NMCPU=1 and set NICPU and NJCPU 
! such that NICPU x NJCPU comes very close to the maximum number of cores
! you have available. 
!
! This is optimal for two reasons.
! 1) With NMCPU=1, there is absolutely no communication between processes 
!    during the assimilation of observations. This prevents MPI communication
!    bottlenecks.
! 2) Each process contains a copy of the slab's ensemble mean state. If you sum
!    the size of the slab's ensemble mean state across all processes, the memory
!    used to store all these copies is equivalent to NMCPU x full domain mean
!    state. 
!    If we set NMCPU = ensemble size, the many copies of the full domain's mean
!    state occupies the same amount of memory as the entire ensemble.
!    Setting NMCPU = 1 should minimize the memory used by ensemble mean state.

km =0. ; kmsend=0.;


timer=MPI_Wtime()
t_update_x = 0.; t_update_y = 0.
t_ab = 0.; t_c = 0.; t_d = 0.; t_d = 0.; t_e = 0.;
t_f = 0.; t_g = 0.; t_h = 0.

obs_assimilate_cycle : do it = 1,obs%num
! 0. get basic information of the obs being assimilated
   iob=ind(it)
   obs%kick_flag(iob)=0
   obstype = obs%type(iob)
   error = obs%err(iob)


   y_hxm = obs%dat(iob) - ya(iob,numbers_en+1)
   if ( my_proc_id==0 ) write(*,'(a,i6,a,f10.2,a,f10.2,a,f8.2,a,f8.2,a,i4,a,i4,a,i4,a,i2)') &
      'No.',iob,' '//obstype//' =',obs%dat(iob), ' ya=', ya(iob,numbers_en+1), ' y-ya=', y_hxm, &
      ' err=',error,' hroi=',obs%roi(iob,1),'(',obs%roi(iob,3),') vroi=',obs%roi(iob,2),'at z',int(obs%position(iob,3))
   if( abs(y_hxm)>(error*5.) .and. &
      .not.(obstype=='min_slp   ' .or. obstype=='longitude ' .or. obstype=='latitude  ' .or. obstype=='slp       '&
       .or. obstype=='Radiance  ' .or. obstype=='Microwave ' ) ) then
       !.or. obstype=='Radiance  ' .or. obstype(1:9)=='Psounding') ) then
      if ( my_proc_id==0 ) write(*,*)' ...kicked off for large error'
      obs%kick_flag(iob)=1
      cycle obs_assimilate_cycle
   endif
   if ( obstype=='Microwave ' .and. (ya(iob,numbers_en+1) .NE. ya(iob,numbers_en+1)) ) then
      if ( my_proc_id==0 ) write(*,*)' ...kicked off for NaN'
      obs%kick_flag(iob)=1
      cycle obs_assimilate_cycle
   endif
   assimilated_obs_num=assimilated_obs_num+1
   ngx = obs%roi(iob,1)
   ngz = obs%roi(iob,2)
! Gaussian error added if using truth/idealized as obs
   if ( use_simulated .or. use_ideal_obs ) then
      call date_and_time(date, time, zone, values) 
      y_hxm = y_hxm + error*gaussdev(sum(values))
   end if
! fac*var=hBh'=<hxa hxa'>=<(ya-yam)(ya-yam)'>
! alpha is for Square-Root filter (Whitaker and Hamill 2002, MWR130,1913-1924)
! d=hBh'+r, for each obs, d reduces to a scalar, r=error^2.
   var = 0.
   do ie = 1, numbers_en
      hxa(ie) = ya(iob,ie) - ya(iob,numbers_en+1)
      var     = var + hxa(ie)*hxa(ie)
   enddo 
   fac  = 1./real(numbers_en-1) 
   d    = fac * var + error * error 
   alpha = 1.0/(1.0+sqrt(error*error/d))

  ! If doing AOEI for a radiance obs, prepare to inflate stuff
  if ( ( trim( obstype ) == 'Radiance'  .and. use_aoei ) .or. &
       ( trim( obstype ) == 'Microwave' .and. aoei_microwave ) ) then
    d = max(fac * var + error * error, y_hxm * y_hxm)
    alpha = 1.0/(1.0+sqrt((d-fac * var)/d))
    if ( my_proc_id == 0 .and. sqrt(d-fac * var) > error ) &
        write(*,*) 'observation-error inflated to ',sqrt(d-fac * var)
  endif


  ! Skipping over update_x_var completely if the slab is not updated by any of
  ! the ROIs.
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

!--------------------------------------------------------------------------------------------------
! cycle through variables to process, in x, 2D variables are stored in 3D form (with values only on
! k=1), when sending them among cpus, only the lowest layer (:,:,1) are sent and received.
  call cpu_time(t0)
  run_update_xvar: if ( flag_run_xvar > 0 ) then
    update_x_var : do m=1,nv

       call cpu_time(t1)

       ! Check if variable needs updating
       varname=enkfvar(m)
       update_flag = 0
       do iv = 1, num_update_var
         if ( varname .eq. updatevar(iv) ) update_flag = 1
       enddo

       call cpu_time(t2)
       t_ab = t_ab + t2-t1

       if ( update_flag==0 ) cycle update_x_var

       ! -------------------------------------------------------------
       ! Identify observation update zone
       ! --------------------------------------------------------------
       call cpu_time(t1)

       ! a. Identify ROI (special handling for Radiance obs)
       ! ---------------------------------------------------
       ngx = obs%roi(iob,1)
       ngz = obs%roi(iob,2)

       ! Special ROI handling for Radiance obs (ROI varies with obs)
       if (obstype=='Radiance  ') then
         ! --- Use Successive Covariance Localization (SCL) Zhang et al. (2009) MWR for BT assimilation
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
         ! --- Use Empirical Localization Functions (ELFs)
         if (use_elf) then
           call corr_elf(varname, obs%sat(iob), obs%ch(iob), ya_ca(iob), ngx, ngz, obs%position(iob,3))
           if ( my_proc_id==0 .and. (varname=='U         '.or.varname=='V         '.or.varname=='QVAPOR    ')) &
               write(*,*) varname, ya_ca(iob), ngx, ngz, obs%position(iob,3)
         endif
         if (ngx == 0 .or. ngz == 0) update_flag = 0
         ! --- ELFs end
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

       call cpu_time(t2)
       t_ab = t_ab + t2-t1


       ! If update zone does not overlap with slab, move to next slab variable
       if ( flag_slab_update == 0 ) cycle update_x_var


 

       call cpu_time(t1)
       ! c. Identify update zone in slab and allocate Kalman gain
       ! --------------------------------------------------------
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
       call cpu_time(t2)

       t_c = t_c + t2 - t1

       ! Allocate Kalman gain arrays
       !if (nmcpu > 1) then
       !allocate( km1send( uied-uist+1, ujed-ujst+1, ked-kst+1 ) )
       !allocate( km1    ( uied-uist+1, ujed-ujst+1, ked-kst+1 ) )
       !km1send = 0.; km1=0.
       !endif

       call cpu_time(t1)
       ! d. Compute Kalman gain matrix
       ! -----------------------------

       ! Compute covariance for members held by process
       do n = 1, nm
         ie = (n-1) * nmcpu + gid + 1
         if ( ie <= numbers_en ) then
           ! Fill up Kalman gain
           ! kmsend = kmsend + x( uist:uied, ujst:ujed, kst:ked, m, n ) * hxa(ie) 
           kmsend( uist:uied, ujst:ujed, kst:ked ) = &
             kmsend( uist:uied, ujst:ujed, kst:ked ) &
              + x( uist:uied, ujst:ujed, kst:ked, m, n ) * hxa(ie)
       !    km1send = km1send + x( uist:uied, ujst:ujed, kst:ked, m, n ) * hxa(ie)

         endif
       enddo

       ! Gather all necessary covariances across processes if nmcpu > 1
       if ( nmcpu > 1 ) then 
         call MPI_AllReduce( kmsend, km, ni*nj*nk, &
                             MPI_REAL, MPI_SUM, g_comm, ierr )
       else
         !write(*,*)'meep'
         km( uist:uied, ujst:ujed, kst:ked ) &
           = kmsend( uist:uied, ujst:ujed,kst:ked )
       !  km1 = km1send
       endif
                         

       ! Compute Kalman gain
       km = km * fac / d
       !km1 = km1 * fac/d
       call cpu_time(t2)
       t_d = t_d + t2 - t1

       call cpu_time(t1)
       ! e. Localize Kalman gain matrix
       ! ----------------------------------
       do k = kst, ked
         do j = ujst, ujed 
           do i = uist, uied
             ! Compute Gaspari-Cohn localization factor
             call corr( (  i + istart - 1 ) - obs%position(iob,1), &
                        (  j + jstart - 1 ) - obs%position(iob,2), &
                        k - obs%position(iob,3), ngx, ngz, loc_fac )
             ! Apply GC localization factor
             km(i, j, k) = km(i, j, k) * loc_fac
           enddo
         enddo
       enddo
       call cpu_time(t2)
       t_e  = t_e + t2 - t1

       !if (my_proc_id == 0) then
       !  write(*,*) 'km1', sum( km1)/( (uied-uist+1)*(ujed-ujst+1)*(ked-kst+1) )
       !  write(*,*) 'km ', sum( km)/( (uied-uist+1)*(ujed-ujst+1)*(ked-kst+1) )
       !  !call exit(ierr)
       !endif

       call cpu_time(t1)
       ! f. Update slab's model perturbations and mean
       ! ----------------------------------------------
       do n = 1, nm
         ie = (n-1) * nmcpu + gid + 1
         if ( ie <= numbers_en ) then
           ! Update to perturbation
           x( uist:uied, ujst:ujed, kst:ked, m, n ) = &
             x( uist:uied, ujst:ujed, kst:ked, m, n )  &
             - km( uist:uied, ujst:ujed, kst:ked) * alpha * hxa(ie)
             !- km1 * alpha * hxa(ie)
             !- km( uist:uied, ujst:ujed, kst:ked) * alpha * hxa(ie)
         endif
       enddo
       ! Update to mean state
       xm( uist:uied, ujst:ujed, kst:ked, m ) = &
         xm( uist:uied, ujst:ujed, kst:ked, m ) &
         + km( uist:uied, ujst:ujed, kst:ked) * y_hxm
         !+ km1*y_hxm
         !+ km( uist:uied, ujst:ujed, kst:ked) * y_hxm
     
       call cpu_time(t2)
       t_f = t_f + t2 - t1

      
       call cpu_time(t1)
       ! g. Special treatment for negative QVAPOR
       ! ----------------------------------------
       if ( trim(varname) .eq. 'QVAPOR' ) then
         do n = 1, nm
           ie = (n-1)*nmcpu+gid+1
           if(ie<=numbers_en) then
             ! Identifying locations where member has neg QVAPOR
             ! and turning it to zero.
             where( xm( uist:uied, ujst:ujed, kst:ked, m ) &
                    + x( uist:uied, ujst:ujed, kst:ked, m, n ) < 0 ) & 
               x( uist:uied, ujst:ujed, kst:ked, m, n ) = &
                 0. -xm( uist:uied, ujst:ujed, kst:ked, m )
           endif
         enddo
         ! Identify locations where mean field has neg QVAPOR
         where( xm( uist:uied, ujst:ujed, kst:ked, m ) < 0. ) &
             xm( uist:uied, ujst:ujed, kst:ked, m ) = 0.
       endif
       call cpu_time(t2)
       t_g = t_g + t2 - t1


       call cpu_time(t1)
       !! h. Deallocate Kalman gain and move to next variable
       !! ---------------------------------------------------
       !deallocate(km1)
       !deallocate(km1send)

       km( uist:uied, ujst:ujed, kst:ked) =0.
       kmsend( uist:uied, ujst:ujed, kst:ked) =0.

       ! Trying to prevent stale procs
       if ( nmcpu > 1 ) &
         call MPI_Barrier( g_comm, ierr)

       call cpu_time(t2)
       t_h = t_h + t2 - t1


    enddo update_x_var
  endif run_update_xvar
  call cpu_time(t2)
  t_update_x = t_update_x + t2 - t0

!--------------------------------------------------------------------------------------------------
! 6. update relevant ya
!    for members ya(:,1:n)=ya+ym (perturbation+mean), and ym=ya(:,n+1)
!    h_i B h_j'= <ya_i ya_j'> = cov(ya_i,ya_j)
!    d   = variance(ya) + obserr variance
!    ym=ym+corr_coef*hBh'(y-ym)/d
!    ya=ya+alpha*corr_coef*hBh'(0-ya)/d
!    --- basically these are the update equations of x left-multiplied by H.
   call cpu_time(t0)

   ngx = obs%roi(iob,1)
   ist = max( update_is, int(obs%position(iob,1))-ngx )
   ied = min( update_ie, int(obs%position(iob,1))+ngx )
   jst = max( update_js, int(obs%position(iob,2))-ngx )
   jed = min( update_je, int(obs%position(iob,2))+ngx )
   kst = max( update_ks, int(obs%position(iob,3))-ngz )
   ked = min( update_ke, int(obs%position(iob,3))+ngz )

   ! Saving copy of obs mean
   ya_mean = ya( iob, numbers_en+1)


   ! Only update obs that will be assimilated later.
   update_y_cycle : do iit = 1, obs%num !iiob=,obs%num

      iiob = ind( iit)
 
      ! skip those ya outside update zone
      if ( obs%position(iiob,1)<ist .or. obs%position(iiob,1)>ied .or. &
           obs%position(iiob,2)<jst .or. obs%position(iiob,2)>jed .or. & 
           obs%position(iiob,3)<kst .or. obs%position(iiob,3)>ked ) cycle update_y_cycle

      call corr(real(obs%position(iiob,1)-obs%position(iob,1)), real(obs%position(iiob,2)-obs%position(iob,2)), &
                real(obs%position(iiob,3)-obs%position(iob,3)), obs%roi(iob,1), obs%roi(iob,2), loc_fac)
      !loc_fac = 1.

      !endif
      cov=0.
      do ie=1,numbers_en
         cov     = cov + hxa(ie)*(ya(iiob,ie)-ya(iiob,numbers_en+1))
      enddo

      do ie=1,numbers_en+1
         if(ie<=numbers_en) &
            ya(iiob,ie)=ya(iiob,ie)-loc_fac*alpha*fac*cov*hxa(ie)/d !perturbation
         ya(iiob,ie)=ya(iiob,ie)+loc_fac*fac*cov*(obs%dat(iob)-ya_mean)/d        !mean
         if(obs%type(iiob)(10:10)=='Q' .and. ya(iiob,ie)<0.) ya(iiob,ie)=0.  ! remove negative values of Q
      enddo

   end do update_y_cycle
   call cpu_time(t2)
   t_update_y = t_update_y + t2 - t0

   !! calculate mean
   !yam_radiance = 0
   !do ie = 1, numbers_en
   !  yam_radiance = yam_radiance + ya(:,ie)/real(numbers_en)
   !enddo

end do obs_assimilate_cycle

! Print statement to tell folks that proc finished.
if ( use_nonlinear_enkf == .false. ) &
  write(*,'(a,x,i5,x,a,x,f8.1,x,a)') 'Proc ',my_proc_id,'finished assimilating obs in',MPI_Wtime()-timer,'seconds'


call MPI_Barrier(comm, ierr)
if ( my_proc_id==0 ) write(*,*)'Number of assimilated obs =',assimilated_obs_num



! Looking at total core time to update stuff
t_buffer = t_update_y
t_update_y = 0.
call MPI_ALLREDUCE( t_buffer, t_update_y, 1, MPI_REAL, MPI_SUM, comm,ierr)
t_buffer = t_update_x
t_update_x = 0.
call MPI_ALLREDUCE( t_buffer, t_update_x, 1, MPI_REAL, MPI_SUM, comm,ierr)

t_buffer = t_ab
t_ab = 0.
call MPI_ALLREDUCE( t_buffer, t_ab, 1, MPI_REAL, MPI_SUM, comm,ierr)
t_buffer = t_c 
t_c = 0.
call MPI_ALLREDUCE( t_buffer, t_c, 1, MPI_REAL, MPI_SUM, comm,ierr)
t_buffer = t_d 
t_d = 0.
call MPI_ALLREDUCE( t_buffer, t_d, 1, MPI_REAL, MPI_SUM, comm,ierr)
t_buffer = t_e                            
t_e = 0.                                  
call MPI_ALLREDUCE( t_buffer, t_e, 1, MPI_REAL, MPI_SUM, comm,ierr)
t_buffer = t_f                            
t_f = 0.                                  
call MPI_ALLREDUCE( t_buffer, t_f, 1, MPI_REAL, MPI_SUM, comm,ierr)
t_buffer = t_g                            
t_g = 0.                                  
call MPI_ALLREDUCE( t_buffer, t_g, 1, MPI_REAL, MPI_SUM, comm,ierr)
t_buffer = t_h                            
t_h = 0.                                  
call MPI_ALLREDUCE( t_buffer, t_h, 1, MPI_REAL, MPI_SUM, comm,ierr)



if (my_proc_id==0) then
  write(*,'(a,x,f10.1,x,a)') 'update_x_var took', t_update_x/nprocs, 'seconds per core.'
  write(*,'(a,x,f10.1,x,a)') 'update_y_cycle took', t_update_y/nprocs, 'seconds per core.'
  open( 17, file='enkf_time.log', status='old',position='append', action='write')
  write(17,'(a30, x,f8.1,x,a,/)')'update_x_var took',  t_update_x/nprocs, 'seconds per core.'
  write(17,'(a36, x,f8.1,x,a,/)')'step ab took',  t_ab/nprocs, 'seconds per core.'
  write(17,'(a36, x,f8.1,x,a,/)')'step  c took',  t_c/nprocs, 'seconds per core.'
  write(17,'(a36, x,f8.1,x,a,/)')'step  d took',  t_d/nprocs, 'seconds per core.'
  write(17,'(a36, x,f8.1,x,a,/)')'step  e took',  t_e/nprocs, 'seconds per core.'
  write(17,'(a36, x,f8.1,x,a,/)')'step  f took',  t_f/nprocs, 'seconds per core.'
  write(17,'(a36, x,f8.1,x,a,/)')'step  g took',  t_g/nprocs, 'seconds per core.'
  write(17,'(a36, x,f8.1,x,a,/)')'step  h took',  t_h/nprocs, 'seconds per core.'
  write(17,'(a30, x,f8.1,x,a,/)')'update_y_cycle took',  t_update_y/nprocs, 'seconds per core.'
  close(17)
endif

!if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' 1 tooks ', t1, ' seconds.'
!if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' 2 tooks ', t2, ' seconds.'
!if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' 3 tooks ', t3, ' seconds.'
!if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' 4 tooks ', t4, ' seconds.'
!if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' 5 tooks ', t5, ' seconds.'
if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' Assimilation tooks ', MPI_Wtime()-timer, ' seconds.'

end subroutine enkf
! --------------------------------------------------------------------------------












! ================================================================================ 
! Subroutine to print out filter diagnostics
subroutine write_filter_diagnostics(proj)

  ! Use modules
  use constants
  use namelist_define
  use wrf_field_indices_define
  use obs_define
  use mpi_module
  use netcdf
  use radar 
  use common_variables

  ! Variable definitions
  ! --------------------
  implicit none

  ! Subroutine arguments
  type ( proj_info ), intent(in) :: proj

  ! Other variables used in subroutine
  real ::var_a, var_b, errpr
  integer :: iob, ie, n, ounit, iost
  real, dimension(3) :: center_xb

  

  ! Iterate thru observations and print out the diagnostics
  write_diagnostics_loop: do iob = 1, obs%num

    ! Figure out geographical location of observation
    call ij_to_latlon(proj, obs%position(iob,1), obs%position(iob,2), center_xb(1), center_xb(2))

    ! Compute prior and posterior variances
    var_a=0.0
    var_b=0.0
    compute_sim_obs_var_loop: do ie = 1, numbers_en
      var_a=var_a+(ya(iob,ie)-ya(iob,numbers_en+1))**2
      var_b=var_b+(yf(iob,ie)-yf(iob,numbers_en+1))**2
    enddo compute_sim_obs_var_loop
    var_a=var_a/real(numbers_en-1)
    var_b=var_b/real(numbers_en-1)

    
    ! Deciding which file to write diagnostic to.
    ! (depends on whether obs was rejected)
    if(obs%kick_flag(iob)==0) then
      ounit=10000
    else
      ounit=10001
    endif


    ! Write obs diagnostic to chosen file.
    if ( my_proc_id == 0 ) then
      write(ounit,'(a,6f11.2,2i4,6f11.5)') obs%type(iob), &   !type of observation
        center_xb(1), center_xb(2), &     !latitude, longitude
        obs%position(iob,4), &            !height
        obs%position(iob,1:3), &          !grid location: i,j,k
        obs%roi(iob,1:2), &               !horizontal ROI (# of grid points), vertical ROI (# of layers)
        obs%dat(iob), &                   !observation (y^o)
        obs%err(iob), &                   !observation error variance (R)
        yf(iob,numbers_en+1), &           !observation prior mean (H\bar{x}^f)
        ya(iob,numbers_en+1), &           !observation posterior mean (H\bar{x}^a)
        sqrt(var_b), &                    !observation prior spread (sqrt{H P^f H^T})
        sqrt(var_a)                       !observation posterior spread (sqrt{H P^f H^T}) (before relaxation)
    endif


  enddo write_diagnostics_loop

end subroutine write_filter_diagnostics

! ================================================================================












! ===============================================================================
! Subroutine to perform covariance relaxation
subroutine covariance_relaxation( xf, wrf_file, g_comm, s_comm, gid, sid, iid, jid )

  ! Use modules
  use constants
  use namelist_define
  use wrf_field_indices_define
  use obs_define
  use mpi_module
  use netcdf
  use radar 
  use common_variables

  ! Variable definitions
  ! --------------------
  implicit none

  ! Subroutine arguments
  real, dimension(ni,nj,nk,nv,nm), intent(in)    :: xf
  character (len=10), intent(in) :: wrf_file 
  integer, intent(in) :: sid,gid,iid,jid,g_comm,s_comm

  ! Other variables used in subroutine
  ! Innovation statistics-related variables
  real :: m_d2, m_var_b, m_var_a, m_dab, var_a, var_b, error

  ! Variables for adaptive inflation/relaxation stuff
  real    :: sigma_infb, sigma_infa, inflate_o, inflate_b, error_tmp,inflate_update,error_update
  real, parameter :: kappa_inf = 1.03
  real, parameter :: sigma_info = 1.0
  real, parameter :: sigma_inf_min = 1.0
  real, parameter :: sigma_inf_default = 1.0

  ! For Vortex RTPP
  character(len= 1) :: ilatc, ilonc
  integer :: ilat,ilon,tcind1_i,tcind1_j
  real :: dx, tclat, tclon, dis0, dis1, r, a, Rmin, Rmax

  ! Variables for counting and ost
  integer :: iob, n, iost, i, j, ii, jj, kk, ie

  character (len=10) :: varname




  ! Variable initialization
  m_d2=0.0                ! Mean squared innovation
  m_var_b=0.0             ! Mean prior variance
  m_var_a=0.0             ! Mean posterior variance
  m_dab=0.0               ! Mean radiance increment
  n=0                     ! Number of assimilated obs.



  ! --------------------------------------------------------------------------
  ! 1) Compute useful filter diagnostics
  ! ------------------------------------
  innovation_stats_loop: do iob=1,obs%num

    ! Estimate ensemble variances
    var_a=0.0
    var_b=0.0
    do ie=1,numbers_en
      var_a=var_a+(ya(iob,ie)-ya(iob,numbers_en+1))**2
      var_b=var_b+(yf(iob,ie)-yf(iob,numbers_en+1))**2
      !if (obs%type(iob)=='Radiance  ') then
      !   var_a=var_a+(ya(iob,ie)-yam_radiance(iob))**2
      !   var_b=var_b+(yf(iob,ie)-yfm_radiance(iob))**2
      !else
      !   var_a=var_a+(ya(iob,ie)-ya(iob,numbers_en+1))**2
      !   var_b=var_b+(yf(iob,ie)-yf(iob,numbers_en+1))**2
      !endif
    enddo
    var_a=var_a/real(numbers_en-1)
    var_b=var_b/real(numbers_en-1)


    ! Compute innovation statistics (used in adaptive relaxation)
    !if(my_proc_id==0) write(*,*) iob,'=',obs%kick_flag(iob)
    if(obs%kick_flag(iob)==0) then
      n=n+1
      m_d2=m_d2+((obs%dat(iob)-yf(iob,numbers_en+1))/obs%err(iob))**2
      m_var_b=m_var_b+var_b/(obs%err(iob)**2)
      m_var_a=m_var_a+var_a/(obs%err(iob)**2)
      if (obs%type(iob)=='Radiance  ') then
         m_dab=m_dab+(obs%dat(iob)-ya(iob,numbers_en+1))*(obs%dat(iob)-yf(iob,numbers_en+1))
         error = obs%err(iob)
      endif
    endif
  
  enddo innovation_stats_loop
  m_d2=m_d2/real(n)
  m_var_b=m_var_b/real(n)
  m_var_a=m_var_a/real(n)
  m_dab=m_dab/real(n)




  ! ---------------------------------------------------------------------
  ! 2) Perform covariance relaxation
  ! --------------------------------
  if ( my_proc_id==0 ) write(*,*)'Performing covariance relaxation...'

  !---option1. relax to prior perturbation (Zhang et al. 2004)
  !x=(1.0-mixing)*x + mixing*xf  (Performed later)
  
  !---option2. relax to prior spread (Whitaker and Hamill 2012), adaptive version (ACR)
  !------inflation factor from innovation statistics
  if(my_proc_id==0) write(*,*) 'lambda=',sqrt((m_d2-1.)/m_var_b), ' kappa=',(sqrt(m_var_b)-sqrt(m_var_a))/sqrt(m_var_a)
  if(my_proc_id==0) write(*,*) 'sigma_o=',sqrt(m_dab)
  inflate_o = ((m_d2-1.)/m_var_b)
  error_update = sqrt(m_dab)
  if (inflate_o .lt. 1.) inflate_o = 1.
  !------ensemble spread in state space
  !std_x =0.0
  !std_xf=0.0
  !call MPI_Allreduce(sum( x**2,5),std_x, ni*nj*nk*nv,MPI_REAL,MPI_SUM,g_comm,ierr)
  !call MPI_Allreduce(sum(xf**2,5),std_xf,ni*nj*nk*nv,MPI_REAL,MPI_SUM,g_comm,ierr)
  !std_x =sqrt( std_x/real(numbers_en-1))
  !std_xf=sqrt(std_xf/real(numbers_en-1))
  !do n=1,nm
  !  ie=(n-1)*nmcpu+gid+1
  !  if(ie==numbers_en+1) &
  !    x(:,:,:,:,n)=x(:,:,:,:,n)*(mixing*(std_xf-std_x)/std_x+1)
  !enddo
  
  !!! adaptively update inflation factor 2017.9.5
  ! reading previous time cycle estimation
  open (11, file = 'parameters_update', status = 'old', form = 'formatted', iostat = iost )
  if( iost .eq. 0 ) then
    read(11,fmt='(f12.4,f12.4,f12.4)')inflate_b, error_tmp, sigma_infb
    inflate_b = inflate_b**2
    sigma_infb = sigma_infb * kappa_inf ! spread increase with forecast
  else
    inflate_b = inflate
    sigma_infb = sigma_inf_default
  endif
  if (sigma_infb .lt. sigma_inf_min) sigma_infb = sigma_inf_min
  ! scalar KF 
  inflate_update = (sigma_info*inflate_b + sigma_infb*inflate_o) / (sigma_info + sigma_infb)
  sigma_infa = (1. - sigma_infb/(sigma_info + sigma_infb))*sigma_infb
  if ( my_proc_id == 0 ) write(*,'(a,f12.4,a,f12.4)')' adaptivelty updated inflation =', inflate_update
  !!!
  
  if ( my_proc_id == 0 ) then
    open(11,file='parameters_update'//times(1:4)//times(6:7)//times(9:10)//times(12:13)//times(15:16))
    write(11,'(f12.4,f12.4,f12.4)')sqrt(inflate_update),error_update,sigma_infa
    close(11)
  endif
  
  
  !Add the mean back to the analysis field
  if ( my_proc_id==0 ) write(*,*)'Add the mean back to the analysis field...'
  if ( use_vortex_RTPP ) then !by Robert Nystrom 2018.4.7
      if ( my_proc_id==0 ) write(*,*)'Using vortex RTPP...'
      ! Get storm center from "tcvitals.dat"
      open(10, file="tcvitals.dat", status='old', form = 'formatted', iostat = iost )
      if ( iost .ne. 0 ) then
          stop 'NO a-deck file found!'
      endif
      read(10,'(32x, i4, a1, i5, a1)')ilat,ilatc,ilon,ilonc
      tclat=ilat/10.
      if (ilatc.eq.'S')tclat=-ilat/10.
      tclon=ilon/10.
      if (ilonc.eq.'W')tclon=-ilon/10.
      if ( my_proc_id==0 ) write(*,*)'TC center:', tclat, tclon
      ! Map the TC center to the WRF domain 
      if ( tclat.gt.xlat(1,1) .and. tclat.lt.xlat(1,jx) .and.      &
         tclon.gt.xlong(1,1) .and. tclon.lt.xlong(ix,1) ) then
          dis0=999999.
          do j = 1, jx
              do i = 1, ix
                  dis1 = (tclat-xlat(i,j))**2 + (tclon-xlong(i,j))**2
                  if ( dis1 <= dis0 ) then
                      dis0 = dis1
                      tcind1_i = i
                      tcind1_j = j
                  end if
              end do
          end do
      else
          tcind1_i = -9999
          tcind1_j = -9999
      end if
      if ( my_proc_id==0 ) write(*,*)'TC center in WRFout:', tcind1_i, tcind1_j 
      ! Perform RTPP relative to vortex center
      Rmin=300 !km
      Rmax=600 !km
      call wrf_var_dimension ( wrf_file, varname, ix, jx, kx, kxs, ii, jj, kk )
      do n = 1, nm
          ie = (n-1)*nmcpu+gid+1
          if ( ie<=numbers_en+1 ) then
              call wrf_var_dimension ( wrf_file, varname, ix, jx, kx, kxs, ii, jj, kk )
              do j = 1, jj
                  do i = 1, ii
                      r = (((i-tcind1_i)**2+(j-tcind1_j)**2)**0.5)*dx
                      if ( r > Rmax ) then
                          a = mixing
                      else if ( r < Rmin ) then
                          a = mixing_core
                      else if ( r <= Rmax .and. r >= Rmin )then
                          a = mixing-(mixing-mixing_core)*((Rmax-r)/(Rmax-Rmin))
                      endif
                      !if ( my_proc_id==0 ) write(*,*)'i,j,r,alpha:',i, j, r, a
                      x(i,j,:,:,n) = (1.0-a)*x(i,j,:,:,n) + a*xf(i,j,:,:,n) + xm(i,j,:,:)
                  end do
              end do
          endif
      enddo    
  else !Regular uniform RTPP
      if (my_proc_id ==0) then
        write(*,*) 'Pre RTPP stats'
        write(*,*) '|U10| xf avg', sum( abs(xf(:,:,1,ind_u10,numbers_en+1)) )/(ni*nj), sum( xf(:,:,1,ind_u10,1:numbers_en) )/(ni*nj*numbers_en)
        write(*,*) '|U10| x avg' , sum( abs( x(:,:,1,ind_u10,numbers_en+1)) )/(ni*nj), sum(  x(:,:,1,ind_u10,1:numbers_en) )/(ni*nj*numbers_en)
        write(*,*) 'lvl 1 |U| xf avg', sum( abs(xf(:,:,1,ind_u,numbers_en+1)) )/(ni*nj), sum( xf(:,:,1,ind_u,1:numbers_en) )/(ni*nj*numbers_en)
        write(*,*) 'lvl 1 |U| x avg' , sum( abs( x(:,:,1,ind_u,numbers_en+1)) )/(ni*nj), sum(  x(:,:,1,ind_u,1:numbers_en) )/(ni*nj*numbers_en)
      endif

      do n = 1, nm
          ie = (n-1)*nmcpu+gid+1
          if( ie<=numbers_en+1 )  &
          x(:,:,:,:,n) = (1.0-mixing)*x(:,:,:,:,n) + mixing*xf(:,:,:,:,n) + xm
      enddo
  endif

end subroutine covariance_relaxation
! ================================================================================











! ================================================================================
! Subroutine to remove negative Q values (originally by Minamide 2015.5.26)
subroutine minamide_neg_q_removal( g_comm,s_comm,gid,sid,iid,jid )

  ! Use modules
  use constants
  use namelist_define
  use wrf_field_indices_define
  use obs_define
  use mpi_module
  use netcdf
  use radar 
  use common_variables

  ! Variable definitions
  ! --------------------
  implicit none

  ! Subroutine arguments
  integer, intent(in) :: sid,gid,iid,jid,g_comm,s_comm

  ! Variables to hold negative and positive variables
  real, dimension(ni,nj,nk)     :: xq_n,xq_p
  real, dimension(ni,nj,nk,nm)  :: xq_nsend,xq_psend

  character (len=10) :: varname

  ! Variables for counting
  integer :: m, n, ie


  if ( my_proc_id==0 ) write(*,*) 'updating negative values'

  if(raw%radiance%num.ne.0) then
   do m=1,nv
     varname=enkfvar(m)
     xq_p = 0.
     xq_n = 0.
     xq_psend = 0.
     xq_nsend = 0.
     if (varname=='QCLOUD    ' .or. varname=='QRAIN     ' .or. varname=='QICE      ' .or. &
         varname=='QGRAUP    ' .or. varname=='QSNOW     ') then
      do n=1,nm
       ie=(n-1)*nmcpu+gid+1
       if(ie==numbers_en+1 .and. print_detail > 0 ) write(*,*)varname, ' original xq value',minval(x(:,:,:,m,n)),'~',maxval(x(:,:,:,m,n))
       if(ie<=numbers_en) then
         where(x(:,:,:,m,n) >= 0.) xq_psend(:,:,:,n) = x(:,:,:,m,n)
         where(x(:,:,:,m,n) < 0.) xq_nsend(:,:,:,n) = x(:,:,:,m,n)
       endif
      enddo
      call MPI_Allreduce(sum(xq_psend,4),xq_p,ni*nj*nk,MPI_REAL,MPI_SUM,g_comm,ierr)
      call MPI_Allreduce(sum(xq_nsend,4),xq_n,ni*nj*nk,MPI_REAL,MPI_SUM,g_comm,ierr)
      if ( my_proc_id==0 ) write(*,*) 'xq_p',minval(xq_p),'~',maxval(xq_p)
      if ( my_proc_id==0 ) write(*,*) 'xq_n',minval(xq_n),'~',maxval(xq_n)
      do n=1,nm
       ie=(n-1)*nmcpu+gid+1
       if(ie<=numbers_en) then
         where(x(:,:,:,m,n) < 0.) x(:,:,:,m,n) = 0.
         where(xq_p >= abs(xq_n).and.xq_p > 0.) x(:,:,:,m,n) = x(:,:,:,m,n)*(xq_p+xq_n)/xq_p
         where(xq_p < abs(xq_n).or. xq_p == 0.) x(:,:,:,m,n) = 0.
       endif
       if(ie<=numbers_en+1) where(xm(:,:,:,m) < 0.) x(:,:,:,m,n) = 0.
       if(ie==numbers_en+1 .and. print_detail > 0) write(*,*)'non-negative xq value',minval(x(:,:,:,m,n)),'~',maxval(x(:,:,:,m,n))
      enddo
     endif
   enddo
  endif

end subroutine minamide_neg_q_removal



