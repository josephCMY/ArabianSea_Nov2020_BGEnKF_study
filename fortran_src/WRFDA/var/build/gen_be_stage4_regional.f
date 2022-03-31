












program gen_be_stage4_regional

   use da_control, only : filename_len,stderr,stdout,do_normalize,use_rf
   use da_reporting, only : da_error
   use da_tools_serial, only : da_get_unit, da_advance_cymdh
   use da_wavelet, only: lf,namw,nb,nij,ws,wsd

   implicit none

   character*10        :: start_date, end_date       
   character*10        :: date, new_date             
   character*10        :: variable                   
   character(len=filename_len)       :: run_dir      
   character(len=filename_len)       :: filename     
   character*1         :: k_label                    
   character*2         :: ck                         
   integer             :: ni, nj                     
   integer             :: k                          
   integer             :: stride                     
   integer             :: nbins                      
   integer             :: ibin                       
   integer             :: nn                         
   integer             :: sdate, cdate, edate        
   integer             :: interval                   
   integer             :: ne                         
   integer             :: member                     
   integer             :: jstart, jend               
   integer             :: i                          
   logical             :: print_wavelets
   real                :: rcount                     
   real                :: fnorm                      
   integer, allocatable:: nr(:)                      
   real, allocatable   :: field(:,:)                 
   real, allocatable   :: grid_box_area(:,:)         
   real, allocatable   :: sd(:,:)                    
   real, allocatable   :: cov(:)                     

   real, allocatable   :: wav(:,:)                   

   namelist / gen_be_stage4_regional_nl / start_date, end_date, interval, variable, &
                                          ne, k, nbins, ibin, stride, run_dir, &
                                          do_normalize, print_wavelets, lf, namw, nb, use_rf

   integer :: ounit,iunit,namelist_unit

   stderr = 0
   stdout = 6

   call da_get_unit(ounit)
   call da_get_unit(iunit)
   call da_get_unit(namelist_unit)


   write(6,'(a)')' [1] Initialize namelist variables and other scalars.'


   start_date = '2004030312'
   end_date = '2004033112'
   interval = 24
   variable = 'psi'
   ne = 1
   k = 1
   if( use_rf )then
      nbins = 1
      ibin = 1
      stride = 1
   endif
   run_dir = '/mmmtmp/dmbarker'
   open(unit=namelist_unit, file='gen_be_stage4_regional_nl.nl', &
        form='formatted', status='old', action='read')
   read(namelist_unit, gen_be_stage4_regional_nl)
   close(namelist_unit)

   read(start_date(1:10), fmt='(i10)')sdate
   read(end_date(1:10), fmt='(i10)')edate
   if( use_rf )then
      write(UNIT=6,FMT='(4a)')' Computing error correlation scales for dates ', start_date, ' to ', end_date
      write(UNIT=6,FMT='(a,i3,a,i3)') '                                    for bin ', ibin, ' of ', nbins
      write(UNIT=6,FMT='(a,i8)')' Stride over which to jump points in correlation calculation = ', stride
   else
      call da_error("gen_be_stage4_regional.b",84,(/"needs qftest_w, compile after setenv WAVELET 1"/))
      write(6,'(" Computing field and wavelet std. devs. for dates ",a," to ",a)')start_date,end_date
      write(6,'("    using ",i8," bands of ",a,i0.2," wavelet.")')nb,namw,lf





      allocate(nij(0:nb,0:1,0:2))
   endif
   write(UNIT=6,FMT='(a,i8,a)')' Interval between dates = ', interval, 'hours.'
   write(UNIT=6,FMT='(a,i8)')' Number of ensemble members at each time = ', ne

   write(UNIT=ck,FMT='(i2)') k
   if ( k < 10 ) ck = '0'//ck(2:2)
   k_label = 'm'

   call initialize_dates(date,start_date,cdate,sdate,rcount)


   if( use_rf )then
      write(6,'(" [2.0] Input fields and calculate correlation as a function of distance between points.")')
   else
      write(6,'(" [2.0] Input fields and calculate wavelet mean and mean-square values.")')
      filename=trim(run_dir)//'/grid_box_area'
      open(iunit,file=filename,form='unformatted')
      read(iunit)ni,nj
      allocate(grid_box_area(ni,nj))
      read(iunit)grid_box_area
      close(iunit)
   endif


   if( do_normalize )then
      do while( cdate<=edate )			
         do member = 1, ne			
            call read_ni_nj(iunit,member,run_dir,variable,date(1:10),ck,ni,nj)
            if( rcount==1.0 )call allocate_field_moments(field,wav,sd,ni,nj)
            read(iunit)field
            close(iunit)
            call update_moments(field,wav,sd,ni,nj,rcount,1.)
         enddo					
         call calculate_next_date(date,interval,new_date,cdate)
      enddo					

      write(6,'(" [2.1] Compute field standard deviations.")')

      call update_moments(field,wav,sd,ni,nj,rcount,0.)


      deallocate(field,wav)			
      call initialize_dates(date,start_date,cdate,sdate,rcount)
   endif					

   do while ( cdate <= edate )
      do member = 1, ne

         write(UNIT=6,FMT='(5a,i4)')'    Date = ', date, ', variable ', trim(variable), &
                           ' and member ', member
         call read_ni_nj(iunit,member,run_dir,variable,date(1:10),ck,ni,nj)
         if( use_rf )then
            jstart = floor(real(ibin-1)*real(nj)/real(nbins))+1
            jend   = floor(real(ibin)  *real(nj)/real(nbins))
            write(UNIT=6,FMT='(a,i4,a,i4)') 'Computing length scale over j range ', jstart, ' to ', jend
         else
            write(6,'("Computing wavelet coefficients")')
         endif
         if ( rcount == 1.0 ) then
            if( use_rf )then
               nn = ni * ni + nj * nj		
               allocate(field(1:ni,1:nj))
               allocate(nr(0:nn))
               allocate(cov(0:nn))
               cov(0:nn) = 0.0
               call get_grid_info( ni, nj, nn, stride, nr, jstart, jend )
            else				
               nij(0,0,0)=nj			
               nij(0,1,0)=ni			
               do i=1,0,-1			
                  call da_error("gen_be_stage4_regional.b",167,(/"needs dwta_partition, compile after setenv WAVELET 1"/))
                  write(6,&
                        '("Direction ",a," has partition nij(0:",i1,",",i1,",0:2)={",99("{",i3,",",i3,",",i3,"}",:))',&
                        ADVANCE="NO")char(106-i),nb,i,transpose(nij(:,i,:))
                  write(6,'("}.")')
               enddo
               call allocate_field_moments(field,wav,wsd,nij(0,1,2),nij(0,0,2))
               allocate(ws(maxval(nij(0,:,2))))
            endif
         endif					

         read(UNIT=iunit)field(1:ni,1:nj)
         close(UNIT=iunit)

         if( do_normalize )field(:ni,:nj)=field(:ni,:nj)/sd
         if( use_rf )then

            call get_covariance( ni, nj, nn, stride, rcount, nr, jstart, jend, field(:ni,:nj), cov )
            rcount=rcount+1.0
         else
            field(:ni,:nj)=field(:ni,:nj)*sqrt(grid_box_area)
            fnorm=sqrt(sum(field(:ni,:nj)**2))
            print'("r=||    ",a,"(",a,",",i0,",:",i0,",:",i0,") ||  : ",ES15.9)', &
                  trim(variable),date,member,ni,nj,fnorm

            call da_error("gen_be_stage4_regional.b",196,(/"needs dwtai2_w, compile after setenv WAVELET 1"/))
            print'("1-||  W[",a,"(",a,",",i0,",:",i0,",:",i0,")]||/r: ",ES15.9)', &
                  trim(variable),date,member,nij(0,1,2),nij(0,0,2),1.-sqrt(sum(field**2))/fnorm
            call update_moments(field,wav,wsd,nij(0,1,2),nij(0,0,2),rcount,1.)
         endif

      enddo					
      call calculate_next_date(date,interval,new_date,cdate)
   end do					

   if( use_rf )then

      write(UNIT=6,FMT='(a)')' [3] Compute fit of correlation to a straight line.'

      filename = 'sl_print.'//trim(variable)//'.'//ck
      open( unit=ounit, file=trim(filename), form='formatted', iostat=i, &
            action='write', access='sequential', status='replace')
      if( i/=0 )print'("OPEN(FILE=",a,",IOSTAT=",i0,")")',trim(filename),i
      call make_scale_length( variable, ck, ounit, nn, nr, cov )
      close(unit=ounit, status='keep')
      deallocate(cov)
      if( do_normalize )then
         write(ck,'(i0)')k
         open(ounit,action="write",file="mom",form="unformatted",status="replace")
         write(ounit)do_normalize
         write(ounit)nj,ni
         write(ounit)sd
         close(ounit)
      endif
   else 

      write(6,'(" [3] Compute wavelet standard deviations.")')

      call update_moments(field,wav,wsd,nij(0,1,2),nij(0,0,2),rcount,0.)
      print'(es9.2,"<wsd[",a,",k=",i0,"]~",es9.2,"<",es9.2)', &
         minval(wsd),trim(variable),k,sqrt(sum(wsd**2)/product(nij(0,:,2))),maxval(wsd)
      open(ounit,action="write",file="momw",form="unformatted",status="replace")

      write(ounit)do_normalize,namw,lf,nb
      write(ounit)nij
      write(ounit)wsd
      if( do_normalize )write(ounit)sd
      close(ounit)
      if( print_wavelets )then
         do i=1,0,-1			
            deallocate(ws)
            allocate(cov(nij(0,i,2)),ws(2*floor(.5*(nij(0,i,0)+lf))+lf-2))
            open(ounit,action="write",file="w_"//char(106-i)//".dat",form="unformatted",status="replace")
            write(ounit)namw,lf,nb,nij(:,i,:)
            do k=1,nij(0,i,2)		
               cov=0.			
               cov(k)=1.		
               write(6,*) "Needs idwtai_w, please compile code with WAVELET setting"
               stop
               write(ounit)real(cov(:nij(0,i,0)),4)
            enddo			
            close(ounit)
            deallocate(cov)
         enddo
      endif 

      deallocate(grid_box_area,nij,wav,ws,wsd)
      if( do_normalize )deallocate(sd)
   endif
   deallocate(field) 

contains

 SUBROUTINE allocate_field_moments(field,av,sd,ni,nj)
 IMPLICIT NONE
 INTEGER,         INTENT(IN   )::ni,nj
 REAL,ALLOCATABLE,INTENT(INOUT)::av(:,:),field(:,:),sd(:,:)
 ALLOCATE(av(ni,nj),field(ni,nj),sd(ni,nj))
 av=0.
 sd=0.
 ENDSUBROUTINE allocate_field_moments

 SUBROUTINE calculate_next_date(date,interval,new_date,cdate)
 IMPLICIT NONE
 CHARACTER*10,INTENT(INOUT)::date
 CHARACTER*10,INTENT(  OUT)::new_date
 INTEGER,     INTENT(IN   )::interval
 INTEGER,     INTENT(  OUT)::cdate
 CALL da_advance_cymdh(date,interval,new_date)
 date=new_date
 READ(date,'(I10)')cdate
 ENDSUBROUTINE calculate_next_date

 SUBROUTINE initialize_dates(date,start_date,cdate,sdate,rcount)
 IMPLICIT NONE
 CHARACTER*10,INTENT(IN )::start_date
 CHARACTER*10,INTENT(OUT)::date
 INTEGER,INTENT(IN )::sdate
 INTEGER,INTENT(OUT)::cdate
 REAL,   INTENT(OUT)::rcount
 date=start_date
 cdate=sdate
 rcount=1.0
 ENDSUBROUTINE initialize_dates

 SUBROUTINE update_moments(field,av,sd,ni,nj,rcount,irc)
 IMPLICIT NONE
 INTEGER,INTENT(IN   ):: ni,nj
 REAL,   INTENT(IN   ):: field(ni,nj),irc
 REAL,   INTENT(INOUT):: av(ni,nj),rcount,sd(ni,nj)
 REAL                 :: c
 IF( irc==1. )THEN
    av=av+field
    sd=sd+field**2			
    rcount=rcount+irc
 ELSE
    rcount=rcount-1.			
    av=av/rcount
    sd=(sd-rcount*av**2)/(rcount-1.)
    c=COUNT(sd<0.)
    IF(c>0.) &
       PRINT'(a,": rcount=",F2.0,", zeroing 0 > sd**2 > ",ES9.2," (",ES8.2,")%.")', &
          "gen_be_stage4_regional.b",rcount,MINVAL(sd),100*c/(ni*nj)

    sd=SQRT(MAX(0.,sd))
 ENDIF
 ENDSUBROUTINE update_moments

 SUBROUTINE read_ni_nj(iunit,member,run_dir,variable,date,ck,ni,nj)
 IMPLICIT NONE 
 CHARACTER(LEN=filename_len),INTENT(IN):: run_dir	
 CHARACTER*2,                INTENT(IN):: ck		
 CHARACTER*10,               INTENT(IN):: date,variable
 INTEGER,                    INTENT(IN):: iunit,member
 INTEGER,                   INTENT(OUT):: ni,nj
 CHARACTER*3                           :: ce		
 CHARACTER(LEN=filename_len)           :: filename	
 INTEGER                               :: kdum
 WRITE(ce,'(i3.3)')member

 filename=TRIM(run_dir)//'/'//TRIM(variable)//'/'//date//'.'//TRIM(variable) &
                                    //'.e'//ce//'.'//ck
 OPEN(iunit,FILE=filename,FORM="UNFORMATTED")
 READ(iunit)ni,nj,kdum
 ENDSUBROUTINE read_ni_nj

subroutine get_grid_info( ni, nj, nn, stride, nr, jstart, jend )

   implicit none

   integer, intent(in)    :: ni, nj                  
   integer, intent(in)    :: nn                      
   integer, intent(in)    :: stride                  
   integer, intent(in)    :: jstart, jend            
   integer, intent(out)   :: nr(0:nn)                

   integer                :: i, j, m, n              
   integer                :: d2                      

   nr(0:nn) = 0



!$OMP PARALLEL DO PRIVATE(I,J,N,D2) REDUCTION(+:NR)
   do j = jstart, jend, stride                         
      do i = 1, ni, stride                             




         n = j
         do m = i, ni
            d2 = (m-i) * (m-i)                         
            nr(d2) = nr(d2) + 1                        
         end do


         do n = j+1, jend
            do m = 1, ni
               d2 = (m-i)*(m-i) + (n-j)*(n-j)          
               nr(d2) = nr(d2) + 1                     
            end do 
         end do 
      end do 
   end do 
!$OMP END PARALLEL DO

end subroutine get_grid_info

subroutine get_covariance( ni, nj, nn, stride, count, nr, jstart, jend, field, cov )
   



   implicit none

   integer, intent(in)    :: ni, nj                  
   integer, intent(in)    :: nn                      
   integer, intent(in)    :: stride                  
   real, intent(in)       :: count                   
   integer, intent(in)    :: nr(0:nn)                
   integer, intent(in)    :: jstart, jend            
   real, intent(in)       :: field(1:ni,1:nj)        
   real, intent(inout)    :: cov(0:nn)               

   integer                :: i, j, m, n              
   integer                :: d2                      
   real                   :: count_inv               
   real                   :: bb(0:nn)                

   count_inv = 1.0 / count

   bb(0:nn) = 0.0

!$OMP PARALLEL DO PRIVATE(I,J,N,D2) REDUCTION(+:BB)
   do j = jstart, jend, stride                         
      do i = 1, ni, stride                             




         n = j
         do m = i, ni
            d2 = (m-i) * (m-i)                         
            bb(d2) = bb(d2) + field(i,j) * field(m,n)
         end do


         do n = j+1, jend
            do m = 1, ni
               d2 = (m-i)*(m-i) + (n-j)*(n-j)          
               bb(d2) = bb(d2) + field(i,j) * field(m,n)
            end do 
         end do 

      end do 
   end do 
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(D2)
   do d2 = 0, nn
      if ( nr(d2) /= 0 ) then
         bb(d2) = bb(d2) / real(nr(d2))


         cov(d2) = ( (count - 1.0) * cov(d2) + bb(d2) ) * count_inv

      end if
   end do
!$OMP END PARALLEL DO

end subroutine get_covariance

subroutine make_scale_length( variable, ck, ounit, nn, nr, cov )





   implicit none

   character*10, intent(in):: variable               
   character*2, intent(in):: ck                      
   integer, intent(in)    :: ounit                   
   integer, intent(in)    :: nn                      
   integer, intent(in )   :: nr(0:nn)                
   real, intent(in)       :: cov(0:nn)               

   real(kind=8)           :: yr(0:nn)                
   real(kind=8)           :: nrr(0:nn)               
   real(kind=8)           :: d(0:nn)                 

   integer                :: n, d2                   
   integer                :: nmax                    
   real                   :: yr_cutoff               
   real                   :: corr_min                
   real(kind=8)           :: coeff1, coeff2          
   real(kind=8)           :: ml, sl                  

   yr(0:nn) = 0.0
   nrr(0:nn) = 0.0
   d(0:nn) = 0.0

   yr_cutoff = 3.0                                   
   corr_min = exp( -0.125 * yr_cutoff * yr_cutoff )  
   write(UNIT=6,FMT='(a,1pe15.5)')' Fit Gaussian curve to data for correlations >= ', corr_min

   write(UNIT=6,FMT='(5a)')'  n  ', ' nr(n) ', ' d(n) ', ' cov(n) ', ' yr(n)  '
   n = 0
   nrr(n) = real(nr(n))
   write(UNIT=6,FMT='(i6,4e13.5)') n, nrr(n), d(0), cov(n), yr(n)

   do d2 = 1, nn

      if ( nr(d2) > 0 .and. cov(d2) < cov(0) ) then 

         if ( cov(d2) / cov(0) < corr_min ) exit 
         n = n + 1
         yr(n) = sqrt( 8.0 * log(cov(0) / cov(d2)) )
         nrr(n) = real(nr(d2))
         d(n) = sqrt(real(d2))                  
         write(UNIT=6,FMT='(i6,4e13.5)') n, nrr(n), d(n), cov(d2), yr(n)
       end if
   end do
   nmax = n



















   coeff1 = 0.0
   coeff2 = 0.0

   do n = 1, nmax
      WRITE(UNIT=6,FMT='("n, nrr, d, yr:",i3,3e15.6)') n, nrr(n), d(n), yr(n)
      coeff1 = coeff1 + nrr(n) * d(n) * yr(n)
      coeff2 = coeff2 + nrr(n) * d(n) * d(n)
   end do

   if (coeff2 > 0.0) then

     ml = coeff1 / coeff2

     sl = 1.0 / ml

   else



     ml = 0.0
     sl = 0.0

   endif

   write(UNIT=6,FMT='(/3a,3e30.8/)') trim(variable), ' scale-length at mode:', ck, ml, sl

   write(UNIT=ounit,FMT='(a,2e20.8)') ck, ml, sl

end subroutine make_scale_length


end program gen_be_stage4_regional
