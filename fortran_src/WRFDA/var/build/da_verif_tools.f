














module da_verif_tools
   
   
   implicit none
   
   
   

   real, private, parameter    :: pi = 3.1415926
   real, private, parameter    :: deg_per_rad = 180.0/pi
   real, private, parameter    :: rad_per_deg = pi / 180.0
   logical, private, parameter :: print_detail_map = .false.

   integer, private, parameter :: stdout = 6

   
   
   
   real, public , parameter    :: earth_radius_m = 6370000.0
   real, public , parameter    :: radians_per_degree = pi / 180.0

   
 
   
   integer, public, parameter  :: PROJ_LATLON = 0
   integer, public, parameter  :: PROJ_MERC = 1
   integer, public, parameter  :: PROJ_LC = 3
   integer, public, parameter  :: PROJ_PS = 5

   
  

   type proj_info
      integer          :: code     
      real             :: lat1     
      real             :: lon1     
      real             :: dx       
                                   
      real             :: dlat     
      real             :: dlon     
      real             :: stdlon   
      real             :: truelat1 
      real             :: truelat2 
      real             :: hemi     
      real             :: cone     
      real             :: polei    
      real             :: polej    
      real             :: rsw      
      real             :: rebydx   
      real             :: knowni   
      real             :: knownj   
      real             :: latinc   
      real             :: loninc   
      logical          :: init     
                                 
   end type proj_info

   type(proj_info) :: map_info

   character(len=10000),private :: message(50)


contains

subroutine da_llxy_latlon(lat, lon, proj, x, y)

   
   
   

   implicit none

   real, intent(in)             :: lat
   real, intent(in)             :: lon
   type(proj_info), intent(in)  :: proj
   real, intent(out)            :: x
   real, intent(out)            :: y

   real                         :: deltalat
   real                         :: deltalon
   real                         :: lon360
   real                         :: latinc
   real                         :: loninc

   
   
   

   latinc = proj%truelat1
   loninc = proj%stdlon

   
   

   deltalat = lat - proj%lat1

   
   
   if (lon < 0) then 
      lon360 = lon + 360.0 
   else 
      lon360 = lon
   end if    
   deltalon = lon360 - proj%lon1      
    
   
   x = deltalon/loninc + 1.0
   y = deltalat/latinc + 1.0

end subroutine da_llxy_latlon


subroutine da_llxy_lc(lat, lon, proj, x, y)

   
   
   
   
    
   implicit none

   real, intent(in)              :: lat      
   real, intent(in)              :: lon      
   type(proj_info),intent(in)    :: proj     

   real, intent(out)             :: x        
   real, intent(out)             :: y        

   real                          :: arg
   real                          :: deltalon
   real                          :: tl1r
   real                          :: rm
   real                          :: ctl1r

   
   
   deltalon = lon - proj%stdlon
   if (deltalon > +180.0) deltalon = deltalon - 360.0
   if (deltalon < -180.0) deltalon = deltalon + 360.0
    
   
   tl1r = proj%truelat1 * rad_per_deg
   ctl1r = COS(tl1r)     
   
   
   rm = proj%rebydx * ctl1r/proj%cone * &
       (TAN((90.0*proj%hemi-lat)*rad_per_deg/2.0) / &
        TAN((90.0*proj%hemi-proj%truelat1)*rad_per_deg/2.0))**proj%cone

   arg = proj%cone*(deltalon*rad_per_deg)
   x = proj%polei + proj%hemi * rm * Sin(arg)
   y = proj%polej - rm * COS(arg)

   
   
   
   
   
   if (proj%hemi == -1.0) then
      x = 2.0 - x  
      y = 2.0 - y
   end if

end subroutine da_llxy_lc


subroutine da_llxy_merc(lat, lon, proj, x, y)

   
   
   
  
   implicit none

   real, intent(in)              :: lat
   real, intent(in)              :: lon
   type(proj_info),intent(in)    :: proj
   real,intent(out)              :: x
   real,intent(out)              :: y
   real                          :: deltalon

   deltalon = lon - proj%lon1
   if (deltalon < -180.0) deltalon = deltalon + 360.0
   if (deltalon > 180.0) deltalon = deltalon - 360.0
   x = 1.0 + (deltalon/(proj%dlon*deg_per_rad))
   y = 1.0 + (ALOG(TAN(0.5*((lat + 90.0) * rad_per_deg)))) / &
           proj%dlon - proj%rsw

end subroutine da_llxy_merc


subroutine da_llxy_ps(lat,lon,proj,x,y)

   
   
   
   
   
   

   implicit none

   real, intent(in)               :: lat
   real, intent(in)               :: lon
   type(proj_info),intent(in)     :: proj

   real, intent(out)              :: x 
   real, intent(out)              :: y 
   
   real                           :: reflon
   real                           :: scale_top
   real                           :: ala
   real                           :: alo
   real                           :: rm

   reflon = proj%stdlon + 90.0
   
   

   scale_top = 1.0 + proj%hemi * Sin(proj%truelat1 * rad_per_deg)

   
   ala = lat * rad_per_deg
   rm = proj%rebydx * COS(ala) * scale_top/(1.0 + proj%hemi *Sin(ala))
   alo = (lon - reflon) * rad_per_deg
   x = proj%polei + rm * COS(alo)
   y = proj%polej + proj%hemi * rm * Sin(alo)

end subroutine da_llxy_ps


subroutine da_llxy_wrf(proj, lat, lon, x, y)

   
   
   
   

   implicit none

   type(proj_info), intent(in)  :: proj
   real,            intent(in)  :: lat
   real,            intent(in)  :: lon
   real,            intent(out) :: x
   real,            intent(out) :: y

   if (.NOT.proj%init) then
      write(stdout,*) "da_llxy_wrf.inc"// &
           "You have not called map_set for this projection!"
      stop
   end if

   select case(proj%code)
 
      case(PROJ_LATLON)
         call da_llxy_latlon(lat,lon,proj,x,y)

      case(PROJ_MERC)
         call da_llxy_merc(lat,lon,proj,x,y)
         x = x + proj%knowni - 1.0
         y = y + proj%knownj - 1.0

      case(PROJ_PS)
         call da_llxy_ps(lat,lon,proj,x,y)
      
      case(PROJ_LC)
         call da_llxy_lc(lat,lon,proj,x,y)
         x = x + proj%knowni - 1.0
         y = y + proj%knownj - 1.0

      case default
         write(unit=message(1),fmt='(A,I2)') &
            'Unrecognized map projection code: ', proj%code
         write(stdout,*) "da_llxy_wrf.inc"//message(1:1)
         stop 
   end select

end subroutine da_llxy_wrf


subroutine da_xyll(proj, xx, yy, lat, lon)

   
   
   
   

   implicit none

   type(proj_info), intent(in)  :: proj
   real,            intent(in)  :: xx
   real,            intent(in)  :: yy
   real,            intent(out) :: lat
   real,            intent(out) :: lon

   real :: x, y

   if (.NOT.proj%init) then
      write(stdout,*) "da_xyll.inc"// &
           "You have not called map_set for this projection!"
      stop
   end if

   x = xx
   y = yy

   select case (proj%code)
      case (PROJ_LATLON)
         call da_xyll_latlon(x, y, proj, lat, lon)

      case (PROJ_MERC)
         x = xx - proj%knowni + 1.0
         y = yy - proj%knownj + 1.0
         call da_xyll_merc(x, y, proj, lat, lon)

      case (PROJ_PS)
         call da_xyll_ps(x, y, proj, lat, lon)

      case (PROJ_LC)

         x = xx - proj%knowni + 1.0
         y = yy - proj%knownj + 1.0
         call da_xyll_lc(x, y, proj, lat, lon)

      case default
         write(unit=message(1),fmt='(A,I2)') &
            "Unrecognized map projection code: ", proj%code
         write(stdout,*) "da_xyll.inc"//message(1:1)
         stop

   end select

end subroutine da_xyll


subroutine da_xyll_latlon(x, y, proj, lat, lon)

   
   
   

   implicit none

   real, intent(in)             :: x
   real, intent(in)             :: y
   type(proj_info), intent(in)  :: proj
   real, intent(out)            :: lat
   real, intent(out)            :: lon

   real                         :: deltalat
   real                         :: deltalon
   real                         :: latinc
   real                         :: loninc

   
   
   

   latinc = proj%truelat1
   loninc = proj%stdlon

   

   deltalat = (x-1.0)*latinc
   deltalon = (y-1.0)*loninc
   lat = proj%lat1 + deltalat
   lon = proj%lon1 + deltalon

   if ((ABS(lat) > 90.0).OR.(ABS(deltalon) > 360.0)) then
      
      lat = -999.0
      lon = -999.0
   else
      lon = lon + 360.0
      lon = MOD(lon,360.0)
      if (lon > 180.0) lon = lon -360.0
   end if

end subroutine da_xyll_latlon


subroutine da_xyll_lc(x, y, proj, lat, lon)

   
   

   implicit none

   real, intent(in)              :: x        
   real, intent(in)              :: y        
   type(proj_info),intent(in)    :: proj     

                
   real, intent(out)             :: lat      
   real, intent(out)             :: lon      

   real                          :: inew
   real                          :: jnew
   real                          :: r
   real                          :: chi,chi1,chi2
   real                          :: r2
   real                          :: xx
   real                          :: yy

   chi1 = (90.0 - proj%hemi*proj%truelat1)*rad_per_deg
   chi2 = (90.0 - proj%hemi*proj%truelat2)*rad_per_deg

   
   
   if (proj%hemi == -1.0) then 
      inew = -x + 2.0
      jnew = -y + 2.0
   else
      inew = x
      jnew = y
   end if

   
   xx = inew - proj%polei
   yy = proj%polej - jnew
   r2 = (xx*xx + yy*yy)
   r = sqrt(r2)/proj%rebydx
   
   
   if (r2 == 0.0) then
      lat = proj%hemi * 90.0
      lon = proj%stdlon
   else
       
      
      lon = proj%stdlon + deg_per_rad * ATAN2(proj%hemi*xx,yy)/proj%cone
      lon = MOD(lon+360.0, 360.0)

      
      
      
      
        
      if (chi1 == chi2) then
         chi = 2.0*ATAN((r/TAN(chi1))**(1.0/proj%cone) * TAN(chi1*0.5))
      else
         chi = 2.0*ATAN((r*proj%cone/Sin(chi1))**(1.0/proj%cone) * TAN(chi1*0.5)) 
      end if
      lat = (90.0-chi*deg_per_rad)*proj%hemi
   end if

   if (lon > +180.0) lon = lon - 360.0
   if (lon < -180.0) lon = lon + 360.0

end subroutine da_xyll_lc


subroutine da_xyll_merc(x, y, proj, lat, lon)

   
   
   

   implicit none

   real,intent(in)               :: x
   real,intent(in)               :: y    
   type(proj_info),intent(in)    :: proj
   real, intent(out)             :: lat
   real, intent(out)             :: lon 

   lat = 2.0*ATAN(EXP(proj%dlon*(proj%rsw + y-1.0)))*deg_per_rad - 90.0
   lon = (x-1.0)*proj%dlon*deg_per_rad + proj%lon1
   if (lon > 180.0) lon = lon - 360.0
   if (lon < -180.0) lon = lon + 360.0

end subroutine da_xyll_merc


subroutine da_xyll_ps(x, y, proj, lat, lon)

   
   
   

   implicit none

   real, intent(in)                    :: x    
   real, intent(in)                    :: y    
   type (proj_info), intent(in)        :: proj
    
   real, intent(out)                   :: lat     
   real, intent(out)                   :: lon     

   real                                :: reflon
   real                                :: scale_top
   real                                :: xx,yy
   real                                :: gi2, r2
   real                                :: arccos

   
   
   reflon = proj%stdlon + 90.0
   
   
   scale_top = 1.0 + proj%hemi * Sin(proj%truelat1 * rad_per_deg)

   
   xx = x - proj%polei
   yy = (y - proj%polej) * proj%hemi
   r2 = xx**2 + yy**2

   
   if (r2 == 0.0) then 
      lat = proj%hemi * 90.0
      lon = reflon
   else
      gi2 = (proj%rebydx * scale_top)**2.0
      lat = deg_per_rad * proj%hemi * ASin((gi2-r2)/(gi2+r2))
      arccos = ACOS(xx/sqrt(r2))
      if (yy > 0) then
         lon = reflon + deg_per_rad * arccos
      else
         lon = reflon - deg_per_rad * arccos
      end if
   end if
  
   
   if (lon > 180.0) lon = lon - 360.0
   if (lon < -180.0) lon = lon + 360.0

end subroutine da_xyll_ps


subroutine da_set_lc(proj)

   
   
   
   

   implicit none
    
   type(proj_info), intent(inout) :: proj

   real :: arg
   real :: deltalon1
   real :: tl1r
   real :: ctl1r

   
   call da_lc_cone(proj%truelat1, proj%truelat2, proj%cone)
   if (print_detail_map) then
      write(unit=stdout, fmt='(A,F8.6)') 'Computed cone factor: ', proj%cone
   end if
   
   
   deltalon1 = proj%lon1 - proj%stdlon
   if (deltalon1 .gt. +180.0) deltalon1 = deltalon1 - 360.0
   if (deltalon1 .lt. -180.0) deltalon1 = deltalon1 + 360.0

   
   tl1r = proj%truelat1 * rad_per_deg
   ctl1r = COS(tl1r)

   
   proj%rsw = proj%rebydx * ctl1r/proj%cone * &
           (TAN((90.0*proj%hemi-proj%lat1)*rad_per_deg/2.0) / &
            TAN((90.0*proj%hemi-proj%truelat1)*rad_per_deg/2.0))**proj%cone

   
   arg = proj%cone*(deltalon1*rad_per_deg)
   proj%polei = 1.0 - proj%hemi * proj%rsw * Sin(arg)
   proj%polej = 1.0 + proj%rsw * COS(arg)  
   if (print_detail_map) then
      write(unit=stdout,fmt='(A,2F10.3)') 'Computed pole i/j = ', proj%polei, proj%polej
   end if

end subroutine da_set_lc                             


subroutine da_set_ps(proj)

   
   
   
   

   implicit none
 
   type(proj_info), intent(inout)    :: proj

   real :: ala1
   real :: alo1
   real :: reflon
   real :: scale_top

   
   proj%cone = 1.0

   reflon = proj%stdlon + 90.0

   
   scale_top = 1.0 + proj%hemi * Sin(proj%truelat1 * rad_per_deg)

   
   ala1 = proj%lat1 * rad_per_deg
   proj%rsw = proj%rebydx*COS(ala1)*scale_top/(1.0+proj%hemi*Sin(ala1))

   
   alo1 = (proj%lon1 - reflon) * rad_per_deg
   proj%polei = proj%knowni - proj%rsw * COS(alo1)
   proj%polej = proj%knownj - proj%hemi * proj%rsw * Sin(alo1)
   if (print_detail_map) then
      write(unit=stdout,fmt='(A,2F10.1)') 'Computed (I,J) of pole point: ',proj%polei,proj%polej
   end if

end subroutine da_set_ps


subroutine da_map_init(proj)

   
   
   

   implicit none

   type(proj_info), intent(inout)  :: proj

   proj%lat1     = -999.9
   proj%lon1     = -999.9
   proj%dx       = -999.9
   proj%stdlon   = -999.9
   proj%truelat1 = -999.9
   proj%truelat2 = -999.9
   proj%hemi     = 0.0
   proj%cone     = -999.9
   proj%polei    = -999.9
   proj%polej    = -999.9
   proj%rsw      = -999.9
   proj%knowni   = -999.9
   proj%knownj   = -999.9
   proj%latinc   = -999.9
   proj%loninc   = -999.9
   proj%init     = .false.

end subroutine da_map_init


subroutine da_map_set(proj_code,lat1,lon1,knowni,knownj,dx,stdlon,truelat1,truelat2,latinc,loninc,proj)
   
   
   
   
   
   
   
   
   

   implicit none
   
   integer,         intent(in)  :: proj_code
   real,            intent(in)  :: lat1
   real,            intent(in)  :: lon1
   real,            intent(in)  :: dx
   real,            intent(in)  :: stdlon
   real,            intent(in)  :: truelat1
   real,            intent(in)  :: truelat2
   real,            intent(in)  :: knowni , knownj
   real,            intent(in)  :: latinc, loninc 
   type(proj_info), intent(out) :: proj

   
   if (ABS(lat1) > 90.0) then
      write(stdout,*) "da_map_set.inc"// &
         "Latitude of origin corner required as follows: -90N <= lat1 < = 90N"
      stop
   end if
   if (ABS(lon1) > 180.0) then
      write(stdout,*) "da_map_set.inc"// &
           "Longitude of origin required as follows: -180E <= lon1 <= 180W" 
      stop
   end if
   if ((dx .LE. 0.0).AND.(proj_code .NE. PROJ_LATLON)) then
      write(stdout,*)  "da_map_set.inc"// &
           "Require grid spacing (dx) in meters be positive!" 
      stop
   end if
   if ((ABS(stdlon) > 180.0).AND.(proj_code .NE. PROJ_MERC)) then
      write(stdout,*) "da_map_set.inc"// &
           "Need orientation longitude (stdlon) as: -180E <= lon1 <= 180W"  
      stop
   end if
   if (proj%code .NE. PROJ_LATLON .and. ABS(truelat1)>90.0) then
      write(stdout,*) "da_map_set.inc"// &
           "Set true latitude 1 for all projections!" 
      stop
   end if
   
   call da_map_init(proj) 
   proj%code  = proj_code
   proj%lat1 = lat1
   proj%lon1 = lon1
   proj%knowni = knowni
   proj%knownj = knownj
   proj%dx    = dx
   proj%stdlon = stdlon
   proj%truelat1 = truelat1
   proj%truelat2 = truelat2
   proj%latinc   = latinc  
   proj%loninc   = loninc  
   if (proj%code .NE. PROJ_LATLON) then
      proj%dx = dx
      if (truelat1 < 0.0) then
         proj%hemi = -1.0 
      else
         proj%hemi = 1.0
      end if
      proj%rebydx = earth_radius_m / dx
   end if
   pick_proj: select case(proj%code)

      case(PROJ_PS)
         if (print_detail_map) then
            write(unit=stdout,fmt='(A)') 'Setting up POLAR STEREOGRAPHIC map...'
         end if
         call da_set_ps(proj)

      case(PROJ_LC)
         if (print_detail_map) then
            write(unit=stdout,fmt='(A)') 'Setting up LAMBERT CONFORMAL map...'
         end if
         if (ABS(proj%truelat2) > 90.0) then
            if (print_detail_map) then
               write(unit=stdout,fmt='(A)') 'Second true latitude not set, assuming a tangent'
               write(unit=stdout,fmt='(A,F10.3)') 'projection at truelat1: ', proj%truelat1
            end if
            proj%truelat2=proj%truelat1
         end if
         call da_set_lc(proj)
   
      case (PROJ_MERC)
         if (print_detail_map) then
            write(unit=stdout,fmt='(A)') 'Setting up MERCATOR map...'
         end if
         call da_set_merc(proj)
   
      case (PROJ_LATLON)
         if (print_detail_map) then
            write(unit=stdout,fmt='(A)') 'Setting up CYLinDRICAL EQUIDISTANT LATLON map...'
         end if
         
         if (proj%lon1 < 0.0) proj%lon1 = proj%lon1 + 360.0
   
         proj%cone = 1.0                                  
      case default
         write(unit=message(1),fmt='(A,I2)') 'Unknown projection code: ', proj%code
         write(stdout,*) "da_map_set.inc"//message(1:1)
         stop

   end select pick_proj
   proj%init = .true.

end subroutine da_map_set


subroutine da_set_merc(proj)
  
   
   
   

   implicit none

   type(proj_info), intent(inout)       :: proj

   real :: clain

   

   clain = COS(rad_per_deg*proj%truelat1)
   proj%dlon = proj%dx / (earth_radius_m * clain)

   
   

   proj%rsw = 0.0
   if (proj%lat1 .NE. 0.0) then
      proj%rsw = (alog(tan(0.5*((proj%lat1+90.)*rad_per_deg))))/proj%dlon
   end if

end subroutine da_set_merc


subroutine da_lc_cone(truelat1, truelat2, cone)

   
   
   

   implicit none
    
   real, intent(in)             :: truelat1  
   real, intent(in)             :: truelat2  
   real, intent(out)            :: cone

   
   
   
   
   
   if (abs(truelat1-truelat2) > 0.1) then
      cone = alog10(cos(truelat1*rad_per_deg)) - &
             alog10(cos(truelat2*rad_per_deg))
      cone = cone /(alog10(tan((45.0 - abs(truelat1)/2.0) * rad_per_deg)) - &
             alog10(tan((45.0 - abs(truelat2)/2.0) * rad_per_deg)))        
   else
      cone = sin(abs(truelat1)*rad_per_deg)  
   end if

end subroutine da_lc_cone

end module da_verif_tools

