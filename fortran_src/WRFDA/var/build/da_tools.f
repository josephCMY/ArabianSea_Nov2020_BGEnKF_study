












module da_tools
   
   
   
   
   
   use module_bc, only : bdyzone
   use module_dm, only : wrf_dm_sum_real
   use module_domain, only : xb_type, domain

   use da_control, only : pi, gravity, gas_constant, ims, ime, jms,jme, &
      kms,kme,its,ite,jts,jte,kts,kte,ids,ide,stdout, &
      trace_use_dull, trace_use, fg_format_kma_global, coarse_ds, coarse_ix, &
      coarse_jy, fg_format, c2, cone_factor, earth_radius, dsm, &
      map_projection, psi1, pole, start_x, phic, start_y, xlonc, ycntr, &
      obs_qc_pointer, anal_type_verify, fg_format_wrf_arw_regional, &
      fg_format_wrf_nmm_regional, fg_format_wrf_arw_global, fg_format_kma_global, &
      set_omb_rand_fac, fails_error_max, fails_buddy_check, no_buddies, &
      missing_r, x_start_sub_domain, global, myproc, comm, &
      x_end_sub_domain, y_end_sub_domain, def_sub_domain, &
      y_start_sub_domain, start_lat, delt_lat, delt_lon, start_lon, cp, &
      missing_data, surface_correction,print_detail_map, use_rad, stderr, &
      t_kelvin, trace_use_frequent, jds, jde, pptop,ppbot,npres_print, &
      rad_to_deg, deg_to_rad, num_procs, print_detail_obs, psfc_from_slp



   use da_define_structures, only : info_type, field_type, x_type,  &
      model_loc_type, synop_type, bad_info_type, da_gauss_noise, &
      iv_type, y_type, da_random_seed, infa_type
   use da_tools_serial, only : da_array_print
   use da_tracing, only : da_trace_entry, da_trace_exit
   use da_reporting, only : da_error, message, da_warning, da_message
   use da_lapack, only : dsyev
   
   implicit none
   
   include 'mpif.h'

   
   !---------------------------------------------------------------------------
   ! Code copied from SI
   !---------------------------------------------------------------------------

   !dis   
   !dis    Open Source License/Disclaimer, Forecast Systems Laboratory
   !dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
   !dis    
   !dis    This software is distributed under the Open Source Definition,
   !dis    which may be found at http://www.opensource.org/osd.html.
   !dis    
   !dis    In particular, redistribution and use in source and binary forms,
   !dis    with or without modification, are permitted provided that the
   !dis    following conditions are met:
   !dis    
   !dis    - Redistributions of source code must retain this notice, this
   !dis    list of conditions and the following disclaimer.
   !dis    
   !dis    - Redistributions in binary form must provide access to this
   !dis    notice, this list of conditions and the following disclaimer, and
   !dis    the underlying source code.
   !dis    
   !dis    - All modifications to this software must be clearly documented,
   !dis    and are solely the responsibility of the agent making the
   !dis    modifications.
   !dis    
   !dis    - If significant modifications or enhancements are made to this
   !dis    software, the FSL Software Policy Manager
   !dis    (softwaremgr@fsl.noaa.gov) should be notified.
   !dis    
   !dis    THIS SOFTWARE AND ITS doCUMENTATION ARE in THE PUBLIC doMAin
   !dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE unitED STATES
   !dis    GOVERNMENT, ITS inSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
   !dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE useFULNESS
   !dis    OF THE SOFTWARE AND doCUMENTATION FOR ANY Purpose.  THEY ASsumE
   !dis    NO RESPONSIBILITY (1) FOR THE use OF THE SOFTWARE AND
   !dis    doCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO useRS.
   !dis   
   !dis 

   !module map_utils

   ! Module that defines constants, data structures, and
   ! subroutines used to convert grid indices to lat/lon
   ! and vice versa.   
   !
   ! SUPPORTED PROJECTIONS
   ! ---------------------
   ! Cylindrical Lat/Lon (code = PROJ_LATLON)
   ! Mercator (code = PROJ_MERC)
   ! Lambert Conformal (code = PROJ_LC)
   ! Polar Stereographic (code = PROJ_PS)
   !
   ! REMARKS
   ! -------
   ! The routines contained within were adapted from routines
   ! obtained from NCEP's w3 library.  The original NCEP routines were less
   ! flexible (e.g., polar-stereo routines only supported truelat of 60N/60S)
   ! than what we needed, so modifications based on equations in Hoke, Hayes, and
   ! Renninger (AFGWC/TN/79-003) were added to improve the flexibility.  
   ! Additionally, coding was improved to F90 standards and the routines were
   ! combined into this module.  
   !
   ! Assumptions
   ! -----------
   !  Grid Definition:
   !    For mercator, lambert conformal, and polar-stereographic projections,
   !    the routines within assume the following:
   !
   !       1.  Grid is dimensioned (i,j) where i is the East-West direction, 
   !           positive toward the east, and j is the north-south direction, 
   !           positive toward the north.  
   !       2.  Origin is at (1,1) and is located at the southwest corner,
   !           regardless of hemispere.
   !       3.  Grid spacing (dx) is always positive.
   !       4.  Values of true latitudes must be positive for NH domains
   !           and negative for SH domains.
   !
   !     For the latlon projection, the grid origin may be at any of the
   !     corners, and the deltalat and deltalon values can be signed to 
   !     account for this using the following convention:
   !       Origin Location        Deltalat Sign      Deltalon Sign
   !       ---------------        -------------      -------------
   !        SW Corner                  +                   +
   !        NE Corner                  -                   -
   !        NW Corner                  -                   +
   !        SE Corner                  +                   -
   !       
   !  Data Definitions:
   !       1. Any arguments that are a latitude value are expressed in 
   !          degrees north with a valid range of -90 -> 90
   !       2. Any arguments that are a longitude value are expressed in
   !          degrees east with a valid range of -180 -> 180
   !       3. Distances are in meters and are always positive.
   !       4. The standard longitude (stdlon) is defined as the longitude
   !          line which is parallel to the grid's y-axis (j-direction), along
   !          which latitude increases (NOT the absolute value of latitude, but
   !          the actual latitude, such that latitude increases continuously
   !          from the south pole to the north pole) as j increases.  
   !       5. One true latitude value is required for polar-stereographic and
   !          mercator projections, and defines at which latitude the 
   !          grid spacing is true.  For lambert conformal, two true latitude
   !          values must be specified, but may be set equal to each other to
   !          specify a tangent projection instead of a secant projection.  
   !       
   ! USAGE
   ! -----
   ! To use the routines in this module, the calling routines must have the 
   ! following statement at the beginning of its declaration block:
   !   use map_utils
   ! 
   ! The use of the module not only provides access to the necessary routines,
   ! but also defines a structure of type (proj_info) that can be used
   ! to declare a variable of the same type to hold your map projection
   ! information.  It also defines some integer parameters that contain
   ! the projection codes so one only has to use those variable names rather
   ! than remembering the acutal code when using them.  The basic steps are
   ! as follows:
   !  
   !   1.  Ensure the "use map_utils" is in your declarations.
   !   2.  Declare the projection information structure as type(proj_info):
   !         type(proj_info) :: proj
   !   3.  Populate your structure by calling the map_set routine:
   !         call map_set(code,lat1,lon1,knowni,knownj,dx,stdlon,truelat1,truelat2,proj)
   !       where:
   !         code (input) = one of PROJ_LATLON, PROJ_MERC, PROJ_LC, or PROJ_PS
   !         lat1 (input) = Latitude of grid origin point (i,j)=(1,1) 
   !                         (see assumptions!)
   !         lon1 (input) = Longitude of grid origin 
   !         knowni (input) = origin point, x-location
   !         knownj (input) = origin point, y-location
   !         dx (input) = grid spacing in meters (ignored for LATLON projections)
   !         stdlon (input) = Standard longitude for PROJ_PS and PROJ_LC, 
   !               deltalon (see assumptions) for PROJ_LATLON, 
   !               ignored for PROJ_MERC
   !         truelat1 (input) = 1st true latitude for PROJ_PS, PROJ_LC, and
   !                PROJ_MERC, deltalat (see assumptions) for PROJ_LATLON
   !         truelat2 (input) = 2nd true latitude for PROJ_LC, 
   !                ignored for all others.
   !         proj (output) = The structure of type (proj_info) that will be fully 
   !                populated after this call
   !
   !   4.  Now that the proj structure is populated, you may call either
   !       of the following routines:
   !       
   !       da_latlon_to_ij(proj, lat, lon, i, j)
   !       da_xyll(proj, i, j, lat, lon)
   !
   !       It is incumbent upon the calling routine to determine whether or
   !       not the values returned are within your domain's bounds.  All values
   !       of i, j, lat, and lon are real values.
   !
   !
   ! REFERENCES
   ! ----------
   !  Hoke, Hayes, and Renninger, "Map Preojections and Grid Systems for
   !       Meteorological Applications." AFGWC/TN-79/003(Rev), Air Weather
   !       Service, 1985.
   !
   !  NCAR MM5v3 Modeling System, REGRIDDER program, module_first_guess_map.F
   !  NCEP routines w3fb06, w3fb07, w3fb08, w3fb09, w3fb11, w3fb12
   !
   ! HISTORY
   ! -------
   ! 27 Mar 2001 - Original Version
   !               Brent L. Shaw, NOAA/FSL (CSU/CIRA)
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define some private constants

   real, private, parameter    :: deg_per_rad = 180.0/pi
   real, private, parameter    :: rad_per_deg = pi / 180.0

   ! Mean Earth Radius in m.  The value below is consistent
   ! with NCEP's routines and grids.
   ! real, public , parameter    :: earth_radius_m = 6371200.0 ! from Brent
   real, public , parameter    :: earth_radius_m = 6370000.0
   real, public , parameter    :: radians_per_degree = pi / 180.0

   ! Define public parameters
 
   ! Projection codes for proj_info structure:
   integer, public, parameter  :: PROJ_LATLON = 0
   integer, public, parameter  :: PROJ_MERC = 1
   integer, public, parameter  :: PROJ_LC = 3
   integer, public, parameter  :: PROJ_PS = 5

   
  ! Define data structures to define various projections

   type proj_info
      integer          :: code     ! integer code for projection type
      real             :: lat1     ! SW latitude (1,1) in degrees (-90->90N)
      real             :: lon1     ! SW longitude (1,1) in degrees (-180->180E)
      real             :: dx       ! Grid spacing in meters at truelats, used
                                   ! only for ps, lc, and merc projections
      real             :: dlat     ! Lat increment for lat/lon grids
      real             :: dlon     ! Lon increment for lat/lon grids
      real             :: stdlon   ! Longitude parallel to y-axis (-180->180E)
      real             :: truelat1 ! First true latitude (all projections)
      real             :: truelat2 ! Second true lat (LC only)
      real             :: hemi     ! 1 for NH, -1 for SH
      real             :: cone     ! Cone factor for LC projections
      real             :: polei    ! Computed i-location of pole point
      real             :: polej    ! Computed j-location of pole point
      real             :: rsw      ! Computed radius to SW corner
      real             :: rebydx   ! Earth radius divided by dx
      real             :: knowni   ! X-location of known lat/lon
      real             :: knownj   ! Y-location of known lat/lon
      real             :: latinc   ! Latitude increments in degrees 
      real             :: loninc   ! Longitude increments in degrees 
      logical          :: init     ! Flag to indicate if this struct is 
                                 ! ready for use
   end type proj_info

   type(proj_info) :: map_info, map_info_ens



contains

subroutine da_llxy (info, loc, outside, outside_all)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   ! This routine converts (lat, lon) into (x,y) coordinates

   implicit none

   type(info_type),       intent(in)    :: info
   type(model_loc_type),  intent(inout) :: loc
   logical      ,         intent(out)   :: outside      !wrt local domain
   logical, optional,     intent(out)   :: outside_all  !wrt all domains

   ! too many return statments to trace
   ! if (trace_use_frequent) call da_trace_entry("da_llxy")

   outside = .false.
   loc % x   = -1.0
   loc % y   = -1.0
   
   ! get the (x, y) coordinates

   if ( fg_format == fg_format_wrf_arw_regional ) then
      call da_llxy_wrf(map_info, info%lat, info%lon, loc%x, loc%y)
   else if (fg_format == fg_format_wrf_nmm_regional) then
      call da_llxy_rotated_latlon(info%lat, info%lon, map_info, loc%x, loc%y)
   else if (global) then
      call da_llxy_global (info%lat, info%lon, loc%x, loc%y)
   else
      call da_llxy_default (info%lat, info%lon, loc%x, loc%y)
   end if

   call da_togrid (loc%x, its-2, ite+2, loc%i, loc%dx, loc%dxm)!

   call da_togrid (loc%y, jts-2, jte+2, loc%j, loc%dy, loc%dym)

   ! refactor to remove this ugly duplication later
   if (present(outside_all)) then
      outside_all = .false.
      ! Do not check for global options 
      if (.not. global) then 
         if ((int(loc%x) < ids) .or. (int(loc%x) >= ide) .or. &
            (int(loc%y) < jds) .or. (int(loc%y) >= jde)) then
            outside_all = .true. 
            outside = .true. 
            return
         end if
         if (def_sub_domain) then
            if (x_start_sub_domain > loc%x .or. y_start_sub_domain > loc%y .or. &
                x_end_sub_domain   < loc%x .or. y_end_sub_domain   < loc%y) then
               outside_all = .true.
            outside = .true. 
            return
            end if
         end if
      end if
   end if

   if (fg_format == fg_format_kma_global) then
      if ((loc%j < jts-1) .or. (loc%j > jte)) then
         outside = .true.
         return
      end if

      if (loc%j == jde) then
         loc%j = loc%j - 1
         loc%dy  = 1.0
         loc%dym = 0.0
      end if

      return
   end if

   ! Check for edge of domain:

   if ((loc%i < ids) .or. (loc%i >= ide) .or. &
      (loc%j < jds) .or. (loc%j >= jde)) then
      outside     = .true. 
      return
   end if

   ! FIX? hack
   if ((loc%i < its-1) .or. (loc%i > ite) .or. &
      (loc%j < jts-1) .or. (loc%j > jte)) then
   ! if ((loc%i < its-1) .or. (loc%i >= ite) .or. &
   !     (loc%j < jts-1) .or. (loc%j >= jte)) then
      outside = .true.
      return

      if (def_sub_domain) then
         if (x_start_sub_domain > loc%x .or. y_start_sub_domain > loc%y .or. &
             x_end_sub_domain   < loc%x .or. y_end_sub_domain   < loc%y) then
             outside = .true.
         end if
      end if
   end if

   ! if (trace_use_frequent) call da_trace_exit("da_llxy")

end subroutine da_llxy


subroutine da_llxy_new (info, outside, outside_all)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   ! This routine converts (lat, lon) into (x,y) coordinates

   implicit none

   type(infa_type),   intent(inout) :: info
   logical,           intent(inout) :: outside(:,:)      ! wrt local domain
   logical, optional, intent(out)   :: outside_all(:,:)  ! wrt all domains

   if (trace_use) call da_trace_entry("da_llxy_new")

   outside(:,:) = .false.
   info%x(:,:)   = -1.0
   info%y(:,:)   = -1.0
   
   ! get the (x, y) coordinates

   if (fg_format == fg_format_wrf_arw_regional) then
      call da_llxy_wrf_new(map_info, info)
   else if (fg_format == fg_format_wrf_nmm_regional) then
      write(unit=message(1),fmt='(A,I5)') &
         "Needs to be developed for fg_format_nmm_regional = ",fg_format
      call da_error("da_llxy_new.inc",28,message(1:1))
   else if (global) then
      call da_llxy_global_new (info)
   else
      call da_llxy_default_new (info)
   end if

   call da_togrid_new (info%x, its-2, ite+2, info%i, info%dx, info%dxm)
   call da_togrid_new (info%y, jts-2, jte+2, info%j, info%dy, info%dym)

   ! refactor to remove this ugly duplication later
   if (present(outside_all)) then
      outside_all(:,:) = .false.
      ! Do not check for global options 
      if (.not. global) then 
         where ((int(info%x(:,:)) < ids) .or. (int(info%x(:,:)) >= ide) .or. &
            (int(info%y(:,:)) < jds) .or. (int(info%y(:,:)) >= jde))
            outside_all(:,:) = .true. 
            outside(:,:) = .true. 
         end where
         if (def_sub_domain) then
            where (x_start_sub_domain > info%x(:,:) .or. y_start_sub_domain > info%y(:,:) .or. &
                x_end_sub_domain   < info%x(:,:) .or. y_end_sub_domain   < info%y(:,:))
               outside_all(:,:) = .true.
               outside(:,:) = .true. 
            end where
         end if
      end if
   end if

   if (fg_format == fg_format_kma_global) then
      where ((info%j(:,:) < jts-1) .or. (info%j(:,:)  > jte))
         outside(:,:) = .true.
      end where

      where (info%j(:,:) == jde)
         info%j(:,:) = info%j(:,:) - 1
         info%dy(:,:)  = 1.0
         info%dym(:,:) = 0.0
      end where

      return
   end if

   ! Check for edge of domain:

   where ((info%i(:,:) < ids) .or. (info%i(:,:) >= ide) .or. &
      (info%j(:,:) < jds) .or. (info%j(:,:) >= jde))
      outside     = .true. 
   end where

   ! FIX? hack
   where ((info%i(:,:) < its-1) .or. (info%i(:,:) > ite) .or. &
      (info%j(:,:) < jts-1) .or. (info%j(:,:) > jte))
      outside(:,:) = .true.
   end where

   if (def_sub_domain) then
      where (x_start_sub_domain > info%x(:,:) .or. y_start_sub_domain > info%y(:,:) .or. &
             x_end_sub_domain   < info%x(:,:) .or. y_end_sub_domain   < info%y(:,:))
         outside = .true.
      end where 
   end if

   if (trace_use) call da_trace_exit("da_llxy_new")

end subroutine da_llxy_new


subroutine da_llxy_default (xlati,xloni,x,y)

   !----------------------------------------------------------------------------
   ! Purpose:  calculates the (x,y) location (dot) in the mesoscale grids
   ! -------   from latitudes and longitudes
   !
   !           for global domain co-ordinates
   !
   !  input:
   !  -----
   !   xlat:    latitudes
   !   xlon:    longitudes
   !
   ! output:
   ! -----
   !   x:        the coordinate in x (i)-direction.
   !   y:        the coordinate in y (j)-direction.
   !
   !----------------------------------------------------------------------------
   
   implicit none
   
   real, intent(in)  :: xlati, xloni
   real, intent(out) :: x, y

   real              :: dxlon
   real              :: xlat, xlon
   real              :: xx, yy, xc, yc
   real              :: cell, psi0, psx, r, flp
   real              :: centri, centrj
   real              :: ratio
   real              :: bb
   real, parameter   :: conv = 180.0 / pi

   if (trace_use_frequent) call da_trace_entry("da_llxy_default")
   
   xlon = xloni
   xlat = xlati

   xlat = max (xlat, -89.95)
   xlat = min (xlat, +89.95)
   
   dxlon = xlon - xlonc
   if (dxlon >  180) dxlon = dxlon - 360.0
   if (dxlon < -180) dxlon = dxlon + 360.0
   
   if (map_projection == 3) then
      xc = 0.0
      yc = YCNTR

      cell = cos(xlat/conv)/(1.0+sin(xlat/conv))
      yy = -c2*alog(cell)
      xx = c2*dxlon/conv
   else
      psi0 = (pole - phic)/conv
      xc = 0.0

      ! calculate x,y coords. relative to pole

      flp = cone_factor*dxlon/conv
   
      psx = (pole - xlat)/conv
   
      if (map_projection == 2) then
         ! Polar stereographics:
         bb = 2.0*(cos(psi1/2.0)**2)
         yc = -earth_radius*bb*tan(psi0/2.0)
          r = -earth_radius*bb*tan(psx/2.0)
      else
         ! Lambert conformal:
         bb = -earth_radius/cone_factor*sin(psi1)
         yc = bb*(tan(psi0/2.0)/tan(psi1/2.0))**cone_factor
          r = bb*(tan(psx /2.0)/tan(psi1/2.0))**cone_factor
      end if

      if (phic < 0.0) then
         xx = r*sin(flp)
         yy = r*cos(flp)
      else
         xx = -r*sin(flp)
         yy =  r*cos(flp)
      end if

   end if

   ! transform (1,1) to the origin
   ! the location of the center in the coarse domain

   centri = real (coarse_ix + 1)/2.0  
   centrj = real (coarse_jy + 1)/2.0  

   ! the (x,y) coordinates in the coarse domain

   x = (xx - xc)/coarse_ds + centri 
   y = (yy - yc)/coarse_ds + centrj  

   ratio = coarse_ds / dsm

   ! only add 0.5 so that x/y is relative to first cross points:

   x = (x - start_x)*ratio + 0.5
   y = (y - start_y)*ratio + 0.5

   if (trace_use_frequent) call da_trace_exit("da_llxy_default")

end subroutine da_llxy_default


subroutine da_llxy_default_new (info)

   !----------------------------------------------------------------------------
   !
   !                 routine llxy
   !                **************
   !
   !
   ! Purpose:  calculates the (x,y) location (dot) in the mesoscale grids
   ! -------   from latitudes and longitudes
   !
   !           for global domain co-ordinates
   !
   !  input:
   !  -----
   !   xlat:    latitudes
   !   xlon:    longitudes
   !
   ! output:
   ! -----
   !   x:        the coordinate in x (i)-direction.
   !   y:        the coordinate in y (j)-direction.
   !
   !----------------------------------------------------------------------------
   
   implicit none

   type(infa_type), intent(inout) :: info

   real              :: dxlon
   real              :: xlat, xlon
   real              :: xx, yy, xc, yc
   real              :: cell, psi0, psx, r, flp
   real              :: centri, centrj
   real              :: ratio
   real              :: bb
   real, parameter   :: conv = 180.0 / pi

   integer :: n

   if (trace_use) call da_trace_entry("da_llxy_default_new")

   ! Slow, but I can't be arsed to do all the temporary arrays

   do n=1,ubound(info%lat,2)
      xlon = info%lon(1,n)
      xlat = info%lat(1,n)

      xlat = max (xlat, -89.95)
      xlat = min (xlat, +89.95)

      dxlon = xlon - xlonc
      if (dxlon >  180) dxlon = dxlon - 360.0
      if (dxlon < -180) dxlon = dxlon + 360.0

      if (map_projection == 3) then
         xc = 0.0
         yc = YCNTR

         cell = cos(xlat/conv)/(1.0+sin(xlat/conv))
         yy = -c2*alog(cell)
         xx = c2*dxlon/conv
      else
         psi0 = (pole - phic)/conv
         xc = 0.0

         ! calculate x,y coords. relative to pole

         flp = cone_factor*dxlon/conv

         psx = (pole - xlat)/conv

         if (map_projection == 2) then
            ! Polar stereographics:
            bb = 2.0*(cos(psi1/2.0)**2)
            yc = -earth_radius*bb*tan(psi0/2.0)
             r = -earth_radius*bb*tan(psx/2.0)
         else
            ! Lambert conformal:
            bb = -earth_radius/cone_factor*sin(psi1)
            yc = bb*(tan(psi0/2.0)/tan(psi1/2.0))**cone_factor
             r = bb*(tan(psx /2.0)/tan(psi1/2.0))**cone_factor
         end if

         if (phic < 0.0) then
            xx = r*sin(flp)
            yy = r*cos(flp)
         else
            xx = -r*sin(flp)
            yy =  r*cos(flp)
         end if

      end if

      ! transform (1,1) to the origin
      ! the location of the center in the coarse domain

      centri = real (coarse_ix + 1)/2.0  
      centrj = real (coarse_jy + 1)/2.0  

      ! the (x,y) coordinates in the coarse domain

      info%x(1,n) = (xx - xc)/coarse_ds + centri 
      info%y(1,n) = (yy - yc)/coarse_ds + centrj  

      ratio = coarse_ds / dsm

      ! only add 0.5 so that x/y is relative to first cross points:

      info%x(:,n) = (info%x(1,n) - start_x)*ratio + 0.5
      info%y(:,n) = (info%y(1,n) - start_y)*ratio + 0.5
   end do

   if (trace_use) call da_trace_exit("da_llxy_default_new")

end subroutine da_llxy_default_new


subroutine da_llxy_kma_global(lat,lon,x,y)

   !----------------------------------------------------------------------------
   ! Purpose:  calculates the(x,y) location(dot) in the global grids
   !           from latitudes and longitudes
   !----------------------------------------------------------------------------
   
   implicit none
   
   real, intent(in)  :: lat, lon
   real, intent(out) :: x, y

   real              :: xlat, xlon

   if (trace_use_frequent) call da_trace_entry("da_llxy_kma_global")
   
   xlat = lat - start_lat
   xlon = lon - start_lon

   if (xlat < 0.0) xlat = xlat + 180.0
   if (xlon < 0.0) xlon = xlon + 360.0

   x = start_x + xlon/delt_lon
   y = start_y + xlat/delt_lat

   if (trace_use_frequent) call da_trace_exit("da_llxy_kma_global")
   
end subroutine da_llxy_kma_global


subroutine da_llxy_kma_global_new(info)

   !----------------------------------------------------------------------------
   ! Purpose:  calculates the(x,y) location(dot) in the global grids
   !           from latitudes and longitudes
   !----------------------------------------------------------------------------
   
   implicit none

   type(infa_type), intent(inout) :: info

   real    :: xlat, xlon
   integer :: n

   if (trace_use) call da_trace_entry("da_llxy_kma_global_new")

!FAST

!   where (lat(:,:) - start_lat < 0)
!      y(:,:) = start_y + (lat(:,:) - start_lat+180.0)/delt_lat
!   else
!      y(:,:) = start_y + (lat(:,:) - start_lat)/delt_lat
!   end where

!   where (lon(:,:) - start_lon < 0.)
!      x(:,:) = start_x + (lon(:,:) - start_lon+360.0)/delt_lon
!   else
!      x(:,:) = start_x + (lon(:,:) - start_lon)/delt_lon
!   end where

! SLOW

   do n=lbound(info%lat,2),ubound(info%lat,2)
      xlat = info%lat(1,n) - start_lat
      xlon = info%lon(1,n) - start_lon
      if (xlat < 0.0) xlat = xlat + 180.0
      if (xlon < 0.0) xlon = xlon + 360.0
      info%x(:,n) = start_x + xlon/delt_lon
      info%y(:,n) = start_y + xlat/delt_lat
   end do

   if (trace_use) call da_trace_exit("da_llxy_kma_global_new")
   
end subroutine da_llxy_kma_global_new


subroutine da_llxy_global(lat,lon,x,y)

   !----------------------------------------------------------------------------
   ! Purpose:  calculates the(x,y) location(dot) in the global grids
   !           from latitudes and longitudes
   !----------------------------------------------------------------------------
   
   implicit none
   
   real, intent(in)  :: lat, lon
   real, intent(out) :: x, y

   real              :: xlat, xlon

   if (trace_use_frequent) call da_trace_entry("da_llxy_global")
   xlat = lat - start_lat
   xlon = lon - start_lon
   if (xlat < 0.0) xlat = xlat + 180.0
   if (xlon < 0.0) xlon = xlon + 360.0

   x = start_x + xlon/delt_lon
   y = start_y + xlat/delt_lat
   if((fg_format == fg_format_wrf_arw_global) .and. (lat.le.start_lat)) y = 1.0 

   if (trace_use_frequent) call da_trace_exit("da_llxy_global")
   
end subroutine da_llxy_global
subroutine da_llxy_global_new(info)

   !----------------------------------------------------------------------------
   ! Purpose:  calculates the(x,y) location(dot) in the global grids
   !           from latitudes and longitudes
   !----------------------------------------------------------------------------
   
   implicit none

   type(infa_type), intent(inout) :: info

   real    :: xlat, xlon
   integer :: n

   if (trace_use) call da_trace_entry("da_llxy_global_new")

!FAST

!   where (lat(:,:) - start_lat < 0)
!      y(:,:) = start_y + (lat(:,:) - start_lat+180.0)/delt_lat
!   else
!      y(:,:) = start_y + (lat(:,:) - start_lat)/delt_lat
!   end where

!   where (lon(:,:) - start_lon < 0.)
!      x(:,:) = start_x + (lon(:,:) - start_lon+360.0)/delt_lon
!   else
!      x(:,:) = start_x + (lon(:,:) - start_lon)/delt_lon
!   end where

! SLOW

   do n=lbound(info%lat,2),ubound(info%lat,2)
      xlat = info%lat(1,n) - start_lat
      xlon = info%lon(1,n) - start_lon
      if (xlat < 0.0) xlat = xlat + 180.0
      if (xlon < 0.0) xlon = xlon + 360.0
      info%x(:,n) = start_x + xlon/delt_lon
      info%y(:,n) = start_y + xlat/delt_lat
      if((fg_format == fg_format_wrf_arw_global) .and. (info%lat(1,n).le.start_lat)) info%y(:,n) = 1.0                       
   end do

   if (trace_use) call da_trace_exit("da_llxy_global_new")
   
end subroutine da_llxy_global_new


subroutine da_llxy_rotated_latlon(lat,lon, proj, x, y)                          

   !----------------------------------------------------------------------- 
   ! Purpose: Compute the x/y location of a lat/lon on a rotated LATLON grid.
   ! Author :  Syed RH Rizvi,     MMM/NCAR
   !           06/01/2008
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)             :: lat
   real, intent(in)             :: lon
   type(proj_info), intent(in)  :: proj
   real, intent(out)            :: x
   real, intent(out)            :: y

   real                         :: rot_lat, rot_lon, deltalat,deltalon, lon360,latinc,loninc 
   real                         :: xlat, xlon, cen_lat, cen_lon


   if (trace_use_frequent) call da_trace_entry("da_llxy_rotated_latlon")
   ! To account for issues around the dateline, convert the incoming
   ! longitudes to be 0->360.0
   if (lon < 0) then 
      lon360 = lon + 360.0 
   else 
      lon360 = lon
   end if    

   xlat = deg_to_rad*lat
   xlon = deg_to_rad*lon360
   cen_lat = deg_to_rad*proj%lat1
   cen_lon = deg_to_rad*proj%lon1
   if (cen_lon < 0.) cen_lon = cen_lon + 360.
  
   latinc = proj%latinc 
   loninc = proj%loninc  

   rot_lon = rad_to_deg*atan( cos(xlat) * sin(xlon-cen_lon)/ &
             (cos(cen_lat)*cos(xlat)*cos(xlon-cen_lon) + sin(cen_lat)*sin(xlat)))
   rot_lat = rad_to_deg*asin(  cos(cen_lat)*sin(xlat) - sin(cen_lat)*cos(xlat)*cos(xlon-cen_lon))


   deltalat = rot_lat 
   deltalon = rot_lon 

    
   ! Compute x/y
   x = proj%knowni + deltalon/loninc + 1.0
   y = proj%knownj + deltalat/latinc + 1.0
   
   if (trace_use_frequent) call da_trace_exit("da_llxy_rotated_latlon")

end subroutine da_llxy_rotated_latlon
subroutine da_llxy_latlon(lat, lon, proj, x, y)

   !----------------------------------------------------------------------- 
   ! Purpose: Compute the x/y location of a lat/lon on a LATLON 
   !          (cylindrical equidistant) grid.
   !-----------------------------------------------------------------------

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

   if (trace_use_frequent) call da_trace_entry("da_llxy_latlon")

   ! To account for issues around the dateline, convert the incoming
   ! longitudes to be 0->360.0
   if (lon < 0) then
      lon360 = lon + 360.0
   else
      lon360 = lon
   end if

   deltalat = lat - proj%lat1
   deltalon = lon360 - proj%lon1

   !For cylindrical equidistant, dx == dy
   loninc = proj%dx*360.0/(2.0*EARTH_RADIUS_M*PI)
   latinc = proj%dx*360.0/(2.0*EARTH_RADIUS_M*PI)

   ! Compute x/y
   x = deltalon/loninc
   y = deltalat/latinc

   x = x + proj%knowni
   y = y + proj%knownj


   if (trace_use_frequent) call da_trace_exit("da_llxy_latlon")

end subroutine da_llxy_latlon


subroutine da_llxy_latlon_new(proj, info)

   !----------------------------------------------------------------------- 
   ! Purpose: Compute the x/y location of a lat/lon on a LATLON grid.
   !-----------------------------------------------------------------------

   implicit none

   type(proj_info), intent(in)    :: proj
   type(infa_type), intent(inout) :: info

   integer :: n

   if (trace_use) call da_trace_entry("da_llxy_latlon_new")

   ! Extract the latitude and longitude increments for this grid
   ! (e.g., 2.5 deg for NCEP reanalysis) from the proj structure, where
   ! loninc is saved in the stdlon tag and latinc is saved in truelat1

   ! Compute the difference between the input lat/lon and the origin lat/lon

   info%y(1,:) = (info%lat(1,:) - proj%lat1)/proj%truelat1 + 1.0

   ! To account for issues around the dateline, convert the incoming
   ! longitudes to be 0->360.0
   where (info%lon(1,:) < 0)
      info%x(1,:) = (info%lon(1,:) + 360.0  - proj%lon1)/proj%stdlon + 1.0
   elsewhere
      info%x(1,:) = (info%lon(1,:) - proj%lon1)/proj%stdlon + 1.0
   end where

   do n=1,ubound(info%x,2)
      info%x(:,n) = info%x(1,n)
      info%y(:,n) = info%y(1,n)
   end do

   if (trace_use) call da_trace_exit("da_llxy_latlon_new")

end subroutine da_llxy_latlon_new


subroutine da_llxy_lc(lat, lon, proj, x, y)

   !-----------------------------------------------------------------------
   ! Purpose: compute the geographical latitude and longitude values
   ! to the cartesian x/y on a Lambert Conformal projection.
   !-----------------------------------------------------------------------
    
   implicit none

   real, intent(in)              :: lat      ! Latitude (-90->90 deg N)
   real, intent(in)              :: lon      ! Longitude (-180->180 E)
   type(proj_info),intent(in)    :: proj     ! Projection info structure

   real, intent(out)             :: x        ! Cartesian X coordinate
   real, intent(out)             :: y        ! Cartesian Y coordinate

   real                          :: arg
   real                          :: deltalon
   real                          :: tl1r
   real                          :: rm
   real                          :: ctl1r

   if (trace_use_dull) call da_trace_entry("da_llxy_lc")
    
   ! Compute deltalon between known longitude and standard lon and ensure
   ! it is not in the cut zone
   deltalon = lon - proj%stdlon
   if (deltalon > +180.0) deltalon = deltalon - 360.0
   if (deltalon < -180.0) deltalon = deltalon + 360.0
    
   ! Convert truelat1 to radian and compute COS for later use
   tl1r = proj%truelat1 * rad_per_deg
   ctl1r = COS(tl1r)     
   
   ! Radius to desired point
   rm = proj%rebydx * ctl1r/proj%cone * &
       (TAN((90.0*proj%hemi-lat)*rad_per_deg/2.0) / &
        TAN((90.0*proj%hemi-proj%truelat1)*rad_per_deg/2.0))**proj%cone

   arg = proj%cone*(deltalon*rad_per_deg)
   x = proj%polei + proj%hemi * rm * Sin(arg)
   y = proj%polej - rm * COS(arg)

   ! Finally, if we are in the southern hemisphere, flip the i/j
   ! values to a coordinate system where (1,1) is the SW corner
   ! (what we assume) which is different than the original NCEP
   ! algorithms which used the NE corner as the origin in the 
   ! southern hemisphere (left-hand vs. right-hand coordinate?)
   if (proj%hemi == -1.0) then
      x = 2.0 - x  
      y = 2.0 - y
   end if

   if (trace_use_dull) call da_trace_exit("da_llxy_lc")

end subroutine da_llxy_lc


subroutine da_llxy_lc_new(proj, info)

   !-----------------------------------------------------------------------
   ! Purpose: compute the geographical latitude and longitude values
   ! to the cartesian x/y on a Lambert Conformal projection.
   !-----------------------------------------------------------------------
    
   implicit none

   type(proj_info), intent(in)    :: proj     ! Projection info structure
   type(infa_type), intent(inout) :: info


   real    :: tl1r
!   real    :: temp1
!   real    :: temp2
   real    :: ctl1r
   real    :: deltalon, rm, arg
   integer :: n

   if (trace_use) call da_trace_entry("da_llxy_lc_new")

! FAST
    
   ! Convert truelat1 to radian and compute COS for later use
!   tl1r = proj%truelat1 * rad_per_deg
!   ctl1r = COS(tl1r)    
!   temp1 = TAN((90.0*proj%hemi-proj%truelat1)*rad_per_deg/2.0)
!   temp2 = proj%rebydx * ctl1r/proj%cone
    
   ! Compute deltalon between known longitude and standard lon and ensure
   ! it is not in the cut zone

!   where (lon - proj%stdlon > +180.0)
!      x = proj%polei + proj%hemi * (temp2 * (TAN((90.0*proj%hemi-lat)*rad_per_deg/2.0) / temp1)**proj%cone &
!         * SIN(proj%cone*((lon - proj%stdlon-360.0)*rad_per_deg))
!      y = proj%polej - (temp2 * (TAN((90.0*proj%hemi-lat)*rad_per_deg/2.0) / temp1)**proj%cone &
!         * COS(proj%cone*((lon - proj%stdlon-360.0)*rad_per_deg))
!   elsewhere  (lon - proj%stdlon - -180.0)
!      x = proj%polei + proj%hemi * (temp2 * (TAN((90.0*proj%hemi-lat)*rad_per_deg/2.0) / temp1) &
!         * SIN(proj%cone*((lon - proj%stdlon+360.0)*rad_per_deg))
!      y = proj%polej - (temp2 * (TAN((90.0*proj%hemi-lat)*rad_per_deg/2.0) / temp1)**proj%cone &
!         * COS(proj%cone*((lon - proj%stdlon+360.0)*rad_per_deg))
!   else
!      x = proj%polei + proj%hemi * (temp2 * (TAN((90.0*proj%hemi-lat)*rad_per_deg/2.0) / temp1) &
!         * SIN(proj%cone*((lon - proj%stdlon)*rad_per_deg))
!      y = proj%polej - (temp2 * (TAN((90.0*proj%hemi-lat)*rad_per_deg/2.0) / temp1)**proj%cone &
!         * COS(proj%cone*((lon - proj%stdlon)*rad_per_deg))
!   end where

   ! Finally, if we are in the southern hemisphere, flip the i/j
   ! values to a coordinate system where (1,1) is the SW corner
   ! (what we assume) which is different than the original NCEP
   ! algorithms which used the NE corner as the origin in the 
   ! southern hemisphere (left-hand vs. right-hand coordinate?)
!   if (proj%hemi == -1.0) then
!      x(:,:) = 2.0 - i(:,:)
!      y(:,:) = 2.0 - j(:,:)
!   end if

! SLOW

   do n=lbound(info%lat,2), ubound(info%lat,2)
      ! Compute deltalon between known longitude and standard lon and ensure
      ! it is not in the cut zone
      deltalon = info%lon(1,n) - proj%stdlon
      if (deltalon > +180.0) deltalon = deltalon - 360.0
      if (deltalon < -180.0) deltalon = deltalon + 360.0

      ! Convert truelat1 to radian and compute COS for later use
      tl1r = proj%truelat1 * rad_per_deg
      ctl1r = COS(tl1r)     

      ! Radius to desired point
      rm = proj%rebydx * ctl1r/proj%cone * &
          (TAN((90.0*proj%hemi-info%lat(1,n))*rad_per_deg/2.0) / &
           TAN((90.0*proj%hemi-proj%truelat1)*rad_per_deg/2.0))**proj%cone

      arg = proj%cone*(deltalon*rad_per_deg)
      info%x(:,n) = proj%polei + proj%hemi * rm * Sin(arg)
      info%y(:,n) = proj%polej - rm * COS(arg)

      ! Finally, if we are in the southern hemisphere, flip the i/j
      ! values to a coordinate system where (1,1) is the SW corner
      ! (what we assume) which is different than the original NCEP
      ! algorithms which used the NE corner as the origin in the 
      ! southern hemisphere (left-hand vs. right-hand coordinate?)
      if (proj%hemi == -1.0) then
         info%x(:,n) = 2.0 - info%x(:,n)  
	 info%y(:,n) = 2.0 - info%y(:,n)
      end if
   end do

   if (trace_use) call da_trace_exit("da_llxy_lc_new")

end subroutine da_llxy_lc_new


subroutine da_llxy_merc(lat, lon, proj, x, y)

   !-----------------------------------------------------------------------
   ! Purpose: Compute x,y coordinate from lat lon for mercator projection
   !-----------------------------------------------------------------------
  
   implicit none

   real, intent(in)              :: lat
   real, intent(in)              :: lon
   type(proj_info),intent(in)    :: proj
   real,intent(out)              :: x
   real,intent(out)              :: y
   real                          :: deltalon

   if (trace_use_frequent) call da_trace_entry("da_llxy_merc")

   deltalon = lon - proj%lon1
   if (deltalon < -180.0) deltalon = deltalon + 360.0
   if (deltalon > 180.0) deltalon = deltalon - 360.0
   x = 1.0 + (deltalon/(proj%dlon*deg_per_rad))
   y = 1.0 + (ALOG(TAN(0.5*((lat + 90.0) * rad_per_deg)))) / &
           proj%dlon - proj%rsw

   if (trace_use_frequent) call da_trace_exit("da_llxy_merc")

end subroutine da_llxy_merc


subroutine da_llxy_merc_new(proj, info)

   !-----------------------------------------------------------------------
   ! Purpose: Compute x/y coordinate from lat lon for mercator projection
   !-----------------------------------------------------------------------
  
   implicit none

   type(proj_info), intent(in)    :: proj
   type(infa_type), intent(inout) :: info

   real :: deltalon
   integer :: n

if (trace_use) call da_trace_entry("da_llxy_merc_new")

! FAST

!   where (lon(:,:) - proj%lon1 < -180.0)
!      x(:,:) = 1.0 + (lon(:,:) - proj%lon1 + 360.0)/(proj%dlon*deg_per_rad))
!   elsewhere (lon(:,:) - proj%lon1 > 180.0)
!      x(:,:) = 1.0 + (lon(:,:) - proj%lon1 - 360.0)/(proj%dlon*deg_per_rad))
!   else
!      x(:,:) = 1.0 + (lon(:,:) - proj%lon1)/(proj%dlon*deg_per_rad))
!   end where

!   y(:,:) = 1.0 + (ALOG(TAN(0.5*((lat(:,: + 90.0) * rad_per_deg)))) / proj%dlon - proj%rsw

! SLOW

   do n=lbound(info%lat,2),ubound(info%lat,2)
      deltalon = info%lon(1,n) - proj%lon1
      if (deltalon < -180.0) deltalon = deltalon + 360.0
      if (deltalon > 180.0) deltalon = deltalon - 360.0
      info%x(:,n) = 1.0 + (deltalon/(proj%dlon*deg_per_rad))
      info%y(:,n) = 1.0 + (ALOG(TAN(0.5*((info%lat(1,n) + 90.0) * rad_per_deg)))) / proj%dlon - proj%rsw
   end do

if (trace_use) call da_trace_exit("da_llxy_merc_new")


end subroutine da_llxy_merc_new


subroutine da_llxy_ps(lat,lon,proj,x,y)

   !-----------------------------------------------------------------------
   ! Purpose: Given latitude (-90 to 90), longitude (-180 to 180), and the
   ! standard polar-stereographic projection information via the 
   ! public proj structure, this routine returns the x/y indices which
   ! if within the domain range from 1->nx and 1->ny, respectively.
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)               :: lat
   real, intent(in)               :: lon
   type(proj_info),intent(in)     :: proj

   real, intent(out)              :: x !(x-index)
   real, intent(out)              :: y !(y-index)
   
   real                           :: reflon
   real                           :: scale_top
   real                           :: ala
   real                           :: alo
   real                           :: rm

   if (trace_use_frequent) call da_trace_entry("da_llxy_ps")

   reflon = proj%stdlon + 90.0
   
   ! Compute numerator term of map scale factor

   scale_top = 1.0 + proj%hemi * Sin(proj%truelat1 * rad_per_deg)

   ! Find radius to desired point
   ala = lat * rad_per_deg
   rm = proj%rebydx * COS(ala) * scale_top/(1.0 + proj%hemi *Sin(ala))
   alo = (lon - reflon) * rad_per_deg
   x = proj%polei + rm * COS(alo)
   y = proj%polej + proj%hemi * rm * Sin(alo)

   if (trace_use_frequent) call da_trace_exit("da_llxy_ps")
 
end subroutine da_llxy_ps


subroutine da_llxy_ps_new(proj,info)

   !-----------------------------------------------------------------------
   ! Purpose: Given latitude (-90 to 90), longitude (-180 to 180), and the
   ! standard polar-stereographic projection information via the 
   ! public proj structure, this routine returns the x/y indices which
   ! if within the domain range from 1->nx and 1->ny, respectively.
   !-----------------------------------------------------------------------

   implicit none

   type(proj_info), intent(in)    :: proj
   type(infa_type), intent(inout) :: info
   
   real    :: reflon
   integer :: n
   real    :: scale_top, ala, rm, alo

   if (trace_use) call da_trace_entry("da_llxy_ps_new")

   reflon = proj%stdlon + 90.0

   ! FAST

!   x(:,:) = proj%polei + proj%rebydx * COS(lat(:,:) * rad_per_deg) &
!      * (1.0 + proj%hemi * SIN(proj%truelat1 * rad_per_deg) &
!      / (1.0 + proj%hemi *SIN(lat(:,:) * rad_per_deg)) &
!      * COS((lon(:,:) - reflon) * rad_per_deg)

!   y(:,:) = proj%polej + proj%hemi * proj%rebydx * COS(lat(:,:) * rad_per_deg) &
!      * (1.0 + proj%hemi * SIN(proj%truelat1 * rad_per_deg) &
!      / (1.0 + proj%hemi *SIN(lat(:,:) * rad_per_deg)) &
!      * SIN((lon(:,:) - reflon) * rad_per_deg)

! SLOW

   do n=lbound(info%lat,2),ubound(info%lat,2)
      ! compute numerator term of map scale factor

      scale_top = 1.0 + proj%hemi * sin(proj%truelat1 * rad_per_deg)

      ! find radius to desired point
      ala = info%lat(1,n) * rad_per_deg
      rm = proj%rebydx * cos(ala) * scale_top/(1.0 + proj%hemi *sin(ala))
      alo = (info%lon(1,n) - reflon) * rad_per_deg
      info%x(:,n) = proj%polei + rm * cos(alo)
      info%y(:,n) = proj%polej + proj%hemi * rm * sin(alo)
   end do

   if (trace_use) call da_trace_exit("da_llxy_ps_new")
 
end subroutine da_llxy_ps_new


subroutine da_llxy_wrf(proj, lat, lon, x, y)

   !-----------------------------------------------------------------------
   ! Purpose: Converts input lat/lon values to the cartesian (x, y) value
   ! for the given projection. 
   !-----------------------------------------------------------------------

   implicit none

   type(proj_info), intent(in)  :: proj
   real,            intent(in)  :: lat
   real,            intent(in)  :: lon
   real,            intent(out) :: x
   real,            intent(out) :: y

   if (trace_use_frequent) call da_trace_entry("da_llxy_wrf")

   if (.NOT.proj%init) then
      call da_error("da_llxy_wrf.inc",19, &
        (/"You have not called map_set for this projection!"/))
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
         call da_error("da_llxy_wrf.inc",44,message(1:1))
   end select

   if (trace_use_frequent) call da_trace_exit("da_llxy_wrf")

end subroutine da_llxy_wrf


subroutine da_llxy_wrf_new(proj, info)

   !-----------------------------------------------------------------------
   ! Purpose: Converts input lat/lon values to the cartesian (x,y) value
   ! for the given projection. 
   !-----------------------------------------------------------------------

   implicit none

   type(proj_info), intent(in)    :: proj
   type(infa_type), intent(inout) :: info

   if (trace_use) call da_trace_entry("da_llxy_wrf_new")

   if (.NOT.proj%init) then
      call da_error("da_llxy_wrf_new.inc",16, &
        (/"You have not called map_set for this projection!"/))
   end if

   select case(proj%code)
 
      case(PROJ_LATLON)
         call da_llxy_latlon_new(proj,info)

      case(PROJ_MERC)
         call da_llxy_merc_new(proj,info)
         info%x(:,:) = info%x(:,:) + proj%knowni - 1.0
         info%y(:,:) = info%y(:,:) + proj%knownj - 1.0

      case(PROJ_PS)
         call da_llxy_ps_new(proj,info)
      
      case(PROJ_LC)
         call da_llxy_lc_new(proj,info)
         info%x(:,:) = info%x(:,:) + proj%knowni - 1.0
         info%y(:,:) = info%y(:,:) + proj%knownj - 1.0

      case default
         write(unit=message(1),fmt='(A,I2)') &
            'Unrecognized map projection code: ', proj%code
         call da_error("da_llxy_wrf_new.inc",41,message(1:1))
   end select

   if (trace_use) call da_trace_exit("da_llxy_wrf_new")

end subroutine da_llxy_wrf_new


subroutine da_xyll(proj, xx, yy, lat, lon)

   !-----------------------------------------------------------------------
   ! Purpose: Computes geographical latitude and longitude for a given (i,j) 
   ! point in a grid with a projection of proj
   !-----------------------------------------------------------------------

   implicit none

   type(proj_info), intent(in)  :: proj
   real,            intent(in)  :: xx
   real,            intent(in)  :: yy
   real,            intent(out) :: lat
   real,            intent(out) :: lon

   real :: x, y

   if (trace_use) call da_trace_entry("da_xyll")

   if (.NOT.proj%init) then
      call da_error("da_xyll.inc",21, &
         (/"You have not called map_set for this projection!"/))
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
         call da_error("da_xyll.inc",49,message(1:1))

   end select

   if (trace_use) call da_trace_exit("da_xyll")

end subroutine da_xyll


subroutine da_xyll_default(XX,YY,XLAT,XLON)

   !-----------------------------------------------------------------------
   ! Purpose: calculates the latitudes and longitudes from the
   !                  (x,y) location (dot) in the mesoscale grids.
   ! on entry     :   
   ! x            : the coordinate in x (j)-direction.
   ! y            : the coordinate in y (i)-direction.
   !
   ! on exit      :                      
   ! xlat         : latitudes 
   ! xlon         : longitudes 
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)  :: XX, YY
   real, intent(out) :: XLAT,XLON
        
   real              :: flp, flpp, r, cell, cel1, cel2
   real              :: rcone_factor, psx, conv
   real              :: cntri, cntrj, x, y, xcntr

   if (trace_use) call da_trace_entry("da_xyll_default")
   
   conv = 180.0 / pi

   ! seperate treatment for global fields
   if (fg_format == fg_format_kma_global) then
      xlat = yy * 180.0 /(coarse_jy-1)  -  90.0    
      xlon = xx * 360.0 /(coarse_ix-1)  - 180.0    
      return
   end if

   cntri = real(coarse_ix+1)/2.0
   cntrj = real(coarse_jy+1)/2.0

   xcntr = 0.0

   ! calculate x and y positions of grid

   x = (xx - 0.5)*dsm/coarse_ds + start_x
   y = (yy - 0.5)*dsm/coarse_ds + start_y
   x = xcntr + (x-cntri)*coarse_ds
   y = ycntr + (y-cntrj)*coarse_ds

   ! now calculate lat and lon of this point

   if (map_projection.ne.3) then
      if(y.eq.0.0) then      
         if (x.ge.0.0) flp =  90.0/conv 
         if (x.lt.0.0) flp = -90.0/conv
      else
         if (phic.lt.0.0)then
            flp = atan2(x,y)   
         else
            flp = atan2(x,-y) 
         end if
      end if 
      flpp = (flp/cone_factor)*conv+xlonc  
      if (flpp.lt.-180.0) flpp = flpp + 360.0    
      if (flpp.gt.180.0)  flpp = flpp - 360.0  
      xlon = flpp 
      ! now solve for latitude
      r = sqrt(x*x+y*y)  
      if (phic.lt.0.0) r = -r  
      if (map_projection.eq.1) then   
         cell = (r*cone_factor)/(earth_radius*sin(psi1))    
         rcone_factor  = 1.0/cone_factor   
         cel1 = tan(psi1/2.0)*(cell)**rcone_factor    
      end if 
      if (map_projection.eq.2) then
         cell = r/earth_radius        
         cel1 = cell/(1.0+cos(psi1))  
      end if 
      cel2 = atan(cel1)    
      psx  = 2.0*cel2*conv
      xlat = pole-psx 
   end if   
   ! calculations for mercator lat,lon    
   if (map_projection.eq.3) then   
      xlon = xlonc + ((x-xcntr)/c2)*conv 
      if (xlon.lt.-180.0) xlon = xlon + 360.0
      if (xlon.gt.180.0)  xlon = xlon - 360.0
      cell = exp(y/c2)  
      xlat = 2.0*(conv*atan(cell))-90.0 
   end if  

   if (trace_use) call da_trace_exit("da_xyll_default")

end subroutine da_xyll_default


subroutine da_xyll_latlon(x, y, proj, lat, lon)

   !-----------------------------------------------------------------------
   ! Purpose: Compute the lat/lon location of an x/y on a LATLON grid.
   !-----------------------------------------------------------------------

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
   real                         :: lon360

   if (trace_use_frequent) call da_trace_entry("da_xyll_latlon")

   loninc = proj%dx*360.0/(2.0*EARTH_RADIUS_M*PI)
   latinc = proj%dx*360.0/(2.0*EARTH_RADIUS_M*PI)

   deltalon = (x - proj%knowni) * loninc
   deltalat = (y - proj%knownj) * latinc

   lon360 = deltalon + proj%lon1
   lat    = deltalat + proj%lat1

   if ((ABS(lat) > 90.0).OR.(ABS(deltalon) > 360.0)) then
      ! Off the earth for this grid; THIS SHOULD NEVER HAPPEN
      lat = -999.0
      lon = -999.0
   else
      lon = MOD(lon360,360.0)
      if (lon > 180.0) lon = lon - 360.0
   end if

   if (trace_use_frequent) call da_trace_exit("da_xyll_latlon")

end subroutine da_xyll_latlon


subroutine da_xyll_lc(x, y, proj, lat, lon)

   ! subroutine da_to convert from the (x,y) cartesian coordinate to the 
   ! geographical latitude and longitude for a Lambert Conformal projection.

   implicit none

   real, intent(in)              :: x        ! Cartesian X coordinate
   real, intent(in)              :: y        ! Cartesian Y coordinate
   type(proj_info),intent(in)    :: proj     ! Projection info structure

                
   real, intent(out)             :: lat      ! Latitude (-90->90 deg N)
   real, intent(out)             :: lon      ! Longitude (-180->180 E)

   real                          :: inew
   real                          :: jnew
   real                          :: r
   real                          :: chi,chi1,chi2
   real                          :: r2
   real                          :: xx
   real                          :: yy

   if (trace_use_dull) call da_trace_entry("da_xyll_lc")


   chi1 = (90.0 - proj%hemi*proj%truelat1)*rad_per_deg
   chi2 = (90.0 - proj%hemi*proj%truelat2)*rad_per_deg

   ! See if we are in the southern hemispere and flip the indices
   ! if we are. 
   if (proj%hemi == -1.0) then 
      inew = -x + 2.0
      jnew = -y + 2.0
   else
      inew = x
      jnew = y
   end if

   ! Compute radius**2 to i/j location
   xx = inew - proj%polei
   yy = proj%polej - jnew
   r2 = (xx*xx + yy*yy)
   r = sqrt(r2)/proj%rebydx
   
   ! Convert to lat/lon
   if (r2 == 0.0) then
      lat = proj%hemi * 90.0
      lon = proj%stdlon
   else
       
      ! Longitude
      lon = proj%stdlon + deg_per_rad * ATAN2(proj%hemi*xx,yy)/proj%cone
      lon = MOD(lon+360.0, 360.0)

      ! Latitude.  Latitude determined by solving an equation adapted 
      ! from:
      !  Maling, D.H., 1973: Coordinate Systems and Map Projections
      ! Equations #20 in Appendix I.  
        
      if (chi1 == chi2) then
         chi = 2.0*ATAN((r/TAN(chi1))**(1.0/proj%cone) * TAN(chi1*0.5))
      else
         chi = 2.0*ATAN((r*proj%cone/Sin(chi1))**(1.0/proj%cone) * TAN(chi1*0.5)) 
      end if
      lat = (90.0-chi*deg_per_rad)*proj%hemi
   end if

   if (lon > +180.0) lon = lon - 360.0
   if (lon < -180.0) lon = lon + 360.0

   if (trace_use_dull) call da_trace_exit("da_xyll_lc")

end subroutine da_xyll_lc


subroutine da_xyll_merc(x, y, proj, lat, lon)

   !-----------------------------------------------------------------------
   ! Compute the lat/lon from i/j for mercator projection
   !-----------------------------------------------------------------------

   implicit none

   real,intent(in)               :: x
   real,intent(in)               :: y    
   type(proj_info),intent(in)    :: proj
   real, intent(out)             :: lat
   real, intent(out)             :: lon 

   if (trace_use_frequent) call da_trace_entry("da_xyll_merc")

   lat = 2.0*ATAN(EXP(proj%dlon*(proj%rsw + y-1.0)))*deg_per_rad - 90.0
   lon = (x-1.0)*proj%dlon*deg_per_rad + proj%lon1
   if (lon > 180.0) lon = lon - 360.0
   if (lon < -180.0) lon = lon + 360.0

   if (trace_use_frequent) call da_trace_exit("da_xyll_merc")

end subroutine da_xyll_merc


subroutine da_xyll_ps(x, y, proj, lat, lon)

   ! This is the inverse subroutine da_of llij_ps.  It returns the 
   ! latitude and longitude of an x/y point given the projection info 
   ! structure.  

   implicit none

   real, intent(in)                    :: x    ! Column
   real, intent(in)                    :: y    ! Row
   type (proj_info), intent(in)        :: proj
    
   real, intent(out)                   :: lat     ! -90 -> 90 North
   real, intent(out)                   :: lon     ! -180 -> 180 East

   real                                :: reflon
   real                                :: scale_top
   real                                :: xx,yy
   real                                :: gi2, r2
   real                                :: arccos

   if (trace_use_frequent) call da_trace_entry("da_xyll_ps")

   ! Compute the reference longitude by rotating 90 degrees to the east
   ! to find the longitude line parallel to the positive x-axis.
   reflon = proj%stdlon + 90.0
   
   ! Compute numerator term of map scale factor
   scale_top = 1.0 + proj%hemi * Sin(proj%truelat1 * rad_per_deg)

   ! Compute radius to point of interest
   xx = x - proj%polei
   yy = (y - proj%polej) * proj%hemi
   r2 = xx**2 + yy**2

   ! Now the magic code
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
  
   ! Convert to a -180 -> 180 East convention
   if (lon > 180.0) lon = lon - 360.0
   if (lon < -180.0) lon = lon + 360.0

   if (trace_use_frequent) call da_trace_exit("da_xyll_ps")
     
end subroutine da_xyll_ps


subroutine da_set_lc(proj)

   !-----------------------------------------------------------------------
   ! Purpose: Initialize the remaining items in the proj structure for a
   ! lambert conformal grid.
   !-----------------------------------------------------------------------

   implicit none
    
   type(proj_info), intent(inout) :: proj

   real :: arg
   real :: deltalon1
   real :: tl1r
   real :: ctl1r

   if (trace_use_dull) call da_trace_entry("da_set_lc")

   ! Compute cone factor
   call da_lc_cone(proj%truelat1, proj%truelat2, proj%cone)
   if (print_detail_map) then
      write(unit=stdout, fmt='(A,F8.6)') 'Computed cone factor: ', proj%cone
   end if
   ! Compute longitude differences and ensure we stay out of the
   ! forbidden "cut zone"
   deltalon1 = proj%lon1 - proj%stdlon
   if (deltalon1 .gt. +180.0) deltalon1 = deltalon1 - 360.0
   if (deltalon1 .lt. -180.0) deltalon1 = deltalon1 + 360.0

   ! Convert truelat1 to radian and compute COS for later use
   tl1r = proj%truelat1 * rad_per_deg
   ctl1r = COS(tl1r)

   ! Compute the radius to our known lower-left (SW) corner
   proj%rsw = proj%rebydx * ctl1r/proj%cone * &
           (TAN((90.0*proj%hemi-proj%lat1)*rad_per_deg/2.0) / &
            TAN((90.0*proj%hemi-proj%truelat1)*rad_per_deg/2.0))**proj%cone

   ! Find pole point
   arg = proj%cone*(deltalon1*rad_per_deg)
   proj%polei = 1.0 - proj%hemi * proj%rsw * Sin(arg)
   proj%polej = 1.0 + proj%rsw * COS(arg)  
   if (print_detail_map) then
      write(unit=stdout,fmt='(A,2F10.3)') 'Computed pole i/j = ', proj%polei, proj%polej
   end if

   if (trace_use_dull) call da_trace_exit("da_set_lc")

end subroutine da_set_lc                             


subroutine da_set_ps(proj)

   ! Initializes a polar-stereographic map projection from the partially
   ! filled proj structure. This routine computes the radius to the
   ! southwest corner and computes the i/j location of the pole for use
   ! in llij_ps and ijll_ps.

   implicit none
 
   type(proj_info), intent(inout)    :: proj

   real :: ala1
   real :: alo1
   real :: reflon
   real :: scale_top

   if (trace_use) call da_trace_entry("da_set_ps")

   ! To define the cone factor for polar stereographic projection 
   proj%cone = 1.0

   reflon = proj%stdlon + 90.0

   ! Compute numerator term of map scale factor
   scale_top = 1.0 + proj%hemi * Sin(proj%truelat1 * rad_per_deg)

   ! Compute radius to lower-left (SW) corner
   ala1 = proj%lat1 * rad_per_deg
   proj%rsw = proj%rebydx*COS(ala1)*scale_top/(1.0+proj%hemi*Sin(ala1))

   ! Find the pole point
   alo1 = (proj%lon1 - reflon) * rad_per_deg
   proj%polei = proj%knowni - proj%rsw * COS(alo1)
   proj%polej = proj%knownj - proj%hemi * proj%rsw * Sin(alo1)
   if (print_detail_map) then
      write(unit=stdout,fmt='(A,2F10.1)') 'Computed (I,J) of pole point: ',proj%polei,proj%polej
   end if

   if (trace_use) call da_trace_exit("da_set_ps")

end subroutine da_set_ps


subroutine da_map_init(proj)

   !-----------------------------------------------------------------------
   ! Purpose: Initializes the map projection structure to missing values
   !-----------------------------------------------------------------------

   implicit none

   type(proj_info), intent(inout)  :: proj

   if (trace_use_dull) call da_trace_entry("da_map_init")

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

   if (trace_use_dull) call da_trace_exit("da_map_init")

end subroutine da_map_init


subroutine da_map_set(proj_code,lat1,lon1,knowni,knownj,dx,stdlon,truelat1,truelat2,latinc,loninc,proj)
   ! Given a partially filled proj_info structure, this routine computes
   ! polei, polej, rsw, and cone (if LC projection) to complete the 
   ! structure.  This allows us to eliminate redundant calculations when
   ! calling the coordinate conversion routines multiple times for the
   ! same map.
   ! This will generally be the first routine called when a user wants
   ! to be able to use the coordinate conversion routines, and it
   ! will call the appropriate subroutines based on the 
   ! proj%code which indicates which projection type  this is.

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

   if (trace_use_dull) call da_trace_entry("da_map_set")

   ! First, check for validity of mandatory variables in proj
   if (ABS(lat1) > 90.0) then
      call da_error("da_map_set.inc",29, &
         (/"Latitude of origin corner required as follows: -90N <= lat1 < = 90N"/))
   end if
   if (ABS(lon1) > 180.0) then
      call da_error("da_map_set.inc",33, &
         (/"Longitude of origin required as follows: -180E <= lon1 <= 180W"/))
   end if
   if ((dx .LE. 0.0).AND.(proj_code .NE. PROJ_LATLON)) then
      call da_error("da_map_set.inc",37, &
         (/"Require grid spacing (dx) in meters be positive!"/))
   end if
   if ((ABS(stdlon) > 180.0).AND.(proj_code .NE. PROJ_MERC)) then
      call da_error("da_map_set.inc",41, &
         (/"Need orientation longitude (stdlon) as: -180E <= lon1 <= 180W"/)) 
   end if
   if (proj%code .NE. PROJ_LATLON .and. ABS(truelat1)>90.0) then
      call da_error("da_map_set.inc",45, &
         (/"Set true latitude 1 for all projections!"/))
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
            write(unit=stdout,fmt='(A)') 'Setting up CYLINDRICAL EQUIDISTANT LATLON map...'
         end if
         ! Convert lon1 to 0->360 notation
         if (proj%lon1 < 0.0) proj%lon1 = proj%lon1 + 360.0
   
         proj%cone = 1.0                                  
      case default
         write(unit=message(1),fmt='(A,I2)') 'Unknown projection code: ', proj%code
         call da_error("da_map_set.inc",107,message(1:1))
   end select pick_proj
   proj%init = .true.

   if (trace_use_dull) call da_trace_exit("da_map_set")

end subroutine da_map_set


subroutine da_set_merc(proj)
  
   !--------------------------------------------------------------------------
   ! Purpose: Sets up the remaining basic elements for the mercator projection
   !--------------------------------------------------------------------------

   implicit none

   type(proj_info), intent(inout)       :: proj

   real :: clain

   if (trace_use) call da_trace_entry("da_set_merc")

   !  Preliminary variables

   clain = COS(rad_per_deg*proj%truelat1)
   proj%dlon = proj%dx / (earth_radius_m * clain)

   ! Compute distance from equator to origin, and store in the 
   ! proj%rsw tag.

   proj%rsw = 0.0
   if (proj%lat1 .NE. 0.0) then
      proj%rsw = (alog(tan(0.5*((proj%lat1+90.)*rad_per_deg))))/proj%dlon
   end if

   if (trace_use) call da_trace_exit("da_set_merc")

end subroutine da_set_merc


subroutine da_lc_cone(truelat1, truelat2, cone)

   !------------------------------------------------------------------------
   ! Purpose: compute the cone factor of a Lambert Conformal projection
   !------------------------------------------------------------------------

   implicit none
    
   real, intent(in)             :: truelat1  ! (-90 -> 90 degrees N)
   real, intent(in)             :: truelat2  !   "   "  "   "     "
   real, intent(out)            :: cone

   if (trace_use_dull) call da_trace_entry("da_lc_cone")

   ! First, see if this is a secant or tangent projection.  For tangent
   ! projections, truelat1 = truelat2 and the cone is tangent to the 
   ! Earth's surface at this latitude.  For secant projections, the cone
   ! intersects the Earth's surface at each of the distinctly different
   ! latitudes
   if (abs(truelat1-truelat2) > 0.1) then
      cone = alog10(cos(truelat1*rad_per_deg)) - &
             alog10(cos(truelat2*rad_per_deg))
      cone = cone /(alog10(tan((45.0 - abs(truelat1)/2.0) * rad_per_deg)) - &
             alog10(tan((45.0 - abs(truelat2)/2.0) * rad_per_deg)))        
   else
      cone = sin(abs(truelat1)*rad_per_deg)  
   end if

   if (trace_use_dull) call da_trace_exit("da_lc_cone")

end subroutine da_lc_cone


subroutine da_convert_zk (info)

   !-----------------------------------------------------------------------
   ! Purpose: Transfer obs. x to grid i and calculate its
   ! distance to grid i and i+1
   !-----------------------------------------------------------------------

   implicit none

   type(infa_type), intent(inout) :: info
   integer :: k, n

   if (trace_use) call da_trace_entry("da_convert_zk")

   do n = info%n1, info%n2

      do k = 1, info%levels(n)

         if ( (info%zk(k,n) > 0.0 .or. anal_type_verify) &
               .and. info%zk(k,n) .ne. missing_r) then 

            info%k(k,n) = int ( info%zk(k,n))
  
            if (info%k(k,n) < kts)  info%k(k,n) = kts
            if (info%k(k,n) >= kte) info%k(k,n) = kte-1

            info%dz(k,n) = info%zk(k,n) - real(info%k(k,n))
            info%dzm(k,n)= 1.0 - info%dz(k,n)

         else
            info%k(k,n) = 0
            info%dz(k,n) = 0.0
            info%dzm(k,n) = 0.0
         endif

      enddo

   enddo

   if (trace_use) call da_trace_exit("da_convert_zk")

end subroutine da_convert_zk



subroutine da_1d_eigendecomposition( bx, e, l )

   !------------------------------------------------------------------------------
   !  Purpose: Compute eigenvectors E and eigenvalues L of vertical covariance matrix
   !           B_{x} defined by equation:  E^{T} B_{x} E = L, given input 3D field of
   !           errors (sum over all horizontal locations).
   !------------------------------------------------------------------------------

   implicit none
   
   real, intent(in)         :: bx(:,:)          ! Global vert. background error.
   
   real, intent(out)        :: e(:,:)           ! Eigenvectors of Bx.
   real, intent(out)        :: l(:)             ! Global eigenvalues of Bx.
   
   integer                  :: kz               ! Size of 3rd dimension.
   integer                  :: m                ! Loop counters
   integer                  :: work             ! Size of work array.
   integer                  :: info             ! Info code.
   
   real*8, allocatable        :: ecopy(:,:)
   real*8, allocatable        :: lcopy(:)
   real*8, allocatable        :: work_array(:)

   if (trace_use_dull) call da_trace_entry("da_1d_eigendecomposition")

   !-------------------------------------------------------------------------
   ! [1.0]: Initialise:
   !-------------------------------------------------------------------------
   
   kz = size(bx, dim=1)
   
   !-------------------------------------------------------------------------
   ! [5.0]: Perform global eigenvalue decomposition using LAPACK software:
   !-------------------------------------------------------------------------
   
   work = 3 * kz - 1
   allocate( work_array(1:work) )
   
   allocate( ecopy(1:kz,1:kz) )
   allocate( lcopy(1:kz) )
   
   ecopy(:,:) = bx(:,:)

   lcopy = 0.0

   call dsyev( 'V', 'U', kz, ecopy, kz, lcopy, &
              work_array, work, info )
   
   if ( info /= 0 ) then
      write(unit=message(1),fmt='(A,I4,A)') &
         ' da_1d_eigendecomposition: info = ', &
         info,' - error in decomposition.'
      call da_error("da_1d_eigendecomposition.inc",54,message(1:1))
   end if
   
   !--Swap order of eigenvalues, vectors so 1st is one with most
   !  variance, etc:
   
   do m=1,kz
      l(m) = lcopy(kz+1-m)
      e(1:kz,m) = ecopy(1:kz,kz+1-m)
   end do
   
   deallocate (work_array)
   deallocate (ecopy)
   deallocate (lcopy)

   if (trace_use_dull) call da_trace_exit("da_1d_eigendecomposition")
   
end subroutine da_1d_eigendecomposition


subroutine da_obs_sfc_correction(info, sfc_obs, n, xb)

   !--------------------------------------------------------------------
   ! Purpose: correct the surface measurements (wind, 
   ! temperature, and pressure) from the observed height to the WRF     
   ! model's lowest half-zeta level before going to the minimization.   
   !                                                                    
   !   Wind       : based on the similarity theory                      
   !   Temperature: Frank Ruggiero's (1996) method                      
   !   Pressure   : Hydrostatic correction                              
   !                                                                    
   ! The order of the vertical index is "kts=1(bottom) and kte(top)".   
   ! With cv_options=2 and sfc_assi_option=1, this procedure must be    
   ! gone through, otherwise unrealistic results may be obtained.   
   !--------------------------------------------------------------------

   implicit none

   type(infa_type),  intent(in)    :: info
   type(synop_type), intent(inout) :: sfc_obs
   integer,          intent(in)    :: n
   type(xb_type),    intent(in)    :: xb

   real    :: roughness, psfc, mslp, dx, dxm, dy, dym, ho, po, to, qo
   real    :: hm, pm, tm, qm, um, vm, correc, val

   integer :: i, j, k
   real    :: t_mdl(kts:kte)
   real    :: q_mdl(kts:kte)
   real    :: u_mdl(kts:kte)
   real    :: v_mdl(kts:kte)
   real    :: height(kts:kte)
   real    :: pressure(kts:kte)

   if (trace_use_dull) call da_trace_entry("da_obs_sfc_correction")

   ! 1. Check if it needs to do the surface correction at the first level
 
   ! 1.1 Surface reports located at far below the lowest model level
 
   ! 2. Model profile at OBS site for surface correction

   i   = info%i(1,n)
   j   = info%j(1,n)
   dx  = info%dx(1,n)
   dy  = info%dy(1,n)
   dxm = info%dxm(1,n)
   dym = info%dym(1,n)

   ! Model roughness at the obs site

   roughness = dym*(dxm*xb%rough(i,j)   + dx*xb%rough(i+1,j)) &
      + dy *(dxm*xb%rough(i,j+1) + dx*xb%rough(i+1,j+1))

   do k = kts, kte
      pressure(k) = dym*(dxm*xb%p(i,j,k) + dx*xb%p(i+1,j,k)) + dy*(dxm*xb%p(i,j+1,k) + dx*xb%p(i+1,j+1,k))
      height(k)   = dym*(dxm*xb%h(i,j,k) + dx*xb%h(i+1,j,k)) + dy*(dxm*xb%h(i,j+1,k) + dx*xb%h(i+1,j+1,k))
      t_mdl(k)    = dym*(dxm*xb%t(i,j,k) + dx*xb%t(i+1,j,k)) + dy*(dxm*xb%t(i,j+1,k) + dx*xb%t(i+1,j+1,k))
      q_mdl(k)    = dym*(dxm*xb%q(i,j,k) + dx*xb%q(i+1,j,k)) + dy*(dxm*xb%q(i,j+1,k) + dx*xb%q(i+1,j+1,k))
      u_mdl(k)    = dym*(dxm*xb%u(i,j,k) + dx*xb%u(i+1,j,k)) + dy*(dxm*xb%u(i,j+1,k) + dx*xb%u(i+1,j+1,k))
      v_mdl(k)    = dym*(dxm*xb%v(i,j,k) + dx*xb%v(i+1,j,k)) + dy*(dxm*xb%v(i,j+1,k) + dx*xb%v(i+1,j+1,k))
   end do 

   ! 3. OBS data and recover the surface pressure from the
   ! mean sea level pressure (mslp)

   ho   = sfc_obs % h
   po   = sfc_obs % p % inv 
   to   = sfc_obs % t % inv
   qo   = sfc_obs % q % inv

   ! 3.1 Compute the surface OBS pressure from mean sea level pressure

   if ( psfc_from_slp ) then
      mslp = info%slp(n)%inv
      if (abs(mslp - missing_r) > 1.0) then
         psfc = missing_r
         if (abs(ho - missing_r) > 1.0) then
            if (abs(to - missing_r) > 1.0) then
               call da_sfcprs (kts, kte, pressure, t_mdl, q_mdl, height, psfc, mslp, ho, to)
            else
               call da_sfcprs (kts, kte, pressure, t_mdl, q_mdl, height, psfc, mslp, ho)
            end if
         end if
         sfc_obs % p % inv = psfc
         ! YRG: to allow assmilate the Psfc from mslp:
         sfc_obs % p % qc  = 0
      end if
      po = sfc_obs % p % inv
   end if

   if (sfc_obs % p % inv < 1.0) then
      sfc_obs % p % qc  = missing_data
   end if

   po = sfc_obs % p % inv

   ! 3.2 Check that obs pressure and height are present
   !     ----------------------------------------------

   if (abs(po - missing_r) < 1.0  .OR. abs(ho - missing_r) < 1.0) then
      if (trace_use_dull) call da_trace_exit("da_obs_sfc_correction")
      return

      ! write(unit=message(1), fmt='(/3(1x,a))') &
      !    'MISSinG HEIGHT OR PRESSURE OBSERVATION ID ', &
      !    trim (sfc_obs%info % id), trim (sfc_obs%info % name)

      ! write(unit=message(2), fmt='(2(A,F12.3,/))') &
      !                         ' height   = ',ho,&
      !                         ' pressure = ',po
      ! call da_error("da_obs_sfc_correction.inc",112,message(1:2))

   end if

   ! 4.  Bring surface observation below model levels with good height quality
   ! ================================================

   if (sfc_obs % h < height(kts)) then

      ! 2.3 Make use of local variables for model  
      !     -------------------------------------

      um = u_mdl(kts)
      vm = v_mdl(kts)
      tm = t_mdl(kts)
      pm = pressure (kts)
      qm = q_mdl(kts)
      hm = height(kts)

      ! 3.2 Correction wind based on similarity laws
      !     -------------------------------------------

      if ((abs(sfc_obs%u%inv - missing_r) > 1.0) .AND. (abs(sfc_obs%v%inv - missing_r) > 1.0)) then
         ! 3.2.1 Correction factor

         ! temperature and moisture at sigma level:

         if (abs(to - missing_r) < 1.0) then
            correc = da_mo_correction(ho, po, tm, qo, hm, pm, tm, qm, um ,vm, roughness)
         else
            correc = da_mo_correction(ho, po, to, qo, hm, pm, tm, qm, um ,vm, roughness)
         end if

         ! 3.2.2 Wind correction 
         !       ---------------

         !  Correct data and replace any previous wind qcs
         !  with surface correction flag

         sfc_obs % u % inv = correc * sfc_obs % u % inv 
         sfc_obs % u % qc  = surface_correction

         sfc_obs % v % inv = correc * sfc_obs % v % inv
         sfc_obs % v % qc  = surface_correction
      end if

      ! 3.4 Correct pressure
      !     ----------------

      if (sfc_obs % p % qc >= 0) then
         !  Correct data
         if (abs(to  - missing_r) > 1.0 .and. abs(qo - missing_r) > 1.0) then
            call da_intpsfc_prs (val, ho, po, hm, tm, qm, to, qo)
         else if (abs(to  - missing_r) > 1.0) then
            call da_intpsfc_prs (val, ho, po, hm, tm, qm, to)
         else
            call da_intpsfc_prs (val, ho, po, hm, tm, qm)
         end if

         !  Replace any previous pressure qc by the surface correction

         sfc_obs % p % inv = val
         sfc_obs % p % qc  = surface_correction
      end if

      ! 3.5 Correct temperature
      !     -------------------

      if (abs(sfc_obs % t % inv - missing_r) > 1.0) then
         !  Correct data
         call da_intpsfc_tem (val, ho, po, to, height, pressure, t_mdl, kts, kte)

         sfc_obs % t % inv = val

         !  Replace any previous temperature qc by the surface correction
         sfc_obs % t % qc  = surface_correction
      end if

      ! 3.6 Assign model lowest level height + 1m to observation
      !      ----------------------------------------------------- 
      ! sfc_obs % h = height(kts) + 1.0

      ! 3.7 Height QC
      !     ---------
      ! sfc_obs % height_qc = surface_correction
   end if

   if (trace_use_dull) call da_trace_exit("da_obs_sfc_correction")

end  subroutine da_obs_sfc_correction


subroutine da_sfcprs (kts, kte, p, t, q, height, psfc, pslv, h_obs, t_obs)
   
   ! Purpose: to compute the pressure at obs height (psfc) from the sea
   !          level pressure (pslv), obs height (h_obs), temperature
   !          (t_obs, optional), and the model profile (p,t,q).
    
   implicit none

   integer,        intent(in)  :: kts, kte
   real,           intent(in)  :: height(kts:kte)
   real,           intent(in)  :: p(kts:kte)
   real,           intent(in)  :: t(kts:kte)
   real,           intent(in)  :: q(kts:kte)
   real,           intent(in)  :: h_obs, pslv    
   real,           intent(out) :: psfc    
   real, optional, intent(in)  :: t_obs    

   real, parameter :: g = gravity
   real, parameter :: r =  gas_constant
   real, parameter :: rov2 = r / 2.0
   real, parameter :: gamma  = 6.5e-3   ! temperature lapse rate
   real, parameter :: gammarg= gamma*gas_constant/g

   real    :: plmb(3)                   

   integer :: k
   integer :: k500
   integer :: k700
   integer :: k850

   logical :: l1
   logical :: l2
   logical :: l3

   real    :: gamma78, gamma57  
   real    :: ht    
   real    :: p1 
   real    :: t1       
   real    :: t500       
   real    :: t700    
   real    :: t850    
   real    :: tc
   real    :: tfixed 
   real    :: tslv, tsfc  

   if (trace_use_dull) call da_trace_entry("da_sfcprs")

   TC   = t_kelvin + 17.5
   PLMB = (/85000.0, 70000.0, 50000.0/)

   !  Find the locations of the 850, 700 and 500 mb levels.

101 continue
   K850 = 0                              ! FinD K AT: P=850
   K700 = 0                              !            P=700
   K500 = 0                              !            P=500

   do K = kte-1, kts, -1
      if      (p(k) > PLMB(1) .and. K850==0) then
         K850 = K + 1
      else if (p(k) > PLMB(2) .and. K700==0) then
         K700 = K + 1
      else if (p(k) > PLMB(3) .and. K500==0) then
         K500 = K + 1
      end if
  
   end do

   if ((K850 == 0) .OR. (K700 == 0) .OR. (K500 == 0)) then
      ! write(unit=message(1),fmt='(A,3f8.0,A)') &
      !    'Error in finding p level for ',PLMB(1), PLMB(2), PLMB(3),' Pa.'
      ! do K = 1, KX
      !    write(unit=message(k+1),fmt='(A,I3,A,F10.2)') 'K = ',K,'  PRESSURE = ',P(K)
      ! end do
      ! write(unit=message(kx+2),fmt='(A,2f8.0,A,f8.0,A)) ','Expected ',    &
      !    PLMB(1), PLMB(2),' and ',PLMB(3),' Pa. values, at least.'
      ! message(kx+3)="not enough levels"
      ! message(kx+4)="Change the pressure levels by -25mb"
      ! call da_error("da_sfcprs.inc",79,message(1:kx+4))
      PLMB(1) = PLMB(1) - 2500.0
      PLMB(2) = PLMB(2) - 2500.0
      PLMB(3) = PLMB(3) - 2500.0
      goto 101
   end if

   !  The 850 hPa level is called something special, and interpolated
   !  to cross points.  And then, we fill those last rows and columns.

   ht = height(k850)

   !  The variable ht is now -h_obs/ht(850 hPa).  The plot thickens.

   ht = -h_obs / ht

   !  Make an isothermal assumption to get a first guess at the surface
   !  pressure.  This is to tell us which levels to use for the lapse
   !  rates in a bit.

   psfc = pslv * (pslv / p(k850)) ** ht

   !  Get a pressure more than 100 hPa above the surface - P1.  The
   !  P1 is the top of the level that we will use for our lapse rate
   !  computations.

   if      ((PSFC - p(k850) - 10000.0) >= 0.0) then
      P1 = p(k850)
   else if ((PSFC - p(k700)) >= 0.0) then
      P1 = PSFC - 10000.0
   else
      P1 = p(k500)
   end if

   !  Compute virtual temperatures for k850, k700, and k500 layers.  Now
   !  you see WHY? we wanted Q on pressure levels, it all is beginning   
   !  to make sense.

   t850 = t(k850) * (1.0 + 0.608 * q(k850))
   t700 = t(k700) * (1.0 + 0.608 * q(k700))
   t500 = t(k500) * (1.0 + 0.608 * q(k500))

   !  Compute two lapse rates between these three levels.  These are
   !  environmental values for each (i,j).

   gamma78 = LOG(t850 / t700)  / LOG (p(k850) / p(k700))
   gamma57 = LOG(t700 / t500)  / LOG (p(k700) / p(k500))

   if ((psfc - p(k850) - 10000.0) >= 0.0) then
      t1 = t850
   else if ((psfc - p(k850)) >= 0.0) then
      t1 = t700 * (p1 / p(k700)) ** gamma78
   else if ((psfc - p(k700)) >= 0.0) then 
      t1 = t500 * (p1 / p(k500)) ** gamma57
   else
      t1 = t500
   end if

   !  From our temperature way up in the air, we extrapolate down to
   !  the sea level to get a guess at the sea level temperature.

   tslv = t1 * (pslv / p1) ** (gammarg)

   !  The new surface temperature is computed from the with new sea level 
   !  temperature, just using the elevation and a lapse rate.  This lapse 
   !  rate is -6.5 K/km.

   if (present (t_obs)) then
     tsfc = t_obs
   else
     tsfc = tslv - gamma * h_obs
   end if

   !  A correction to the sea-level temperature, in case it is too warm.

   TFIXED = TC - 0.005 * (TSFC - TC) ** 2

   L1 = tslv .LT. tc
   L2 = tsfc .LE. tc
   L3 = .NOT. L1
   if (L2 .AND. L3) then
      tslv = tc
   else if ((.NOT. L2) .AND. L3) then
      tslv = tfixed
   end if

   !  Finally, we can get to the surface pressure.

   p1 = -h_obs * g / (rov2 * (tsfc + tslv))
   psfc = pslv * EXP(p1)

   !  Surface pressure and sea-level pressure are the same at sea level.

   if (ABS (h_obs) < 0.1) psfc = pslv

   if (trace_use_dull) call da_trace_exit("da_sfcprs")

end subroutine da_sfcprs


subroutine da_intpsfc_prs (val, ho, po, hm, tm, qm, to, qo)

   !----------------------------------------------------------------------------
   ! Purpose: Correct pressure between two levels. 
   !
   ! Reference: make use of the hydrosatic equation:
   !
   !  P2 = P1 * exp [-G/R * (z2-z1) / (tv1 + tv2)/2)
   !
   ! Where:
   !  z1  = height at level 1
   !  z1  = height at level 2
   !  tv1 = temperature at level 1
   !  tv2 = temperature at level 2
   !  P1  = Pressure at level 1
   !  P2  = Pressure at level 2
   !----------------------------------------------------------------------------

   implicit none

   real, intent (out)          :: val
   real, intent (in)           :: ho, po
   real, intent (in)           :: hm, tm, qm    
   real, intent (in), optional :: to, qo

   real :: tvo, tvm, tv, dz, arg

   if (trace_use) call da_trace_entry("da_intpsfc_prs")

   ! 1.  model and observation virtual temperature
   ! ---------------------------------------------

   tvm = tm  * (1.0 + 0.608 * qm)
   if (present(to) .and. present(qo)) then
      tvo = to  * (1.0 + 0.608 * qo)
   else if (present(to) .and. .not.present(qo)) then
      tvo = to
   else
      tvo = tvm
   end if

   tv  = 0.5 * (tvm + tvo)

   ! 2. height difference bewteen model surface and observations
   ! ------------------------------------------------------------

   dz = hm - ho
 
   ! 3.  extrapolate pressure obs to model surface
   ! ---------------------------------------------

   arg = dz * gravity / gas_constant
   arg = arg    / tv 

   val = po * exp (-arg)

   if (trace_use) call da_trace_exit("da_intpsfc_prs")

end subroutine da_intpsfc_prs


subroutine da_intpsfc_tem (val, ho, po, to, height, pressure, temp, kts, kte)

   !-----------------------------------------------------------------------
   ! Purpose: Correct temperature between two levels.
   !
   ! Reference: 
   ! ---------
   ! The use of surface observations in four dimensional data assmiilation
   !  Using a mesoscale model.
   !  Ruggiero et al., 1996, Monthly Weather Review, Volume 124, 1018-1033
   !
   !----------------------------------------------------------------------------

   implicit none

   real,    intent (out) :: val
   integer, intent (in)  :: kts, kte
   real,    intent (in)  :: ho, po, to
   real,    intent (in)  :: height(kts:kte)
   real,    intent (in)  :: pressure(kts:kte)
   real,    intent (in)  :: temp(kts:kte)

   real    :: prs_mb(kts:kte)
   ! calculated but never used
   !      real    :: dth_12, dth_21, dth_sfc, dth_obs
   !     real    :: dhe_12, dhe_21, dhe_sfc1, dhe_obs1, dhe_sfc2, dhe_obs2
   real    :: dth_21, dth_sfc, dth_obs
   real    :: dhe_12, dhe_sfc1, dhe_obs1, dhe_sfc2, dhe_obs2
   real    :: th_100mb, th_200mb, th_obs, th_sfc
   real    :: th_obs_int, th_sfc_int
   real    :: pdif, rcp
   integer :: k_100mb, k_200mb, jk
   integer :: inc_100mb, inc_200mb

   if (trace_use_dull) call da_trace_entry("da_intpsfc_tem")

   rcp = gas_constant/cp

   ! 1.Find levels: model surface + 100hpa and model surface + 200hpa ar obs loc
   ! ===========================================================================

   ! 1.1 Convert model pressure profile from Pa to hPa  

   prs_mb = pressure / 100.0

   ! 1.2  Find levels surface + 100hPA

   inc_100mb = 100.0
   k_100mb   = kts

   do jk =  kts+1, kte
      pdif = prs_mb (kts) - prs_mb (jk)
      if (pdif .GE. inc_100mb) then
         k_100mb = jk
         exit
      end if
   end do

   ! 1.2  Find levels surface + 200hPA

   inc_200mb = 200.0
   k_200mb   = kts

   do jk =  kts+1, kte
      pdif = prs_mb (kts) - prs_mb (jk)
      if (pdif .GE. inc_200mb) then
         k_200mb = jk
         exit
      end if
   end do

   ! 1.3  Check consistency 

   if ((k_100mb .LE. kts) .OR. (k_200mb .LE. kts) .OR. &
        (k_200mb .LT. k_100mb)) then
      write (unit=message(1),fmt='(A)') ' Cannot find sfc + 200hPa and sfc + 100hPa' 
      write (unit=message(2),fmt='(A,I2,A,F10.3)') ' P (',k_200mb,') = ',prs_mb (k_200mb) 
      write (unit=message(3),fmt='(A,I2,A,F10.3)') ' P (',k_100mb,') = ',prs_mb (k_100mb) 
      write (unit=message(4),fmt='(A,F10.3)')      ' P_SFC  = ',         prs_mb (kts) 
      call da_warning("da_intpsfc_tem.inc",80,message(1:4))
      val = missing_r
      if (trace_use_dull) call da_trace_exit("da_intpsfc_tem")
      return
   end if

   ! 2.  potential temperature 
   ! =========================

   ! 2.1 Potential temperature at 100hPa above model surface

   th_100mb = temp (k_100mb) * (1000.0 / prs_mb (k_100mb))**rcp

   ! 2.2 Potential temperature at 200hPa above model surface

   th_200mb = temp (K_200mb) * (1000.0 / prs_mb (k_200mb))**rcp

   ! 2.3 Potential temperature at observation location

   th_obs   = to * (1000.0 / (po/100.0)) ** rcp


   ! 3.  lapse rate between surface+200hpa and surface+100hpa
   ! =========================================================

   ! 3.1 Potential temperature increment
    
   dth_21 = th_100mb - th_200mb
   ! never used
   ! dth_12 = th_200mb - th_100mb

   ! 3.1 Height increments
    
   ! never used
   ! dhe_21   = height (k_100mb)- height (k_200mb)
   dhe_sfc1 = height (k_100mb)- height (kts)
   dhe_obs1 = height (k_100mb)- ho

   dhe_12   = height (k_200mb)- height (k_100mb)
   dhe_sfc2 = height (k_200mb)- height (kts)
   dhe_obs2 = height (k_200mb)- ho

   ! 3.2 Extrapolated potential temperature at model surface and observation loc

   th_sfc_int = th_100mb + (dth_21/dhe_12) * dhe_sfc1 
   th_obs_int = th_100mb + (dth_21/dhe_12) * dhe_obs1 


   ! 4.  Bring temperature onto model surface
   ! ========================================

   ! 4.1 Difference at observation locations

   dth_obs = th_obs_int - th_obs

   ! 4.2 Difference at model surface

   dth_sfc = (dhe_sfc2/dhe_obs2) * dth_obs

   ! 4.3 Potentiel temperature brought to model surface

   th_sfc = th_sfc_int - dth_sfc

   ! 4.3 Corresponding Temperature

   val  = th_sfc * (prs_mb (kts) / 1000.0)**rcp

   if (trace_use_dull) call da_trace_exit("da_intpsfc_tem")

end subroutine da_intpsfc_tem


function da_mo_correction (ho, po, to, qo, &
                        hm, pm, tm, qm, um ,vm, &
                        roughness)   RESULT (correc)

   !----------------------------------------------------------------------------
   ! Purpose: Compute the correction factor to convert surface wind into 40m 
   ! wind using similarity laws.
   !
   ! Reference:
   ! ---------
   !  A description of the fifth generation Penn State/NCAR Mesoscale Model
   !  Grell et al. 1994, page 29-30 and 80-83.
   !----------------------------------------------------------------------------

   implicit none

   real, intent (in)                :: ho, po, to, qo
   real, intent (in)                :: hm, pm, tm, qm, um, vm
   real, intent (in)                :: roughness

   real                             :: correc, winfac

   real :: thm , tho, tvm, tvo, thvm, thvo, rcp
   real :: za, Vc2, Va2, V2 
   ! FIX? real :: hdif, rib, rll, hll, zint
   real :: hdif, rib, hll, zint

   ! Height difference (in m) above wich correction is applied

   real, parameter :: hmax = 10.0, hs_max = 40.0  

   ! Default Roughness length im m (for land, set to 0 if on sea)

   real, parameter :: zint0 = 0.2

   rcp = gas_constant/cp

   ! 0.0  initialize correction factor to 1
   ! =====================================

   correc = 1.0
   winfac = 1.0

   ! 1.  height difference and roughness length
   ! ==========================================

   ! 1.1 Height difference
   !     -----------------

   hdif = hm - ho

   ! 1.2 No correction if height difference is less than hmax. 
   !     -----------------------------------------------------

   if (hdif <= hmax) then
      return
   end if

   ! too many
   ! if (trace_use) call da_trace_entry("da_mo_correction")

   ! 1.3 Compute the roughness length based upon season and land use 
   !     -----------------------------------------------------------

   zint = roughness

   if (zint < 0.0001) zint = 0.0001

   ! 2.  potential temperature at model surface and observation location
   ! ===================================================================

   ! 2.1 potential temperature on model surface
   !     --------------------------------------

   thm  = tm * (1000.0 / (pm/100.0)) ** rcp

   ! 2.2 potential temperature at observation location
   !     ---------------------------------------------

   tho  = to * (1000.0 / (po/100.0)) ** rcp

   ! 3.  virtual temperature at model surface and observation location
   ! ===================================================================

   ! 3.1 Virtual temperature on model surface
   !     -------------------------------------

   tvm  = tm * (1.0 + 0.608 * qm)

   ! 3.2 Virtual temperature at observation location
   !     -------------------------------------------

      tvo  = to * (1.0 + 0.608 * qo)

   ! 4.  potential virtual temperature at model surface and observation location 
   ! ===========================================================================

   ! 4.1 potential virtual temperature on model surface
   !     ----------------------------------------------

   thvm = tvm * (1000.0 / (pm/100.0)) ** rcp

   ! 4.2 potential virtual temperature at observation location
   !     -----------------------------------------------------

   thvo = tvo * (1000.0 / (po/100.0)) ** rcp


   ! 5.  bulk richardson number and moni-obukov length
   ! =================================================

   ! 5.1 Pre-calculations
   !     ----------------

   za  = hm - ho

   ! Because the following formula is derived based on
   ! the surface layer height is hs_max=40m. Under
   ! free convection, we assume that atmospheric state
   ! above the surface layer is fully mixed, and the
   ! wind at the lowest sigma half level is same as the
   ! wind at top of the surface layer. 

   if (za > hs_max) za = hs_max
      
   Va2 =   um*um + vm*vm
   Vc2 = 4.0 * MAX ((tho - thm), 0.0)

   V2  = Va2 + Vc2

   ! 5.2 Bulk richardson number
   !     ----------------------

   rib = (gravity * za / thm) * (thvm - thvo) / V2

   ! 5.3 Monin-obukov length
   !     -------------------

      hll = rib * LOG (za/zint)

   ! FIX? is this right?
   ! rll = za / hll

   ! 5.4 Ratio PBL height/Monin Obukov length: za/rll
   !     ------------------------------------

   hll =  ABS (hll)


   ! 6.  CORRECTION FACTOR BASED UPON REGIME
   ! =======================================

   ! 6.1 Stable conditions (REGIME 1)
   !     ---------------------------
 
   ! correc = 1.0      !  rib > 0.2

   ! 6.2 Mechanically driven turbulence (REGIME 2)
   !     ------------------------------------------

   ! correc = 1.0      !  0.0 =< rib <= 0.2

   correc = 1.0

   if (rib < 0.0) then

      ! 6.3 Unstable Forced convection (REGIME 3)
      !     -------------------------------------

      ! correc = 1.0  !   hll <= 1.5


      ! 6.4 Free convection (REGIME 4)
      !     --------------------------

      if (hll > 1.5) then
         if (zint < 0.2) then
            correc = 1.000 + 0.320 * zint ** 0.200
         else if (zint < 0.0) then
            correc = 1.169 + 0.315 * zint
         end if

         ! 6.4.1 The factor depended on Za (MM5, param.F)
      
         winfac = 0.5*(log(za/0.05)/log(40.0/0.05) &
                       + log(za)/log(40.0))

         correc =  correc * winfac
      end if
   end if

   ! if (trace_use) call da_trace_exit("da_mo_correction")

end function da_mo_correction


real function da_diff_seconds (date_char_1, date_char_2)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   character(len=24), intent(in) :: date_char_1, date_char_2
  
   integer :: ccyy_1,mo_1,dd_1,hh_1,mi_1,ss_1,jd_1
   integer :: ccyy_2,mo_2,dd_2,hh_2,mi_2,ss_2,jd_2
   integer :: i, year, diff_days
   integer :: start_year, end_year
   integer :: mmday(12)
 
   real    :: s_1, s_2

   if (trace_use_dull) call da_trace_entry("da_diff_seconds")
  
   mmday=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  
   read(date_char_1(1:19), fmt='(i4,1x,4(i2,1x),i2)') &
        ccyy_1, &
          mo_1, &
          dd_1, &
          hh_1, &
          mi_1, &
          ss_1
  
   read(date_char_2(1:19), fmt='(i4,1x,4(i2,1x),i2)') &
        ccyy_2, &
          mo_2, &
          dd_2, &
          hh_2, &
          mi_2, &
          ss_2

   if (ccyy_2 >= ccyy_1) then
      start_year = ccyy_1
      end_year   = ccyy_2
   else
      start_year = ccyy_2
      end_year   = ccyy_1
   end if

   diff_days = 0
  
   do year=start_year,end_year-1
      diff_days = diff_days + 365
      if (mod(year,4) == 0) then
         diff_days = diff_days + 1

         if ((mod(year,100) == 0) .and. (mod(year,400) /= 0)) then
            diff_days = diff_days - 1
         end if
      end if
   end do

   if (mod(ccyy_1,4) == 0) then
      mmday(2) = 29

      if((mod(ccyy_1,100) == 0) .and. (mod(ccyy_1,400) /= 0)) then
         mmday(2) = 28
      end if
   end if

   jd_1 = dd_1

   do i=1,mo_1-1
      jd_1=jd_1+mmday(i)
   end do

   s_1 = real(ss_1) &
       + 60.0*(real(mi_1) &
       + 60.0*(real(hh_1) &
       + 24.0* real(jd_1)))

   if (mod(ccyy_2,4) == 0) then
      mmday(2) = 29

      if((mod(ccyy_2,100) == 0) .and. (mod(ccyy_2,400) /= 0)) then
         mmday(2) = 28
      end if
   end if

   if (ccyy_2 >= ccyy_1) then
      jd_2 = dd_2 + diff_days
   else
      jd_2 = dd_2 - diff_days
   end if

   do i=1,mo_2-1
      jd_2=jd_2+mmday(i)
   end do

   s_2 = real(ss_2) &
       + 60.0*(real(mi_2) &
       + 60.0*(real(hh_2) &
       + 24.0* real(jd_2)))

   da_diff_seconds = abs(s_1 - s_2)

   if (trace_use_dull) call da_trace_exit("da_diff_seconds")

end function da_diff_seconds


function da_residual(n, k, yy, ox, n_bad) result (rr)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,              intent(in)    :: n, k
   type (field_type),    intent(in)    :: ox
   real,                 intent(in)    :: yy
   type (bad_info_type), intent(inout) :: n_bad

   real :: rr

   if (trace_use_frequent) call da_trace_entry("da_residual")

   ! WHY?
   if (k ==0 .OR. n==0) then
   end if

   rr = 0.0

   if (ox % qc >= obs_qc_pointer) then
      ! WHY?
      ! n_bad % num % use = n_bad % num % use + 1
      rr = ox % inv - yy
      ! else
      !    if (ox % qc /= -88) then
      !       n_bad % num % bad = n_bad % num % bad + 1
      !       n_bad % nn(n_bad % num % bad) = n
      !       n_bad % kk(n_bad % num % bad) = k
      !     else
      !       n_bad % num % miss = n_bad % num % miss + 1
      !    end if
   end if

   if (trace_use_frequent) call da_trace_exit("da_residual")

end function da_residual


subroutine da_residual_new(yy, qc, inv, rr)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)  :: qc(:,:)
   real,    intent(in)  :: inv(:,:)
   real,    intent(in)  :: yy(:,:)
   real,    intent(out) :: rr(:,:)

   if (trace_use) call da_trace_entry("da_residual_new")
  
   where (qc(:,:) >= obs_qc_pointer)
      rr(:,:) = inv(:,:) - yy(:,:)
   elsewhere
      rr(:,:) = 0.0
   endwhere

   if (trace_use) call da_trace_exit("da_residual_new")

end subroutine da_residual_new


subroutine da_add_noise (fld, yo, z)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (field_type), intent(inout) :: fld     ! O-B sub-type.
   real, intent(inout)              :: yo      ! Observation.
   real, intent(out)                :: z       ! Random number.

   real                             :: noise

   if (trace_use) call da_trace_entry("da_add_noise")

   z = missing_r

   if (fld % qc >= obs_qc_pointer) then
      ! [1] Calculate scaled Gaussian noise:

      call da_gauss_noise (z)      
      noise = fld % error * z
      
      ! [3] Recalculate corresponding O and O-B:
      yo = yo + noise
      fld % inv = fld % inv + noise
   end if

   if (trace_use) call da_trace_exit("da_add_noise")

end subroutine da_add_noise


subroutine da_add_noise_new (qc, error, inv, yo, z)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)              :: qc
   real, intent(in)                 :: error
   real, intent(inout)              :: inv
   real, intent(inout)              :: yo      ! Observation.
   real, intent(out)                :: z       ! Random number.

   real                             :: noise

   if (trace_use) call da_trace_entry("da_add_noise_new")

   z = missing_r

   if (qc >= obs_qc_pointer) then
      ! [1] Calculate scaled Gaussian noise:

      call da_gauss_noise (z)      
      noise = error * z
      
      ! [3] Recalculate corresponding O and O-B:
      yo = yo + noise
      inv = inv + noise
   end if

   if (trace_use) call da_trace_exit("da_add_noise_new")

end subroutine da_add_noise_new


subroutine da_max_error_qc (it, info, n, field, max_error,failed)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in)     :: it
   type(infa_type), intent(in)      :: info
   integer,           intent(in)    :: n
   type(field_type),  intent(inout) :: field
   real,              intent(in)    :: max_error
   logical,           intent(out)   :: failed

   real                               :: err, err_max
   integer                            :: qc_flag

   if (trace_use_frequent) call da_trace_entry("da_max_error_qc")

   failed = .false.

   qc_flag = field % qc
   err_max = field % error * max_error
   err     = field % inv
   err     = ABS (err)

   if (err > err_max) then
      field % qc = fails_error_max 
      failed = .true.
      field % inv = 0.0
   end if

   if (trace_use_frequent) call da_trace_exit("da_max_error_qc")

end subroutine da_max_error_qc


subroutine da_random_omb(std_dev, yo, qc_flag, omb)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none
   
   real, intent(in)     :: std_dev      ! Standard deviation scaling.
   real, intent(inout)  :: yo           ! Observatopn.
   integer, intent(in)  :: qc_flag      ! Ob QC flag.
   real, intent(inout)  :: omb          ! O-B.

   real :: xb           ! Background value.
   real :: z            ! Gaussian noise.

   if (trace_use) call da_trace_entry("da_random_omb")

   if (qc_flag >= obs_qc_pointer) then
      ! [1] Calculate background value from old O, O-B:
      xb = yo - omb
   
      ! [2] Calculate new O-B as scaled Gaussian noise:

      call da_gauss_noise(z)   
      omb = set_omb_rand_fac * std_dev * z

      ! [3] Recalculate corresponding observation:
      yo = xb + omb
   end if

   if (trace_use) call da_trace_exit("da_random_omb")

end subroutine da_random_omb


subroutine da_set_randomcv(cv_size, rcv)

   !---------------------------------------------------------------------------
   ! Purpose: Allocate components of control variable.
   !
   ! Method:  Allocate component in turn.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: cv_size
   real, intent(inout) :: rcv(1:cv_size)     ! Control variable structure.

   integer             :: i
   real                :: mean_cv, rms_cv, std_dev_cv
   real                :: sum_cv, sum_cv2

   if (trace_use) call da_trace_entry("da_set_randomcv")

   ! [1] Initialize random number generator and scalars:

   call da_random_seed

   sum_cv = 0.0
   sum_cv2 = 0.0

   ! [2] Calculate random numbers with Gaussian distribution:

   do i = 1, cv_size
      call da_gauss_noise(rcv(i))

      sum_cv = sum_cv + rcv(i)
      sum_cv2 = sum_cv2 + rcv(i) * rcv(i)
   end do

   mean_cv = sum_cv / real(cv_size)
   rms_cv = sqrt(sum_cv2 / real(cv_size))
   std_dev_cv = sqrt(rms_cv * rms_cv - mean_cv * mean_cv)

   write(unit=stdout,fmt=*)
   write(unit=stdout,fmt='(a)')' Gaussian (Normal) noise statistics:'
   write(unit=stdout,fmt='(a,f15.5)')' Mean = ',mean_cv
   write(unit=stdout,fmt='(a,f15.5)')' RMS = ', rms_cv
   write(unit=stdout,fmt='(a,f15.5)')' STD DEV = ', std_dev_cv

   if (trace_use) call da_trace_exit("da_set_randomcv")

end subroutine da_set_randomcv


real function da_gaus_noise(KDUM)

   !---------------------------------------------------------------------
   ! Purpose: return a normally(gaussian) distributed random variable ctions
   !     with 0 mean and 1 standard deviation
   !
   !   method :
   !   ------
   !
   !   input :
   !   -----
   !      argument
   !         kdum             : seed for random number generator 
   !                           (set to a large negative number on entry)
   !      common              : none
   !
   !   output :
   !   ------
   !      argument
   !         gauss_noise      : gaussian random betwen 
   !      common              : none
   !
   !   workspace :              none
   !   ---------
   !      local :
   !
   !   external :  
   !   --------
   !      unifva     
   !
   !   reference :
   !   ---------
   !      numerical recipes in fortran. the art of scientific computing.
   !      second edition.  cambridge university press.
   !      press et. al., 1986.
   !
   !   modifications:
   !   --------------
   !       original  : 95-01(f. vandenberghe)
   !       addition  : 96-06(a. weaver)
   !---------------------------------------------------------------------
 
   implicit none
 
   integer, intent(inout) :: kdum

   integer             :: niset
   real                :: gset

   save niset,gset
   data niset/0/

   real zfac,zrsq,zv1,zv2

   if (trace_use) call da_trace_entry("da_gaus_noise")

   if (niset.eq.0) then

1000  continue
      zv1   = 2.0*da_unifva(kdum) - 1.0
      zv2   = 2.0*da_unifva(kdum) - 1.0
      zrsq  = zv1**2 + zv2**2
          
      if ((zrsq>=1.0).or.(zrsq==0.0)) goto 1000

      zfac  = sqrt(-2.0*log(zrsq)/zrsq)
      gset  = zv1*zfac
      da_gaus_noise = zv2*zfac
      niset = 1
   else
      da_gaus_noise = gset
      niset = 0
   end if

   if (trace_use) call da_trace_exit("da_gaus_noise")

end function da_gaus_noise


subroutine da_openfile (kunit, cdfile, cdaccess, cdstatus, cdform)

   !------------------------------------------------------------------------
   ! Purpose: open an file and rewind
   !
   !   method:
   !   ------
   !
   !   input:
   !   -----
   !      kunit:          logical unit
   !      cdfile:         name of file for output of discared obs.
   !      cda!ess:       a!ess (sequential or),
   !      cdstatus:       status (old, new or unknown)
   !      cdform:         form (formatted or unformatted)
   !
   !   output:
   !   ------
   !      opened file
   !
   !   common:                      no
   !   -------
   !   external:                    no
   !   --------
   !   references:                  no
   !   ----------
   !
   !   modifications:
   !   --------------
   !       original :  98-07 (f. vandenberghe)
   !       additions : 98-11 norm doctor (f. vandenberghe)
   !--------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: Kunit
   character*(*), intent(in) :: CDFILE, CDACCESS, CDSTATUS, CDFORM

   integer :: iost

   ! 1.  open FILE
   ! -------------

   IOST = 0

   open (unit   = Kunit, &
          FILE   = CDFILE, &
          ACCESS = CDACCESS, &
          STATUS = CDSTATUS, &
          FORM   = CDFORM, &
          ERR    =  2000, &
          iostat = IOST       )
 
   rewind (Kunit)
 
   return
 
  ! 2.  ERROR PROCESSinG
   ! --------------------
 
2000 continue
   call da_error("da_openfile.inc",62, &
      (/"Cannot open file"//CDFILE/))

end subroutine da_openfile


subroutine da_smooth_anl(slab,imx,jmx,kx,npass,icrsdot)

   !-----------------------------------------------------------------------
   ! Purpose: spatially smooth (usually slab) to remove high
   ! frequency waves
   !-----------------------------------------------------------------------

   implicit none
   
   real,    intent(inout) :: SLAB(:,:,:)
   integer, intent(in)    :: imx, jmx, kx
   integer, intent(in)    :: npass
   integer, intent(in)    :: icrsdot
   
   real, allocatable :: SLABNEW(:,:)
   real              :: XNU(1:2)
   integer           :: ie, je, k 
   integer           :: loop, n, i, j

   if (trace_use) call da_trace_entry("da_smooth_anl")
   
   allocate (slabnew(imx,jmx))

   ie=imx-1-icrsdot
   je=jmx-1-icrsdot
   xnu(1)=0.50
   xnu(2)=-0.52
   do k=1,kx
      do loop=1,npass*2
         n=2-mod(loop,2)
 
         ! first smooth in the imx direction
 
         do i=2,ie
            do j=2,je
               slabnew(i,j)=slab(i,j,k)+xnu(n) * &
               ((slab(i,j+1,k)+slab(i,j-1,k))*0.5-slab(i,j,k))
            end do
         end do
         do i=2,ie
            do j=2,je
               slab(i,j,k)=slabnew(i,j)
            end do
         end do
 
         ! now smooth in the jmx direction
 
         do j=2,je
            do i=2,ie
               slabnew(i,j)=slab(i,j,k)+xnu(n) * &
               ((slab(i+1,j,k)+slab(i-1,j,k))*0.5-slab(i,j,k))
            end do
         end do

         do i=2,ie
            do j=2,je
               slab(i,j,k)=slabnew(i,j)
            end do
         end do
      end do
   end do

   deallocate (slabnew)

   if (trace_use) call da_trace_exit("da_smooth_anl")

end subroutine da_smooth_anl


subroutine da_togrid_new (x, ib, ie, i, dx, dxm)

   !-----------------------------------------------------------------------
   ! Purpose: Transfer obs. x to grid i and calculate its
   ! distance to grid i and i+1
   !-----------------------------------------------------------------------

   implicit none

   real,    intent(in)  :: x(:,:)
   integer, intent(in)  :: ib, ie
   real,    intent(out) :: dx(:,:), dxm(:,:)
   integer, intent(out) :: i(:,:)

   if (trace_use) call da_trace_entry("da_togrid_new")

   where (x(:,:) > 0.0) 
      i(:,:) = int (x(:,:))

      where(i(:,:) < ib)  i(:,:) = ib
      where(i(:,:) >= ie) i(:,:) = ie-1

      dx(:,:) = x(:,:) - real(i(:,:))
      dxm(:,:)= 1.0 - dx(:,:)
   elsewhere 
      i(:,:) = 0
      dx(:,:) = 0.0
      dxm(:,:) = 0.0
   end where

   if (trace_use) call da_trace_exit("da_togrid_new")

end subroutine da_togrid_new


subroutine da_togrid (x, ib, ie, i, dx, dxm)

   !-----------------------------------------------------------------------
   ! Purpose: Transfer obs. x to grid i and calculate its
   ! distance to grid i and i+1
   !-----------------------------------------------------------------------

   implicit none

   real,                     intent(in)  :: x
   integer,                  intent(in)  :: ib, ie
   real,                     intent(out) :: dx, dxm
   integer,                  intent(out) :: i

   if (trace_use_dull) call da_trace_entry("da_togrid")
   
   i = int (x)

   if (i < ib) then
      i = ib
   else if (i >= ie) then
      i = ie - 1
   end if

   dx = x - real(i)
   dxm= 1.0 - dx

   if (trace_use_dull) call da_trace_exit("da_togrid")

end subroutine da_togrid


real function da_unifva (kdum) 

   !--------------------------------------------------------------------
   ! Purpose: Minimal random number generator of Park and Miller with 
   ! Bays-Durham shuffle and added safeguards.
   ! Returns a uniform random deviate between 0.0. and 1.0 (exclusive 
   ! of the endpoint values). Call with kdum a negative integer to 
   ! initialize; thereafter, do not alter kdum between successive 
   ! deviates in sequence. rnmx should approximate the largest 
   ! floating value less than 1. 
   !
   ! See descripiton of function 'ran1', pg. 271.
   !--------------------------------------------------------------------
 
   implicit none
 
   integer, intent(inout) ::   KDUM

   integer JPIA,JPIM,JPIQ,JPIR,JPNTAB,JPNDIV
   real PPAM,PPEPS,PPRNMX

   parameter(JPIA=16807,JPIM=2147483647,JPIQ=127773,JPIR=2836, &
             JPNTAB=32,JPNDIV=1+(JPIM-1)/JPNTAB, &
             PPAM=1./JPIM,PPEPS=1.2E-07,PPRNMX=1.-PPEPS)

   integer JJ
   integer IJJ,IK
 
   integer NIV(JPNTAB),NIY
   save NIV,NIY
   DATA NIV /JPNTAB*0/, NIY /0/

   if (trace_use_frequent) call da_trace_entry("da_unifva")

   ! begin main
   ! ----------

   if ((KDUM.LE.0).OR.(NIY.EQ.0)) then
      KDUM = MAX(-KDUM , 1)

       do JJ = JPNTAB+8,1,-1
          IK   = KDUM/JPIQ
          KDUM = JPIA*(KDUM - IK*JPIQ) - JPIR*IK
 
          if (KDUM.lt.0) KDUM = KDUM + JPIM
          if (JJ.LE.JPNTAB) NIV(JJ) = KDUM
 
       end do

       NIY = NIV(1)
   end if
  
   IK   = KDUM/JPIQ
   KDUM = JPIA*(KDUM - IK*JPIQ) - JPIR*IK
     
   if (KDUM.LT.0) KDUM = KDUM + JPIM

   IJJ      = 1 + NIY/JPNDIV
   NIY      = NIV(IJJ)
   NIV(IJJ) = KDUM
   DA_UNifVA   = Min(PPAM*NIY , PPRNMX)

   if (trace_use_frequent) call da_trace_exit("da_unifva")

end function da_unifva


SUBROUTINE da_buddy_qc ( numobs, m_max, station_id, xob, yob, obs, qc_flag_small,    &
                         name, pressure, dx, buddy_weight , buddy_difference, &
                         iunit, print_buddy, xob_g, yob_g, obs_g, qc_flag_small_g, num_recv)

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Imported from BUDDY in RAWINS

!   Notes:
!      This routine performs error checks by comparing the difference
!         between the first guess field
!         and observations at the observing locations (xob, yob)
!         with the averaged difference at the nearby stations within
!         a given distance. If the difference value at the observation
!         point differs from the average of the nearby difference values
!         by more than a certain maximum (dependent on variable types,
!         pressure levels, and BUDWGT parameter in namelist), the
!         observation is flagged as suspect, and add the flag message
!         to the observation records.

!   *** Note that x and y are in natural coordinates now.
!         Need to check if the x and y are used correctly...

!   other variables:
!    err:       difference between observations and first guess fields
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!  This subroutine da_buddy_qc is adapted by Y.-R.Guo, NCAR/MMM, based on
!  the WRFV3/OBSGRID/src/qc2_module.F90. 
!
!  The algorithm for buddy check used here could be named as "Two Passes 
!  Innovation Buddy Check": the buddy members to a specific obs are
!  determined by going through two times of checks for innovation deprture 
!  from the buddy member's mean values. 
!
!  In OBSGRID, the buddy check is completed at each of  the pressure levels 
!  for all obs, i.e. all obs types are mixed to gether. In order to use this
!  buddy check in wrfvar, first we should sort the obs into the different 
!  pressure bins for 3-D observations,such as SOUND, then apply the buddy 
!  check for each of the obs types in da_get_innov_vector_?????.inc. For the 
!  2-D observations, such as SYNOP, no binning is needed.
!  
!  Most of the tolerances for each of variables at the different pressure 
!  levels are adapted from WRFV3/OBSGRID. The tolerances for specific humidity
!  used here are converted from rh and specified temperature. The tolerance
!  for PMSL is reduced to 350 Pa defined in da_control/da_control.f90.
!
!                                  Yong-Run Guo, 10/10/2008
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

  IMPLICIT NONE

  INTEGER,                       intent(in)    :: numobs, m_max, iunit
  REAL,    DIMENSION ( m_max ), intent(in)    :: obs , xob , yob
  INTEGER, DIMENSION ( m_max ), intent(inout) :: qc_flag_small
  REAL,    DIMENSION ( m_max, 0:num_procs-1 ), intent(in), optional :: obs_g , xob_g , yob_g
  INTEGER, DIMENSION ( m_max, 0:num_procs-1 ), intent(inout), optional :: qc_flag_small_g
  INTEGER, DIMENSION ( 0:num_procs-1 ),        intent(in), optional :: num_recv
  REAL,                    intent(in)          :: buddy_weight,   &
                                                  buddy_difference , &
                                                  dx          , &
                                                  pressure
  CHARACTER(LEN = 5), DIMENSION(m_max), intent(in) :: station_id
  CHARACTER(LEN = 8),                    intent(in) :: name
  LOGICAL           ,                    intent(in) :: print_buddy

  !  err:         difference between observations and first guess fields

  INTEGER                                    :: num, num3, num4, numj, &
                                                kob, n, ierr
  REAL                                       :: range,              &
                                                distance,           &
                                                sum,                &
                                                average,            &
                                                x, y,               &
                                                diff_numj,          &
                                                diff_check_1,       &
                                                diff_check_2
  REAL                                          plevel,             &
                             pp, tt, rh, t_error, rh_error, q_error
  REAL  , DIMENSION ( numobs )               :: err,                &
                                                difmax,             &
                                                diff_at_obs
  REAL  , DIMENSION ( m_max*num_procs )      :: diff
  INTEGER, DIMENSION (3)                     :: i_dummy
  INTEGER , DIMENSION ( numobs )             :: buddy_num,          &
                                                buddy_num_final
  REAL                                       :: scale
  INTEGER, DIMENSION ( m_max )               :: qc_flag


  real :: f_qv_from_rh
  external f_qv_from_rh

! [1.0] Store original qc flags:
  qc_flag = qc_flag_small

  plevel = pressure

! [2.0] Define distance within which buddy check is to find neighboring stations:

 ! 2.1 Define the horizontal range for buddy check:
  IF ( plevel .LE. 201.0 ) THEN
     range = 650.
  ELSE IF ( plevel .LE. 1000.1 ) THEN
     range = 300.
  ELSE
     range = 500.
  END IF

 ! 2.2 Define tolerance value for comparison:

  error_allowance_by_field : SELECT CASE ( name )

     CASE ( 'UU      ' , 'VV      ' )
        IF      ( plevel .LT. 401. ) THEN
           difmax = buddy_difference * 1.7
        ELSE IF ( plevel .LT. 701. ) THEN
           difmax = buddy_difference * 1.4
        ELSE IF ( plevel .LT. 851. ) THEN
           difmax = buddy_difference * 1.1
        ELSE
           difmax = buddy_difference
        END IF

     CASE ( 'PMSL    ' )
        difmax = buddy_difference

     CASE ( 'TT      ' )
        IF      ( plevel .LT. 401. ) THEN
           difmax = buddy_difference * 0.6
        ELSE IF ( plevel .LT. 701. ) THEN
           difmax = buddy_difference * 0.8
        ELSE
           difmax = buddy_difference
        END IF

     CASE ( 'QQ      ' )
!         difmax = buddy_difference

    ! Convert the rh error to q error:

        t_error = 2.0
        rh_error = buddy_difference

        pp = pressure * 100.0
        if ( pp > 85000.0 ) then
          rh = 90.0
          tt = 300.0
        else if ( pp > 70000.0 ) then
          rh = 80.0
          tt = 290.0
        else if ( pp > 35000.0 ) then
          rh = 70.0
          tt = 280.0
        else if ( pp > 10000.0 ) then
          rh = 60.0
          tt = 260.0
        else
          rh = 50.0
          tt = 240.0
        endif

        q_error = f_qv_from_rh( rh_error, t_error, rh, tt, pp)
        difmax = q_error

!         print '("p, t, rh, t_error, rh_error, difmax:",5f12.2,e15.5)', &
!            pp, tt, rh, t_error, rh_error, q_error

  END SELECT error_allowance_by_field

  difmax = difmax * buddy_weight

! [3.0] Compute the errors for buddy check:  

 ! 3.1 Loop through station locations (numobs) to compute buddy average
 !      at each observation locations.

!   station_loop_2 : DO num = numobs, 1, -1
  station_loop_2 : DO num = 1, numobs

     average = 0.0
 ! No buddy check is applied to a bad data (qc < 00:
     if ( qc_flag_small(num) < 0 ) cycle  station_loop_2

     !  Compute distance between all pairs of observations, 
     !  find number of observations within a given distance: buddy_num,
     !  and compute the difference between first guess and obs.

     buddy_num(num) = 0
     if (present(xob_g) .and. present(yob_g) .and. present(obs_g) .and. present(qc_flag_small_g)) then
     else
         call da_error("da_buddy_qc.inc",198, &
                      (/"Buddy check in parallel mode needs global observations"/))
     endif
     pe_loop : DO n = 0, num_procs-1
     station_loop_3 : DO num3 = 1, num_recv(n)

        x = xob(num) - xob_g(num3,n)
        y = yob(num) - yob_g(num3,n)
        IF ( ABS(x) .GT. 1.E-20 .OR. ABS(y) .GT. 1.E-20 ) THEN
           distance = SQRT ( x*x + y*y ) * dx
           IF ( distance .LE. range .AND. distance .GE. 0.1 &
                .and. qc_flag_small_g(num3,n) >= 0 ) THEN
! Only good data are included in buddy check group:
              buddy_num(num) = buddy_num(num) + 1
! innovation (O-B) ==> diff:
              diff(buddy_num(num)) = obs_g(num3,n)
           END IF
        END IF

     END DO station_loop_3
     END DO pe_loop

! innovation at the specific obs(num):
     diff_at_obs(num) = obs(num)
!
! Summation of innovations over the surrounding obs:
     sum = 0.
     DO numj = 1, buddy_num(num)
        sum = sum + diff(numj)
     END DO

     !  Check to see if there are any buddies, compute the mean innovation.

     IF ( buddy_num(num) .GT. 0 ) average = sum / buddy_num(num)

 ! 3.2 Check if there is any bad obs among the obs located within the 
 !  the radius surrounding the test ob.

     diff_check_1 = difmax (1) * 1.25
     diff_check_2 = difmax (1)

     check_bad_ob : IF ( buddy_num(num) .GE. 2 ) THEN

        kob = buddy_num(num)
        remove_bad_ob : DO numj = 1, buddy_num(num)
           diff_numj = ABS ( diff ( numj ) - average )

           IF ( diff ( numj ) .GT. diff_check_1  &
               .AND. diff_numj .GT. diff_check_2 ) THEN
    ! Bad obs:                 Innovation itself: diff(numj) > diff_check_1
    !        The distance between the innovation and average > diff_check_2
              kob = kob - 1
              sum = sum - diff ( numj )
           END IF
        END DO remove_bad_ob

! The final number of buddies:
        buddy_num_final(num) = kob

        !  We may have removed too many observations.

        IF ( kob .GE. 2 ) THEN
        !  Information for buddy check for specific obs(num)
           average = sum / kob
           err(num) = diff_at_obs(num) - average
        ELSE
        !  No buddy check completed for this specific obs(num):
           err(num)     = 0.
! Not change the flag (YRG, 02/10/2009)
!            qc_flag_small(num) = qc_flag_small(num) + no_buddies
! Not change the flag (YRG, 02/10/2009)
!            qc_flag_small_g(num,myproc) = qc_flag_small_g(num,myproc) + no_buddies
           i_dummy(1) = qc_flag_small_g(num,myproc)
           i_dummy(2) = num
           i_dummy(3) = myproc
           call mpi_bcast ( i_dummy, 3, mpi_integer, myproc, comm, ierr )
           qc_flag_small_g(i_dummy(2),i_dummy(3)) = i_dummy(1)
        END IF

     ELSE check_bad_ob
        ! Too less buddy numbers < 2, no buddy check completed:
        err(num)     = 0.
! Not change the flag (YRG, 02/10/20090
!         qc_flag_small(num) = qc_flag_small(num) + no_buddies
! Not change the flag (YRG, 02/10/20090
!         qc_flag_small_g(num,myproc) = qc_flag_small_g(num,myproc) + no_buddies
         i_dummy(1) = qc_flag_small_g(num,myproc)
         i_dummy(2) = num
         i_dummy(3) = myproc
         call mpi_bcast ( i_dummy, 3, mpi_integer, myproc, comm, ierr )
         qc_flag_small_g(i_dummy(2),i_dummy(3)) = i_dummy(1)

     END IF check_bad_ob

 ! 3.3 If the buddy number is ONLY 2, increase the tolerance value:
     IF ( buddy_num_final(num) .EQ. 2 ) difmax(num) = 1.2 * difmax(num)

  END DO station_loop_2

! [4.0] Reset the qc flags according to the buddy check errors:
!
!       If an observation has been flagged as failing the buddy
!       check criteria (error [err] larger than the allowable 
!       difference [difmax]), then qc_flag for this variable, 
!       level, time has to be modified.  

   station_loop_4 : do n = 1, numobs
   ! if the difference at station: n from the buddy-average > difmax,
   ! failed the buddy check:
     if ( abs(err(n)) > difmax(n) ) then
        qc_flag_small(n) = fails_buddy_check
     endif
   enddo station_loop_4

! [5.0] This notifies the user which observations have been flagged as
!       suspect by the buddy check routine.

  IF ( print_buddy ) THEN
     n    = 0
     numj = 0
     num3 = 0
     num4 = 0
     station_loop_5 : DO num = 1 , numobs

        IF ( name(1:8) .EQ. 'PMSL    ' ) THEN
           scale = 100.
        ELSE IF ( name(1:8) .EQ. 'QQ      ' ) THEN
           scale = 0.001
        ELSE
           scale = 1.0
        END IF

        IF ( buddy_num_final(num) .GE. 2 ) THEN
          n = n + 1
          IF ( ABS ( err(num) ) .GT. difmax(num) ) THEN
             numj = numj + 1
             WRITE ( unit=iunit, FMT = '(3I8,&
                   &  " ID=",a5,&    
                   &  ",NAME="  ,a8, &
                   &  ",BUDDY TOLERANCE=",f5.1,&
                   &  ",LOC=(" ,f6.1,",",f6.1,")",&
                   &  ",OBS INV="       ,f7.2,&
                   &  ",DIFF ERR="      ,f7.2,&
                   &  ",NOBS BUDDY="    ,I4,&
                   &  ",QC_FLAG=",I3,&
                   &  "  QC_FLAG_NEW=",I3,"  FAILED")' ) num, n, numj,&
                  station_id(num),name,difmax(num)/scale, &
                  xob(num),yob(num),obs(num)/scale,&
                  err(num)/scale, buddy_num_final(num), &
                  qc_flag(num), qc_flag_small(num)
          ELSE
             num3 = num3 + 1
!              WRITE ( unit=iunit, FMT = '(3I8,&
!                    &  " ID=",a5,&    
!                    &  ",NAME="  ,a8, &
!                    &  ",BUDDY TOLERANCE=",f5.1,&
!                    &  ",LOC=(" ,f6.1,",",f6.1,")",&
!                    &  ",OBS INV="       ,f7.2,&
!                    &  ",DIFF ERR="      ,f7.2,&
!                    &  ",NOBS BUDDY="    ,I4,&
!                    &  ",QC_FLAG=",I3,&
!                    &  "  QC_FLAG_NEW=",I3,"  PASSED")' ) num, n, num3,&
!                   station_id(num),name,difmax(num)/scale, &
!                   xob(num),yob(num),obs(num)/scale,&
!                   err(num)/scale, buddy_num_final(num), &
!                   qc_flag(num), qc_flag_small(num)
          ENDIF
        ELSE
          num4 = num4 + 1
!              WRITE ( unit=iunit, FMT = '(3I8,&
!                    &  " ID=",a5,&    
!                    &  ",NAME="  ,a8, &
!                    &  ",BUDDY TOLERANCE=",f5.1,&
!                    &  ",LOC=(" ,f6.1,",",f6.1,")",&
!                    &  ",OBS INV="       ,f7.2,&
!                    &  ",DIFF ERR="      ,f7.2,&
!                    &  ",NOBS BUDDY="    ,I4,&
!                    &  ",QC_FLAG=",I3,&
!                    &  "  QC_FLAG_NEW=",I3,"NO BUDDY")' ) num, n, num4,&
!                   station_id(num),name,difmax(num)/scale, &
!                   xob(num),yob(num),obs(num)/scale,&
!                   err(num)/scale, buddy_num_final(num), &
!                   qc_flag(num), qc_flag_small(num)
        END IF
     END DO station_loop_5

     write(unit=iunit,fmt = '(5x,"NOB=",i6,&
          & "  Toltal N_buddy-checked =",i6,&
          & "  N_Passed =",i6,"  N_Failed=",i6,"  NO_Buddy=",i6, &
          & "  RANGE=",f10.2,"  DIFMAX=",f10.2)') &
          numobs, n, num3, numj, num4, range, difmax(1)/scale
  END IF

END SUBROUTINE da_buddy_qc


subroutine da_eof_decomposition_test (kz, bx, e, l)
   
   !------------------------------------------------------------------------------
   ! Purpose: 
   ! [1] Print eigenvalues:
   ! [2] Test orthogonality of eigenvectors - sum_k (e_m(k) e_n(k)) = delta_mn:
   ! [3] Test eigenvectors completeness - sum_m (e_m(k1) e_m(k2)) = delta_k1k2:
   ! [4] Check B correctness: B = E*L*E^T
   !------------------------------------------------------------------------------
   
   implicit none

   integer, intent(in) :: kz               ! Dimension of BE matrix   
   real,    intent(in) :: bx(1:kz,1:kz)    ! Global vert. background error.
   real*8,  intent(in) :: e(1:kz,1:kz)     ! Eigenvectors of Bx.
   real*8,  intent(in) :: l(1:kz)          ! Eigenvalues of Bx.
   
   integer                  :: k, k1, k2, m     ! Loop counters
   real                     :: tot_variance     ! Total variance.
   real                     :: cum_variance     ! Cumulative variance.
   real                     :: max_off_diag     ! Maximum off-diagonal.

   real                     :: work(1:kz,1:kz)  ! 2D Work matrix.
   real                     :: bc(1:kz,1:kz)    ! 2D Work matrix.
   logical                  :: array_mask(1:kz) ! Array mask for MAXVAL.

   if (trace_use) call da_trace_entry("da_eof_decomposition_test")

   !------------------------------------------------------------------------- 
   ! [1] Print eigenvalues:
   !-------------------------------------------------------------------------

   tot_variance = sum(l(1:kz))
   cum_variance = 0.0
   
   write(unit=stdout,fmt='(A)')'  Mode    Eigenvalue     Cumulative Variance      e(k,k)'

   do k = 1, kz
      cum_variance = cum_variance + l(k)
      write(unit=stdout,fmt='(I4,4x,e12.4,10x,f8.4,4x,e12.4)') &
            k, l(k), cum_variance / tot_variance, e(k,k)
   end do

   write(unit=stdout,fmt=*)
   
   call da_array_print( 1, e, 'Global Eigenvectors' )

   !-------------------------------------------------------------------------
   ! [2] Test orthogonality of eigenvectors - sum_k (e_m(k) e_n(k)) = delta_mn:
   !-------------------------------------------------------------------------
   
   write(unit=stdout,fmt='(A)')' Eigenvector orthogonality check:'
   write(unit=stdout,fmt='(A)')' Mode     Diagonal         Maximum off-diagonal'

   do k1 = 1, kz
      do k2 = 1, kz
         work(k1,k2) = sum(e(1:kz,k1) * e(1:kz,k2))
      end do
   
      array_mask(1:kz) =.true.
      array_mask(k1) = .false.
      max_off_diag = maxval(abs(work(k1,:)),mask=array_mask(:))
      write(unit=stdout,fmt='(I4,4x,1pe12.4,10x,1pe12.4)')k1, work(k1,k1), max_off_diag
   end do
   write(unit=stdout,fmt=*)

   !-------------------------------------------------------------------------   
   ! [3] Test eigenvectors completeness - sum_m (e_m(k1) e_m(k2)) = delta_k1k2:
   !-------------------------------------------------------------------------   
   
   write(unit=stdout,fmt='(A)')' Eigenvector completeness check:'
   write(unit=stdout,fmt='(A)')' Level    Diagonal         Maximum off-diagonal'

   do k1 = 1, kz
      do k2 = 1, kz
         work(k1,k2) = sum(e(k1,1:kz) * e(k2,1:kz))
      end do
   
      array_mask(1:kz) =.true.
      array_mask(k1) = .false.
      max_off_diag = maxval(abs(work(k1,:)),mask=array_mask(:))
      write(unit=stdout,fmt='(I4,4x,1pe12.4,10x,1pe12.4)')k1, work(k1,k1), max_off_diag
   end do
   write(unit=stdout,fmt=*)

   !-------------------------------------------------------------------------
   ! [4]  check B correctness: B = E*L*E^T
   !-------------------------------------------------------------------------

   write(unit=stdout,fmt='(a/a)') &
        'real and Calculated B (diagonal)', &
        'lvl                 real-B                    Calculated-B'

   do k=1,kz
      do m=1,kz
         work(k,m)=l(k)*e(m,k)
         bc(k,m)=0.0
      end do
   end do
   
   do k1=1,kz
      do k2=1,kz
         do m=1,kz
            bc(k1,k2)=bc(k1,k2)+e(k1,m)*work(m,k2)
         end do
      end do

      write(unit=stdout,fmt='(I5,2F20.5)') k1, bx(k1,k1), bc(k1,k1)
   end do

   do k2=1,kz
      write(unit=stdout, fmt='(a,i4/a)') &
           'real and Calculated B (off diagonal):', k2, &
           'lvl                 real-B                    Calculated-B'

      do k1=1,kz
        write(unit=stdout,fmt='(I5,2F20.5)') k1, bx(k1,k2), bc(k1,k2)
      end do
   end do

   if (trace_use) call da_trace_exit("da_eof_decomposition_test")
   
end subroutine da_eof_decomposition_test


subroutine da_eof_decomposition (kz, bx, e, l)
   
   !---------------------------------------------------------------------------
   ! Purpose: Compute eigenvectors E and eigenvalues L of vertical covariance 
   !          matrix
   !          B_{x} defined by equation:  E^{T} B_{x} E = L, given input kz x kz 
   !          BE field.
   !---------------------------------------------------------------------------
   
   implicit none

   integer, intent(in)  :: kz               ! Dimension of error matrix. 
   real,    intent(in)  :: bx(1:kz,1:kz)    ! Vert. background error.
   real*8,  intent(out) :: e(1:kz,1:kz)     ! Eigenvectors of Bx.
   real*8,  intent(out) :: l(1:kz)          ! Eigenvalues of Bx.

   integer :: work             ! Size of work array.
   integer :: m                ! Loop counters
   integer :: info             ! Info code.

   real*8  :: work_array(1:3*kz-1)
   real*8  :: ecopy(1:kz,1:kz)
   real*8  :: lcopy(1:kz)   

   if (trace_use) call da_trace_entry("da_eof_decomposition")    

   !-------------------------------------------------------------------------
   ! [5.0]: Perform global eigenvalue decomposition using LAPACK software:
   !-------------------------------------------------------------------------
   
   work = 3 * kz - 1   
   ecopy(1:kz,1:kz) = bx(1:kz,1:kz)
   lcopy(1:kz) = 0.0

   call dsyev( 'V', 'U', kz, ecopy, kz, lcopy, work_array, work, info )
   
   if ( info /= 0 ) then
      write(unit=message(1),fmt='(A,I4)') &
         "Error in decomposition, info = ", info
      call da_error("da_eof_decomposition.inc",40,message(1:1))
   end if
   
   ! Swap order of eigenvalues, vectors so 1st is one with most variance:
   
   do m = 1, kz
      l(m) = lcopy(kz+1-m)
      e(1:kz,m) = ecopy(1:kz,kz+1-m)
   end do  

   if (trace_use) call da_trace_exit("da_eof_decomposition")    
   
end subroutine da_eof_decomposition


subroutine da_lubksb(n, np, indx, a, b)

   !-----------------------------------------------------------------------
   ! Purpose: Adapted Numerical Recipes routine to solve the set of n linear 
   ! equations A.X=B.
   ! Routine takes in to account possibility that B will begin with many zero elements, 
   ! so it is efficient for matrix inversion.
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: n              ! Logical size of array.
   integer, intent(in)    :: np             ! Physical size of array.
   integer, intent(in)    :: indx(1:n)      ! Permutation vector returned by LUDCMP. 
   real,    intent(in)    :: a(1:np,1:np)   ! LU decomposition of matrix A in A.x=B.
   real,    intent(inout) :: b(1:n)         ! On input = B, on output = x.

   integer :: i , ii , j , ll
   real    :: sum

   if (trace_use) call da_trace_entry("da_lubksb")

   ii = 0
   do i = 1 , n
      ll = indx(i)
      sum = b(ll)
      b(ll) = b(i)

      if (ii /= 0) then
         do j = ii , i - 1
            sum = sum - a(i,j) * b(j)
         end do
      else if (sum /= 0.0) then
         ii = i
      end if
      b(i) = sum
   end do

   do i = n , 1 , -1
      sum = b(i)
      if (i < n) then
         do j = i + 1 , n
            sum = sum - a(i,j) * b(j)
         end do
      end if
      b(i) = sum / a(i,i)
   end do

   if (trace_use) call da_trace_exit("da_lubksb")

end subroutine da_lubksb


subroutine da_ludcmp(n, np, indx, a, d)

   !-----------------------------------------------------------------------
   ! Purpose: Adapted Numerical Recipes routine to solve the set of n linear 
   ! equations 
   ! A.X=B. Routine takes in to account possibility that B will begin with many 
   ! zero elements, so it is efficient for matrix inversion.
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: n           ! Logical size of array.
   integer, intent(in)    :: np          ! Physical size of array.
   integer, intent(out)   :: indx(1:n)   ! Permutation vector returned by LUDCMP 
   real,    intent(inout) :: a(1:np,1:np)! LU decomposition of matrix A in A.x=B
   real,    intent(out)   :: d           ! On input = B, on output = x.

   real, parameter      :: tiny = 1.0e-20
   real                 :: aamax , dum , sum
   integer              :: i , imax , j , k
   real                 :: vv(1:np)

   if (trace_use) call da_trace_entry("da_ludcmp")

   d = 1.0
   do i = 1 , n
      aamax = 0.0

      do j = 1 , n
         if (abs(a(i,j)) > aamax) aamax = abs(a(i,j))
      end do
      if (aamax == 0.0) then
         call da_error("da_ludcmp.inc",33,(/"Singular matrix"/))
      end if
      vv(i) = 1.0 / aamax
   end do

   do j = 1 , n
      if (j > 1) then
         do i = 1 , j - 1
            sum = a(i,j)
            if (i > 1) then
               do k = 1 , i - 1
                  sum = sum - a(i,k) * a(k,j)
               end do
               a(i,j) = sum
            end if
         end do
      end if

      aamax = 0.0
      do i = j , n
         sum = a(i,j)
         if (j > 1) then
            do k = 1 , j - 1
               sum = sum - a(i,k) * a(k,j)
            end do
            a(i,j) = sum
         end if
         dum = vv(i) * abs(sum)
         if (dum >= aamax) then
            imax = i
            aamax = dum
         end if
      end do

      if (j /= imax) then
         do k = 1 , n
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
         end do
         d = -d
         vv(imax) = vv(j)
      end if

      indx(j) = imax
      if (j /= n) then
         if (a(j,j) == 0.0) a(j,j) = tiny
         dum = 1.0 / a(j,j)
         do i = j + 1 , n
            a(i,j) = a(i,j) * dum
         end do
      end if
   end do

   if (a(n,n) == 0.0) a(n,n) = tiny

   if (trace_use) call da_trace_exit("da_ludcmp")

end subroutine da_ludcmp


subroutine da_set_boundary_xa(grid)

   !------------------------------------------------------------------------
   !  Purpose: 
   !
   !  Merge East-West boundary values for the desired grid%xa-type variables
   !------------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid

   integer :: n, j, k

   if ((its /= ids) .or. (ite /= ide)) return

   if (trace_use) call da_trace_entry("da_set_boundary_xa")

   ! 2d

   k = kte + 1
   do j=jms, jme
      do n=1,bdyzone
         grid%xa%psfc(ids-n,j) = grid%xa%psfc(ide+1-n,j)
         grid%xa%w(ids-n,j,k) = grid%xa%w(ide+1-n,j,k)

         grid%xa%psfc(ide+n,j) = grid%xa%psfc(ids-1+n,j)
         grid%xa%w(ide+n,j,k) = grid%xa%w(ids-1+n,j,k)
      end do
   end do

   ! 3d

   do k=kts, kte
      do j=jms, jme
         do n=1,bdyzone
            grid%xa%u(ids-n,j,k) = grid%xa%u(ide+1-n,j,k)
            grid%xa%v(ids-n,j,k) = grid%xa%v(ide+1-n,j,k)
            grid%xa%t(ids-n,j,k) = grid%xa%t(ide+1-n,j,k)
            grid%xa%p(ids-n,j,k) = grid%xa%p(ide+1-n,j,k)
            grid%xa%q(ids-n,j,k) = grid%xa%q(ide+1-n,j,k)
            grid%xa%w(ids-n,j,k) = grid%xa%w(ide+1-n,j,k)

            grid%xa%u(ide+n,j,k) = grid%xa%u(ids-1+n,j,k)
            grid%xa%v(ide+n,j,k) = grid%xa%v(ids-1+n,j,k)
            grid%xa%t(ide+n,j,k) = grid%xa%t(ids-1+n,j,k)
            grid%xa%p(ide+n,j,k) = grid%xa%p(ids-1+n,j,k)
            grid%xa%q(ide+n,j,k) = grid%xa%q(ids-1+n,j,k)
            grid%xa%w(ide+n,j,k) = grid%xa%w(ids-1+n,j,k)
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_set_boundary_xa")

end subroutine da_set_boundary_xa


subroutine da_set_boundary_xb(grid)

   !------------------------------------------------------------------------
   ! Purpose:  
   !
   ! Merge East-West boundary values for the desired grid%xb-type variables 
   !------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid       ! first guess state.

   integer :: n, j, k

   if ((its /= ids) .or. (ite /= ide)) return

   if (trace_use) call da_trace_entry("da_set_boundary_xb")

   ! 2d
   k = kte + 1
   do j=jms, jme
      do n=1,bdyzone
         grid%xb%lat(ids-n,j)        = grid%xb%lat(ide+1-n,j)
         grid%xb%lon(ids-n,j)        = grid%xb%lon(ide+1-n,j)
         grid%xb%cori(ids-n,j)       = grid%xb%cori(ide+1-n,j)
         grid%xb%terr(ids-n,j)       = grid%xb%terr(ide+1-n,j)
         grid%xb%psfc(ids-n,j)       = grid%xb%psfc(ide+1-n,j)
         grid%xb%map_factor(ids-n,j) = grid%xb%map_factor(ide+1-n,j)
         grid%xb%coefx(ids-n,j)      = grid%xb%coefx(ide+1-n,j)
         grid%xb%coefy(ids-n,j)      = grid%xb%coefy(ide+1-n,j)
         grid%xb%w(ids-n,j,k)        = grid%xb%w(ide+1-n,j,k)
         ! grid%xb%h(ids-n,j,k)        = grid%xb%h(ide+1-n,j,k)

         grid%xb%lat(ide+n,j)        = grid%xb%lat(ids-1+n,j)
         grid%xb%lon(ide+n,j)        = grid%xb%lon(ids-1+n,j)
         grid%xb%cori(ide+n,j)       = grid%xb%cori(ids-1+n,j)
         grid%xb%terr(ide+n,j)       = grid%xb%terr(ids-1+n,j)
         grid%xb%psfc(ide+n,j)       = grid%xb%psfc(ids-1+n,j)
         grid%xb%map_factor(ide+n,j) = grid%xb%map_factor(ids-1+n,j)
         grid%xb%coefx(ide+n,j)      = grid%xb%coefx(ids-1+n,j)
         grid%xb%coefy(ide+n,j)      = grid%xb%coefy(ids-1+n,j)
         grid%xb%w(ide+n,j,k)        = grid%xb%w(ids-1+n,j,k)
         ! grid%xb%h(ide+n,j,k)        = grid%xb%h(ids-1+n,j,k)

         ! Zhiquan Liu add some RTTOV variables
         !--------------------------------------
         grid%xb%t2(ide+n,j) = grid%xb%t2(ids-1+n,j)
         grid%xb%q2(ide+n,j) = grid%xb%q2(ids-1+n,j)
         grid%xb%u10(ide+n,j) = grid%xb%u10(ids-1+n,j)
         grid%xb%v10(ide+n,j) = grid%xb%v10(ids-1+n,j)
         grid%xb%tsk(ide+n,j) = grid%xb%tsk(ids-1+n,j)
         grid%xb%tgrn(ide+n,j) = grid%xb%tgrn(ids-1+n,j)
         grid%xb%landmask(ide+n,j) = grid%xb%landmask(ids-1+n,j)
         grid%xb%snow(ide+n,j) = grid%xb%snow(ids-1+n,j)
         grid%xb%xland(ide+n,j) = grid%xb%xland(ids-1+n,j)

         grid%xb%smois(ide+n,j) = grid%xb%smois(ids-1+n,j)
         grid%xb%tslb(ide+n,j) = grid%xb%tslb(ids-1+n,j)
         grid%xb%xice(ide+n,j) = grid%xb%xice(ids-1+n,j)
         grid%xb%ivgtyp(ide+n,j) = grid%xb%ivgtyp(ids-1+n,j)
         grid%xb%isltyp(ide+n,j) = grid%xb%isltyp(ids-1+n,j)
         grid%xb%vegfra(ide+n,j) = grid%xb%vegfra(ids-1+n,j)
         grid%xb%snowh(ide+n,j) = grid%xb%snowh(ids-1+n,j)

      end do
   end do

   ! 3d
   do k=kts, kte
      do j=jms, jme
         do n=1,bdyzone
            grid%xb%h(ids-n,j,k) = grid%xb%h(ide+1-n,j,k)
            grid%xb%u(ids-n,j,k) = grid%xb%u(ide+1-n,j,k)
            grid%xb%v(ids-n,j,k) = grid%xb%v(ide+1-n,j,k)
            grid%xb%w(ids-n,j,k) = grid%xb%w(ide+1-n,j,k)
            grid%xb%t(ids-n,j,k) = grid%xb%t(ide+1-n,j,k)
            grid%xb%p(ids-n,j,k) = grid%xb%p(ide+1-n,j,k)
            grid%xb%q(ids-n,j,k) = grid%xb%q(ide+1-n,j,k)
            grid%xb%qs(ids-n,j,k) = grid%xb%qs(ide+1-n,j,k)
            grid%xb%es(ids-n,j,k) = grid%xb%es(ide+1-n,j,k)
            grid%xb%rh(ids-n,j,k) = grid%xb%rh(ide+1-n,j,k)
            grid%xb%td(ids-n,j,k) = grid%xb%td(ide+1-n,j,k)
            grid%xb%rho(ids-n,j,k)= grid%xb%rho(ide+1-n,j,k)

            grid%xb%h(ide+n,j,k) = grid%xb%h(ids-1+n,j,k)
            grid%xb%u(ide+n,j,k) = grid%xb%u(ids-1+n,j,k)
            grid%xb%v(ide+n,j,k) = grid%xb%v(ids-1+n,j,k)
            grid%xb%w(ide+n,j,k) = grid%xb%w(ids-1+n,j,k)
            grid%xb%t(ide+n,j,k) = grid%xb%t(ids-1+n,j,k)
            grid%xb%p(ide+n,j,k) = grid%xb%p(ids-1+n,j,k)
            grid%xb%q(ide+n,j,k) = grid%xb%q(ids-1+n,j,k)
            grid%xb%qs(ide+n,j,k) = grid%xb%qs(ids-1+n,j,k)
            grid%xb%es(ide+n,j,k) = grid%xb%es(ids-1+n,j,k)
            grid%xb%rh(ide+n,j,k) = grid%xb%rh(ids-1+n,j,k)
            grid%xb%td(ide+n,j,k) = grid%xb%td(ids-1+n,j,k)
            grid%xb%rho(ide+n,j,k) = grid%xb%rho(ids-1+n,j,k)
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_set_boundary_xb")

end subroutine da_set_boundary_xb


subroutine da_set_boundary_3d(var)
   !------------------------------------------------------------------------
   !  Purpose: 
   !
   !  Merge East-West boundary values for input 3d-array (var)
   !------------------------------------------------------------------------

   implicit none

   real, intent(inout) :: var(ims:ime, jms:jme, kms:kme)

   integer :: n, j, k

   if ((its /= ids) .or. (ite /= ide)) return

   if (trace_use) call da_trace_entry("da_set_boundary_3d")

   do k=kts, kte
      do j=jts, jte
         do n=1,bdyzone
            var(ids-n,j,k) = var(ide+1-n,j,k)
            var(ide+n,j,k) = var(ids-1+n,j,k)
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_set_boundary_3d")

end subroutine da_set_boundary_3d



subroutine da_get_2d_sum(var, name)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real,             intent(in) :: var(ims:ime, jms:jme)
   character(len=*), intent(in) :: name

   real :: partial, total

   if (trace_use) call da_trace_entry("da_get_2d_sum")

   partial = sum(var(its:ite,jts:jte)*var(its:ite,jts:jte))

   total = wrf_dm_sum_real( partial)

   write(unit=stdout, fmt='(3a, e24.14)') 'Square sum of <', trim(name), '>=', total

   if (trace_use) call da_trace_exit("da_get_2d_sum")

end subroutine da_get_2d_sum


subroutine da_get_3d_sum(var, name)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real,             intent(in) :: var(ims:ime, jms:jme, kms:kme)
   character(len=*), intent(in) :: name

   real :: partial, total

   if (trace_use) call da_trace_entry("da_get_3d_sum")

   partial = sum(var(its:ite,jts:jte,kts:kte)*var(its:ite,jts:jte,kts:kte))

   total = wrf_dm_sum_real (partial)

   write(unit=stdout, fmt='(3a, e24.14)') 'Square sum of <', trim(name), '>=', total

   if (trace_use) call da_trace_exit("da_get_3d_sum")

end subroutine da_get_3d_sum


subroutine da_get_print_lvl(prs,ipr)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real,          intent(in)  :: prs                    
   integer,       intent(out) :: ipr                    

   integer       :: k                      

   if (trace_use) call da_trace_entry("da_get_print_lvl")
   ipr = 1
   do k =2, npres_print
 ! Avoid the pressure is unavailable (ex. -8888.88) Shu-Ya Chen (Dec. 2011)
     if( (prs*.01 .ge. pptop(k-1)) .or. (prs < 0.)) exit
     ipr = k
   end do

   if (trace_use) call da_trace_exit("da_get_print_lvl")

end subroutine da_get_print_lvl



subroutine da_get_julian_time(year,month,day,hour,minute,gstime)

   !------------------------------------------------------------------------------
   ! Purpose: Calculate Julian time from year/month/day/hour/minute.
   !------------------------------------------------------------------------------

   implicit none

   integer, intent(in)  :: year
   integer, intent(in)  :: month
   integer, intent(in)  :: day
   integer, intent(in)  :: hour
   integer, intent(in)  :: minute
   real*8,  intent(out) :: gstime

   integer    :: iw3jdn, ndays, nmind

   if (trace_use) call da_trace_entry("da_get_julian_time")

   iw3jdn  =    day - 32075 &
              + 1461 * (year + 4800 + (month - 14) / 12) / 4 &
              + 367 * (month - 2 - (month - 14) / 12 * 12) / 12 &
              - 3 * ((year + 4900 + (month - 14) / 12) / 100) / 4
   ndays = iw3jdn - 2443510

   nmind = ndays*1440 + hour * 60 + minute
   gstime = float(nmind)

   if (trace_use) call da_trace_exit("da_get_julian_time")

end subroutine da_get_julian_time
subroutine da_get_time_slots(nt,tmin,tana,tmax,time_slots)

   !------------------------------------------------------------------------------
   ! Purpose: Calculate time slots for FGAT option.
   !------------------------------------------------------------------------------

   implicit none

   integer,           intent (in)      :: nt    ! number of time slots
   character(len=19), intent (in)      :: tmin  ! begin of time window
   character(len=19), intent (in)      :: tana  ! center of time window
   character(len=19), intent (in)      :: tmax  ! end of time window
   real*8,           intent (out)      :: time_slots(0:nt) !

   integer   :: min_yyyy,min_mm,min_dd,min_hh,min_mn,min_ss
   integer   :: max_yyyy,max_mm,max_dd,max_hh,max_mn,max_ss
   character :: s
   real      :: dt
   integer   :: it

   if (trace_use) call da_trace_entry("da_get_time_slots")

   read(unit=tmin,fmt='(i4,5(a1,i2))') min_yyyy,s,min_mm,s,min_dd,s,min_hh,s,min_mn,s,min_ss
   read(unit=tmax,fmt='(i4,5(a1,i2))') max_yyyy,s,max_mm,s,max_dd,s,max_hh,s,max_mn,s,max_ss

   if (print_detail_obs) then
      write(unit=stdout,fmt='(3X,A,I4,5(1X,I2))') 'Start julian time : ',min_yyyy,min_mm,min_dd,min_hh,min_mn,min_ss
      write(unit=stdout,fmt='(3X,A,I4,5(1X,I2))') 'End julian time   : ',max_yyyy,max_mm,max_dd,max_hh,max_mn,max_ss
   end if

   call da_get_julian_time(min_yyyy,min_mm,min_dd,min_hh,min_mn,time_slots(0))
   call da_get_julian_time(max_yyyy,max_mm,max_dd,max_hh,max_mn,time_slots(nt))

   if (nt > 1) then
      dt = (time_slots(nt)-time_slots(0))/float(nt-1)
      time_slots(1)  = time_slots(0)+dt*0.5
      do it=2,nt-1
         time_slots(it) = time_slots(it-1)+dt
      end do
   end if

   if (print_detail_obs) then
      write(unit=stdout,fmt='(3x,a,240f10.0)') 'Time_slots ', time_slots(0:nt)
      write (unit=stdout,fmt='(A)') " "
   end if

   if (trace_use) call da_trace_exit("da_get_time_slots")

end subroutine da_get_time_slots



!
subroutine da_msl2geo1 (z, lat, lon, h)
      implicit none
      real*8 h, lat, lon, z     ! Note lon is currently not used
      real*8 pi, latr, zm
      real*8 semi_major_axis, semi_minor_axis, grav_polar, grav_equator
      real*8 earth_omega, grav_constant, flattening, somigliana
      real*8 grav_ratio, sin2, termg, termr, grav, eccentricity

!     Parameters below from WGS-84 model software inside GPS receivers.
      parameter(semi_major_axis = 6378.1370d3)    ! (m)
      parameter(semi_minor_axis = 6356.7523142d3) ! (m)
      parameter(grav_polar = 9.8321849378)        ! (m/s2)
      parameter(grav_equator = 9.7803253359)      ! (m/s2)
      parameter(earth_omega = 7.292115d-5)        ! (rad/s)
      parameter(grav_constant = 3.986004418d14)   ! (m3/s2)
      parameter(grav = 9.80665d0)                 ! (m/s2) WMO std g at 45 deg lat
      parameter(eccentricity = 0.081819d0)        ! unitless

!     Derived geophysical constants
      parameter(flattening = (semi_major_axis-semi_minor_axis) /  &
                              semi_major_axis)

      parameter(somigliana =  &
       (semi_minor_axis/semi_major_axis)*(grav_polar/grav_equator)-1.d0)

      parameter(grav_ratio = (earth_omega*earth_omega *  &
       semi_major_axis*semi_major_axis * semi_minor_axis)/grav_constant)

      pi   = 3.14159265358979d0
      latr = lat * (pi/180.d0)        ! in radians
      zm   = z*1000d0                 ! in meters

      if (z.eq.-999.d0) then
         h = -999.d0
      else 
         sin2  = sin(latr) * sin(latr)
         termg = grav_equator *  &
          ( (1.d0+somigliana*sin2) /  &
            sqrt(1.d0-eccentricity*eccentricity*sin2) )
         termr = semi_major_axis   /  &
          (1.d0 + flattening + grav_ratio - 2.d0*flattening*sin2)

!     geopotential height 

         h=((termg/grav)*((termr*zm)/(termr+zm)))*0.001 ! in km

      endif
  
end  subroutine da_msl2geo1 
!
subroutine da_geo2msl1 (h, lat, lon, z)
      implicit none
      real*8 h, lat, lon, z     ! Note lon is currently not used
      real*8 pi, latr, hm
      real*8 semi_major_axis, semi_minor_axis, grav_polar, grav_equator
      real*8 earth_omega, grav_constant, flattening, somigliana
      real*8 grav_ratio, sin2, termg, termr, grav, eccentricity

!     Parameters below from WGS-84 model software inside GPS receivers.
      parameter(semi_major_axis = 6378.1370d3)    ! (m)
      parameter(semi_minor_axis = 6356.7523142d3) ! (m)
      parameter(grav_polar = 9.8321849378)        ! (m/s2)
      parameter(grav_equator = 9.7803253359)      ! (m/s2)
      parameter(earth_omega = 7.292115d-5)        ! (rad/s)
      parameter(grav_constant = 3.986004418d14)   ! (m3/s2)
      parameter(grav = 9.80665d0)                 ! (m/s2) WMO std g at 45 deg lat
      parameter(eccentricity = 0.081819d0)        ! unitless

!     Derived geophysical constants
      parameter(flattening = (semi_major_axis-semi_minor_axis) /  &
                              semi_major_axis)

      parameter(somigliana = &
       (semi_minor_axis/semi_major_axis)*(grav_polar/grav_equator)-1.d0)

      parameter(grav_ratio = (earth_omega*earth_omega * &
       semi_major_axis*semi_major_axis * semi_minor_axis)/grav_constant)

      pi   = 3.14159265358979d0
      latr = lat * (pi/180.d0)        ! in radians
      hm   = h * 1000d0               ! in meters

      if (h.eq.-999.d0) then
         z = -999.d0
      else 
         sin2  = sin(latr) * sin(latr)
         termg = grav_equator *  &
          ( (1.d0+somigliana*sin2) /  &
            sqrt(1.d0-eccentricity*eccentricity*sin2) )
         termr = semi_major_axis   /  &
          (1.d0 + flattening + grav_ratio - 2.d0*flattening*sin2)

!     geometric height 

         z = ( (termr*hm) / ( (termg/grav) * termr - hm ) ) * 0.001 ! in km

      endif
  
end subroutine da_geo2msl1 

end module da_tools

