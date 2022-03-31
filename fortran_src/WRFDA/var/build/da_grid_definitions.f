












module da_grid_definitions
   
   
   
   
   
   use da_control, only : gravity, convert_fd2uv, gas_constant, convert_uv2fd, &
      pi, cone_factor, map_projection, phic,psi1,earth_radius,c2,pole, trace_use, &
      ycntr,truelat1_3dv,xlonc,ptop,t0,base_pres,base_lapse,base_temp, trace_use_dull, &
      trace_use_frequent,missing_r,missing_data
   use da_reporting, only : da_error,message
   use da_tracing, only : da_trace_entry, da_trace_exit
   
   implicit none
   
   contains
   
subroutine da_ref_height(pres, ref_height)

   !---------------------------------------------------------------------------
   ! Purpose: To calculate the height from the reference pressure
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)  :: pres
   real, intent(out) :: ref_height

   real              :: aa, bb, cc
   real              :: p0iso, ziso

   real, parameter   :: rovg = gas_constant/gravity

   if (trace_use) call da_trace_entry("da_ref_height")

   aa = 0.5 * rovg * base_lapse
   bb = rovg * t0

   cc = log(pres/base_pres)
   ref_height = -(bb + aa * cc) * cc

   if (base_temp > 0.0) then
      p0iso=base_pres*exp((base_temp-t0)/base_lapse)
      cc=log(p0iso/base_pres)
      ziso=-(aa*cc+bb)*cc

      if (ref_height > ziso) then
         ref_height=ziso+rovg*base_temp*log(p0iso/pres)
      end if
   end if

   if (trace_use) call da_trace_exit("da_ref_height")

end subroutine da_ref_height


subroutine da_ref_pres(height, ref_pres)

   !---------------------------------------------------------------------------
   ! Purpose: To calculate the reference pressure from the height.
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)  :: height
   real, intent(out) :: ref_pres

   real, parameter :: rovg = gas_constant/gravity

   real :: aa, bb, cc, hh, htop, ziso, p0iso

   if (trace_use_frequent) call da_trace_entry("da_ref_pres")

   call da_ref_height(ptop, htop)

   bb = rovg * t0
   aa = rovg * base_lapse * 0.5

   hh = min(htop, height)
   cc = (-bb + sqrt(bb * bb - 4.0 * aa * hh))/(2.0*aa)
   ref_pres = base_pres * exp(cc)

   if (base_temp > 0.0) then
      p0iso=base_pres*exp((base_temp-t0)/base_lapse)
      cc=log(p0iso/base_pres)
      ziso=-(aa*cc+bb)*cc

      if (height > ziso) then
         ref_pres = p0iso/exp((height-ziso)/(rovg*base_temp))
      end if
   end if

   if (trace_use_frequent) call da_trace_exit("da_ref_pres")

end subroutine da_ref_pres


subroutine da_ffdduv (F,D,U,V,YLON,ID)

   !-------------------------------------------------------------------------
   ! Purpose: TBD
   ! When ID =  1
   ! Convert wind speed (F in m/s) and direction (D in degree 0-360) into
   ! wind (U-V in m/s) components
   !
   ! When ID = -1
   ! Convert wind (U-V in m/s) components into wind speed (F in m/s) and 
   ! direction (D in degree 0-360)
   !
   ! Need map projection parameters from module da_control
   !
   ! PHIC:  Central latitude 
   ! XLONC: Central longitude
   ! XN:    Cone projection
   ! CONV:  180/Pi
   !
   !-------------------------------------------------------------------------

   implicit none

   real,    intent (inout) :: f,d
   real,    intent (inout) :: u, v
   real,    intent (in)    :: ylon
   integer, intent (in)    :: id

   real :: aearth, uearth, vearth
   real :: xlonrt, ang, conv

   if (trace_use_frequent) call da_trace_entry("da_ffdduv")

   conv = 180.0 / pi

   select case (ID)

      case (convert_fd2uv);

         ! convert wind module/direction into u/v wind components on earth,
         ! then convert u/v wind components on earth into lambert conformal or
         ! polar stereographic projection u/v wind components.

         ! projections change requires only a change of the cone constant, xn
         ! equations remain the same.

         AEARTH = D/CONV

         UEARTH = -F*Sin(AEARTH)
         VEARTH = -F*COS(AEARTH)

         ! for conversion to grid coordinates,
         ! see program datamap, subr vect, and
         ! ANTHES METEO. 597 NOTES, EQUA. 2.23, 2.25, 2.28.

         XLONRT = XLONC-YLON

         if (XLONRT .GT. 180.0) XLONRT=XLONRT-360.0
         if (XLONRT .LT.-180.0) XLONRT=XLONRT+360.0

         ANG=XLONRT*CONE_FACTOR/CONV

         ! for mercator projection, the winds are as in earth coordinates

         if (map_projection.EQ.3) ANG=0.0

         if (PHIC.LT.0.0) ANG=-ANG

         U = VEARTH*Sin(ANG) + UEARTH*COS(ANG)
         V = VEARTH*COS(ANG) - UEARTH*Sin(ANG)


         ! CONVERT LAMBERT CONFORMAL OR POLAR STEREOGRAPHIC PROJECTION U/V
         ! WinD COMPONENTS inTO U/V WinD COMPONENTS ON EART
         ! then CONVERT U/V WinD COMPONENTS ON EARTH inTO WinD module/DIRECTION

         ! PROJECTIONS CHANGE REQUIRES ONLY A CHANGE OF THE CONE_FACTOR

      case (convert_uv2fd);

         XLONRT = XLONC-YLON

         if (XLONRT .GT. 180.0) XLONRT=XLONRT-360.0
         if (XLONRT .LT.-180.0) XLONRT=XLONRT+360.0

         ANG=XLONRT*CONE_FACTOR/CONV

         ! FOR MERCATOR PROJECTION, THE WinDS ARE AS in EARTH COORDinATES

         if (map_projection .EQ.  3) ANG = 0.0
         if (PHIC  .LT. 0.0) ANG = -ANG

         UEARTH = U*COS(ANG) - V*Sin(ANG)
         VEARTH = U*Sin(ANG) + V*COS(ANG)

         F = sqrt(UEARTH*UEARTH + VEARTH*VEARTH)

         if (F .EQ. 0.0) then
            D = 0.0
            if (trace_use_frequent) call da_trace_exit("da_ffdduv")
            return
         end if

         if (VEARTH .EQ. 0.0) then
            if (UEARTH .GT. 0.0) D = 270.0
            if (UEARTH .LT. 0.0) D =  90.0
         else
            AEARTH = ATAN (UEARTH/VEARTH)*CONV

            if (UEARTH .LE. 0.0 .AND. VEARTH .LE. 0.0) D = AEARTH
            if (UEARTH .LE. 0.0 .AND. VEARTH .GE. 0.0) D = AEARTH + 180.0
            if (UEARTH .GE. 0.0 .AND. VEARTH .GE. 0.0) D = AEARTH + 180.0
            if (UEARTH .GE. 0.0 .AND. VEARTH .LE. 0.0) D = AEARTH + 360.0

         end if

      case default
         write(unit=message(1),fmt='(A,I2)') ' UNKNOWN OPTION ',ID
         call da_error("da_ffdduv.inc",119,message(1:1))

   end select

   if (trace_use_frequent) call da_trace_exit("da_ffdduv")

end subroutine da_ffdduv


subroutine da_ffdduv_model (F,D,U,V,ID)

   !-------------------------------------------------------------------------
   ! Purpose: TBD

   ! When ID =  1
   ! Convert wind speed (F in m/s) and direction (D in degree 0-360) into
   ! wind (U-V in m/s) components on grid coordinates
   !
   ! When ID = -1
   ! Convert wind (U-V in m/s) components into wind speed (F in m/s) and 
   ! direction (D in degree 0-360) on grid coordinates 
   !
   !-------------------------------------------------------------------------

   implicit none

   real,    intent (inout) :: F,D
   real,    intent (inout) :: U,V
   integer, intent (in)    :: ID
   real                    :: CONV,A

   if (trace_use_frequent) call da_trace_entry("da_ffdduv_model")

   CONV = 180.0 / pi

   select case (ID)

      case (convert_fd2uv);

         U = -F*SIN(D/CONV)
         V = -F*COS(D/CONV)

      case (convert_uv2fd);

         F = sqrt(U*U + V*V)

         if (F .EQ. 0.0) then
            D = 0.0
            if (trace_use_frequent) call da_trace_exit("da_ffdduv_model")
            return
         end if

         if (V .EQ. 0.0) then
            if (U .GT. 0.0) D = 270.0
            if (U .LT. 0.0) D =  90.0
         else
            A = ATAN (U/V)*CONV

            if (U .LE. 0.0 .AND. V .LE. 0.0) D = A
            if (U .LE. 0.0 .AND. V .GE. 0.0) D = A + 180.0
            if (U .GE. 0.0 .AND. V .GE. 0.0) D = A + 180.0
            if (U .GE. 0.0 .AND. V .LE. 0.0) D = A + 360.0

         end if

      case default
         write(unit=message(1),fmt='(A,I2)') ' UNKNOWN OPTION ',ID
         call da_error("da_ffdduv_model.inc",59,message(1:1))

   end select

   if (trace_use_frequent) call da_trace_exit("da_ffdduv_model")

end subroutine da_ffdduv_model


subroutine da_earth_2_model_wind(eu,ev,mu,mv,lon)

   !---------------------------------------------------------------------------
   ! Purpose: Convert earth wind to model wind.
   !
   ! Need map projection parameters.
   !
   ! IPROJ: Projection type
   ! PHIC:  Central latitude 
   ! XLONC: Central longitude
   ! XN:    Cone projection
   ! CONV:  180/Pi
   !---------------------------------------------------------------------------

   implicit none

   real*8, intent(in)  :: eu, ev
   real, intent(out) :: mu, mv
   real, intent(in)  :: lon

   real :: xlonrt, ang

   if (trace_use) call da_trace_entry("da_earth_2_model_wind")

   ! for mercator projection, the winds are as in earth coordinates

   if (map_projection == 3) then
      mu = eu
      mv = ev
      if (trace_use) call da_trace_exit("da_earth_2_model_wind")
      return
   end if

   ! for conversion to grid coordinates,
   ! see program datamap, subr vect, and
   ! ANTHES METEO. 597 NOTES, EQUA. 2.23, 2.25, 2.28.

   xlonrt = xlonc-lon

   if (xlonrt > 180.0) xlonrt=xlonrt-360.0
   if (xlonrt <-180.0) xlonrt=xlonrt+360.0

   ang=xlonrt*cone_factor*pi/180.0

   if (phic < 0.0) ang=-ang

   mu = ev*sin(ang) + eu*cos(ang)
   mv = ev*cos(ang) - eu*sin(ang)

   if (trace_use) call da_trace_exit("da_earth_2_model_wind")

end subroutine da_earth_2_model_wind

subroutine da_ffdduv_diagnose(obs1, inv1, inc1, obs2, inv2, inc2, qc1, qc2, ID) 

   implicit none

   integer ,intent (in)    :: qc1, qc2, ID
   real    ,intent (inout) :: obs1, inv1, inc1, obs2, inv2, inc2

   real   :: uobs, uinv, uinc, uan, ufg, &
             vobs, vinv, vinc, van, vfg, & 
             spdobs, spdinv, spdinc, spdfg, spdan, &
             dirobs, dirinv, dirinc, dirfg, diran 

   if (trace_use) call da_trace_entry("da_ffdduv_diagnose")

   if (obs1 .gt. missing_r .and. obs2 .gt. missing_r) then

      select case (ID)

      case (convert_fd2uv);

         spdobs = obs1
         spdinv = inv1
         spdinc = inc1
         dirobs = obs2
         dirinv = inv2
         dirinc = inc2

         spdfg  = spdobs - spdinv
         spdan  = spdobs - spdinc
         dirfg  = dirobs - dirinv
         diran  = dirobs - dirinc

         call da_ffdduv_model(spdobs, dirobs, uobs, vobs, convert_fd2uv)
         call da_ffdduv_model(spdfg,  dirfg,  ufg,  vfg,  convert_fd2uv)
         call da_ffdduv_model(spdan,  diran,  uan,  van,  convert_fd2uv)

         uinv = uobs - ufg
         uinc = uobs - uan
         vinv = vobs - vfg
         vinc = vobs - van

         obs1 = uobs
         inv1 = uinv
         inc1 = uinc
         obs2 = vobs
         inv2 = vinv
         inc2 = vinc

      case (convert_uv2fd);

         uobs = obs1
         uinv = inv1
         uinc = inc1
         vobs = obs2
         vinv = inv2
         vinc = inc2

         ufg  = uobs - uinv
         uan  = uobs - uinc
         vfg  = vobs - vinv
         van  = vobs - vinc

         call da_ffdduv_model(spdobs, dirobs, uobs, vobs, convert_uv2fd)
         call da_ffdduv_model(spdfg,  dirfg,  ufg,  vfg,  convert_uv2fd)
         call da_ffdduv_model(spdan,  diran,  uan,  van,  convert_uv2fd)

         spdinv = spdobs - spdfg
         spdinc = spdobs - spdan
         dirinv = dirobs - dirfg
         dirinc = dirobs - diran
             
         if (dirinv > 180.0) dirinv = dirinv - 360.0
         if (dirinc > 180.0) dirinc = dirinc - 360.0
         if (dirinv < -180.0) dirinv = dirinv + 360.0
         if (dirinc < -180.0) dirinc = dirinc + 360.0

         obs1 = spdobs
         inv1 = spdinv
         inc1 = spdinc
         obs2 = dirobs
         inv2 = dirinv
         inc2 = dirinc

      case default
         write(unit=message(1),fmt='(A,I2)') ' UNKNOWN OPTION ',ID
         call da_error("da_ffdduv_diagnose.inc",86,message(1:1))

      end select

   else

      obs1 = missing_r
      inv1 = 0.0
      inc1 = 0.0
      obs2 = missing_r
      inv2 = 0.0
      inc2 = 0.0

   end if

   if (trace_use) call da_trace_exit("da_ffdduv_diagnose")

end subroutine da_ffdduv_diagnose

end module da_grid_definitions
