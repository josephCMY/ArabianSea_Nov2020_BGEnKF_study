












module gsi_constants




























  use gsi_kinds, only: r_single,r_kind,i_kind
  implicit none


  integer(i_kind) izero,ione
  real(r_kind) rearth,grav,omega,rd,rv,cp,cv,cvap,cliq
  real(r_kind) csol,hvap,hfus,psat,t0c,ttp,jcal,cp_mass
  real(r_kind) fv,deg2rad,rad2deg,pi,tiny_r_kind,huge_r_kind
  real(r_kind) ozcon,rozcon,tpwcon,rd_over_g,rd_over_cp,g_over_rd
  real(r_kind) amsua_clw_d1,amsua_clw_d2,constoz,zero,one,two,four
  real(r_kind) one_tenth,quarter,three,five
  real(r_kind) rearth_equator,stndrd_atmos_ps
  real(r_kind) semi_major_axis,semi_minor_axis,n_a,n_b
  real(r_kind) eccentricity,grav_polar,grav_ratio
  real(r_kind) grav_equator,earth_omega,grav_constant
  real(r_kind) flattening,eccentricity_linear,somigliana
  real(r_kind) dldt,dldti,hsub,psatk,tmix,xa,xai,xb,xbi
  real(r_kind) eps,epsm1,omeps
  real(r_kind) elocp,cpr,el2orc,cclimit,climit,epsq
  real(r_kind) pcpeff0,pcpeff1,pcpeff2,pcpeff3,rcp,c0,delta
  real(r_kind) h1000,factor1,factor2,rhcbot,rhctop,dx_max,dx_min,dx_inv
  real(r_kind) h300,half,cmr,cws,ke2,row,rrow
  real(r_single) zero_single,tiny_single,huge_single





  parameter(rearth_equator= 6.37813662e6_r_kind) 
  parameter(omega  = 7.2921e-5_r_kind)  
  parameter(cp     = 1.0046e+3_r_kind)  
  parameter(cvap   = 1.8460e+3_r_kind)  
  parameter(csol   = 2.1060e+3_r_kind)  
  parameter(hvap   = 2.5000e+6_r_kind)  
  parameter(hfus   = 3.3358e+5_r_kind)  
  parameter(psat   = 6.1078e+2_r_kind)  
  parameter(t0c    = 2.7315e+2_r_kind)  
  parameter(ttp    = 2.7316e+2_r_kind)  
  parameter(jcal   = 4.1855e+0_r_kind)  
  parameter(stndrd_atmos_ps = 1013.25e2_r_kind) 


  parameter(izero  = 0)
  parameter(ione   = 1)
  parameter(zero_single = 0.0_r_single)
  parameter(zero   = 0.0_r_kind)
  parameter(one_tenth  = 0.10_r_kind)
  parameter(quarter= 0.25_r_kind)
  parameter(one    = 1.0_r_kind)
  parameter(two    = 2.0_r_kind)
  parameter(three  = 3.0_r_kind)
  parameter(four   = 4.0_r_kind)
  parameter(five   = 5.0_r_kind)


  parameter(n_a=77.6_r_kind) 
  parameter(n_b=3.73e+5_r_kind) 


  parameter(semi_major_axis = 6378.1370e3_r_kind)    
  parameter(semi_minor_axis = 6356.7523142e3_r_kind) 
  parameter(grav_polar = 9.8321849378_r_kind)        
  parameter(grav_equator = 9.7803253359_r_kind)      
  parameter(earth_omega = 7.292115e-5_r_kind)        
  parameter(grav_constant = 3.986004418e14_r_kind)   


  parameter(flattening = (semi_major_axis-semi_minor_axis)/semi_major_axis)
  parameter(somigliana = &
       (semi_minor_axis/semi_major_axis) * (grav_polar/grav_equator) - one)
  parameter(grav_ratio = (earth_omega*earth_omega * &
       semi_major_axis*semi_major_axis * semi_minor_axis) / grav_constant) 


  parameter ( dldti = cvap-csol )
  parameter ( hsub = hvap+hfus )
  parameter ( psatk = psat*0.001_r_kind )
  parameter ( tmix = ttp-20._r_kind )
  parameter ( elocp = hvap/cp )
  parameter ( rcp  = one/cp )


  parameter ( h300 = 300._r_kind )
  parameter ( half = 0.5_r_kind )
  parameter ( cclimit = 0.001_r_kind )
  parameter ( climit = 1.e-20_r_kind)
  parameter ( epsq = 2.e-12_r_kind )
  parameter ( h1000 = 1000.0_r_kind)
  parameter ( rhcbot=0.85_r_kind )
  parameter ( rhctop=0.85_r_kind )
  parameter ( dx_max=-8.8818363_r_kind )
  parameter ( dx_min=-5.2574954_r_kind )
  parameter ( dx_inv=one/(dx_max-dx_min) )
  parameter ( c0=0.002_r_kind )
  parameter ( delta=0.6077338_r_kind )
  parameter ( pcpeff0=1.591_r_kind )
  parameter ( pcpeff1=-0.639_r_kind )
  parameter ( pcpeff2=0.0953_r_kind )
  parameter ( pcpeff3=-0.00496_r_kind )
  parameter ( cmr = one/0.0003_r_kind )
  parameter ( cws = 0.025_r_kind )
  parameter ( ke2 = 0.00002_r_kind )
  parameter ( row = 1000._r_kind )
  parameter ( rrow = one/row )


  parameter ( constoz = 604229.0_r_kind)



  parameter ( amsua_clw_d1 = 0.754_r_kind )
  parameter ( amsua_clw_d2 = -2.265_r_kind )

contains
  subroutine init_constants_derived
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_constants_derived          set derived constants
!     prgmmr:    treadon          org: np23           date: 2004-12-02
!
! abstract:  This routine sets derived constants
!
! program history log:
!   2004-12-02  treadon
!   2005-03-03  treadon - add implicit none
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!$$$
    implicit none

!   Trigonometric constants
    pi      = acos(-one)
    deg2rad = pi/180.0_r_kind
    rad2deg = one/deg2rad
    tiny_r_kind = tiny(zero)
    huge_r_kind = huge(zero)
    tiny_single = tiny(zero_single)
    huge_single = huge(zero_single)

!   Geophysical parameters used in conversion of geopotential to
!   geometric height
    eccentricity_linear = sqrt(semi_major_axis**2 - semi_minor_axis**2)
    eccentricity = eccentricity_linear / semi_major_axis

    return
  end subroutine init_constants_derived



end module gsi_constants
