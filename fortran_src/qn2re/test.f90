program test
  use qn2re_wsm6
  use qn2re_gfdlfv3

  real :: q1 = 1e-4

  real :: q2(10,5)
  integer :: i, j

  real :: Re(10,5)
  real :: rhoi(10,5)
  real :: mass(10,5)
  real :: Ni(10,5)

  print *, qn2re_WSM6_cloud(q1, 1.)
  print *, qn2re_WSM6_ice(q1, 1.)
  print *, qn2re_WSM6_Rain(q1, 1.)
  print *, qn2re_WSM6_Snow(q1, 1., 273.15)
  print *, qn2re_WSM6_Graup(q1, 1.)

  
  do i =1, 10
    do j=1, 5
      q2(i,j) = 1e-10 * j * 10**i * 1.00001 * 2
    end do
  end do

  Ni = 5.38e7 * q2**0.75
  RE  = qn2re_WSM6_ice_unpacked(q2, 1.)
  rhoi = qn2re_WSM6_ice_rhoi(q2, 1.)
  mass = 4./3.*3.1415926 * Re**3 * rhoi 
  print *, '---rhoaqi---------'
  do i =1, 10
!    print *, q2(i,:)
  end do
  print *, '---Ni---------'
  do i =1, 10
!    print *, Ni(i,:)
  end do
  print *, '------------'
  do i =1, 10
!    print *, Re(i,:)
  end do
  print *, '------------'
  do i =1, 10
!    print *, rhoi(i,:)
  end do
  print *, '----mass-------'
  do i =1, 10
!    print *, mass(i,:)
  end do

  print *, '---- GFDL -------'
  print *, qn2re_GFDLFV3_cloud(q1, 1.)
  print *, qn2re_GFDLFV3_ice  (q1, 1., 250.)
  print *, qn2re_GFDLFV3_rain (q1, 1.)
  print *, qn2re_GFDLFV3_snow (q1, 1.)
  print *, qn2re_GFDLFV3_graup(q1, 1.)


  print *, '---- example ----'
  print *, 'should output 58, 105, 96, 121'
  print *, nint (qn2re_GFDLFV3_ice( 0.001e-3, 1., 273.16-60)*1e6)
  print *, nint (qn2re_GFDLFV3_ice( 0.001e-3, 1., 273.16-20)*1e6)
  print *, nint (qn2re_GFDLFV3_ice( 0.1e-3, 1., 273.16-60)*1e6)
  print *, nint (qn2re_GFDLFV3_ice( 0.1e-3, 1., 273.16-20)*1e6)

end program test
