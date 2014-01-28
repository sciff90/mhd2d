subroutine iterate()
  use global
  implicit none

  !Initialize Arrays
  complex*16,dimension(0:num_u1-1,0:num_u3-1) :: e1,e2,e3,b1,b2,b3
  complex*16,dimension(0:num_u1-1,0:num_u3-1) :: esup1,esup2,bsup1,bsup2,bsup3
  double precision,dimension(0:num_u1-1,0:num_u3-1) :: cfac
  complex*16,dimension(0:num_u1-1) :: odds_u1,evens_u1
  complex*16,dimension(0:num_u3-1) :: odds_u3,evens_u3
  double precision,dimension(0:num_u3-1) :: u3,driver
  integer,dimension(0:num_u1/2) :: Eidx
  integer,dimension(0:num_u1/2-1) :: Oidx

  complex*16,dimension(0:num_u1-1,0:num_u3-1) :: e1p3,e1m3,e2p1,e2m1,e2p3,&
    e2m3,e3p1,e3m1,av_bsup1,av_bsup3
  complex*16,dimension(0:num_u1-1,0:num_u3-1) :: b1p3,b1m3,b2p1,b2m1,b2p3,&
                                                 b2m3,b3p1,b3m1
  !Counters
  integer*4 :: ii,jj,kk,tpnext
  integer*4 :: iplot,nplots
  double precision :: tt


  do ii = 0, num_u3/2-1
  odds_u3(2*ii+1) = 1.0
  evens_u3(2*ii) = 1.0
  end do
  evens_u3(num_u3-1) = 1.0

  do ii = 0, num_u1/2-1
  odds_u1(2*ii+1) = 1.0
  evens_u1(2*ii) = 1.0
  end do
  evens_u1(0) = 0.0
  evens_u1(num_u1-1) = 0.0

  do ii = 0, num_u1/2
  Eidx(ii) = ii*2
  end do
  do ii = 0, num_u1/2-1
  Oidx(ii) = ii*2+1
  end do
  do ii = 0, num_u3-1
  u3(ii) = dble(ii)
  end do

  cfac = dt*sqrt(va2_arr)*0.5/h3
  nplots = tmax/dtplot
  tpnext = dtplot
  iplot = 0
  tt=0.0

  write(*,*)'Starting Time Loop'
  do while (iplot .lt. nplots)
    tt = tt+dt
    call get_driver(tt,u3,driver)

    !shift ei's
    e1p3 = cshift(e1,shift=1,dim=2)
    e1m3 = cshift(e1,shift=-1,dim=2)
    e1p3(:,num_u3-1) = 0.0
    e1m3(:,0) = 0.0

    e2p1 = cshift(e2,shift=1,dim=1)
    e2m1 = cshift(e2,shift=-1,dim=1)
    e2p1(num_u1-1,:) = 0.0
    e2m1(0,:) = 0.0

    e2p3 = cshift(e1,shift=1,dim=2)
    e2m3 = cshift(e1,shift=-1,dim=2)
    e2p3(:,num_u3-1) = 0.0
    e2m3(:,0) = 0.0

    e3p1 = cshift(e1,shift=1,dim=2)
    e3m1 = cshift(e1,shift=-1,dim=2)
    e3p1(num_u1-1,:) = 0.0
    e3m1(0,:) = 0.0

    !Advance bsupi
    bsup1 = bsup1 + b1e2*(e2p3-e2m3) + b1e3*e3
    bsup2 = bsup2 + b2e1*(e1p3-e1m3) + b2e3*(e3p1-e3m1)
    bsup3 = bsup3 + b3e1*e1 + b3e2*(e2p1-e2m1)


  !  !Inner L shell Perfectly reflecting
    !if (ixbc_0 eq 0) then
      !bsup3(0,*) = 0.d0
      !bsup3(0,0) = 0.d0
      !bsup3(0,num_u3-1) = 0.d0
    !else
      !bsup2(0,*) = (4.0*bsup2(2,*) - bsup2(4,*))/3.0
      !bsup3(0,*) = (4.0*bsup3(2,*) - bsup3(4,*))/3.0
      !!Corners
      bsup3(0,0) = (4.0*bsup3(2,0) - bsup3(4,0))/3.0
      bsup3(0,num_u3-1) = (4.0*bsup3(2,num_u3-1) - bsup3(4,num_u3-1))/3.0
    !end if

    !Outer L shell Perfectly reflecting???????
  !  if (ixbc_n eq 0) then
      !bsup3[num_u1-1,*] = 0.d0
      !bsup3[num_u1-1,0] = 0.d0
      !bsup3[num_u1-1,num_u3-1]= 0.d0
    !else
      !bsup2[num_u1-1,*] = (4.0*bsup2(num_u1-3,*) - bsup2(num_u1-5,*))/3.0
      !bsup3[num_u1-1,*] = (4.0*bsup3(num_u1-3,*) - bsup3(num_u1-5,*))/3.0
      !!Corners
      bsup3(num_u1-1,0) = (4.0*bsup3(num_u1-3,0) - bsup3(num_u1-5,0))/3.0
      bsup3(num_u1-1,num_u3-1)= (4.0*bsup3(num_u1-3,num_u3-1) &
                                - bsup3(num_u1-5,num_u3-1))/3.0
    !end if

    !average bsups
    av_bsup1 = (cshift(cshift(bsup1,shift=-1,dim=2),shift=-1,dim=1) + &
                cshift(cshift(bsup1,shift= 1,dim=2),shift=-1,dim=1) + &
                cshift(cshift(bsup1,shift=-1,dim=2),shift= 1,dim=1) + &
                cshift(cshift(bsup1,shift= 1,dim=2),shift= 1,dim=1))/4.0

    !Clean up after shifts due to odd number of points
    Av_bsup1(:,0) = 0.0
    Av_bsup1(:,Num_u3-1) = 0.0
    Av_bsup1(0,:) = 0.0
    Av_bsup1(Num_u1-1,:) = 0.0

    !for b3 field aligned at North Ionosphere
    av_bsup1(:,0) = (cshift(bsup1(:,1),shift=-1) + cshift(bsup1(:,1),shift=1))/2.0
    Av_bsup1(0,0) = 0.0
    Av_bsup1(Num_u1-1,0) = 0.0

    !for b3 field aligned at Southern Ionosphere
    av_bsup1(:,num_u3-1) = (cshift(bsup1(:,num_u3-2),shift=-1) +&
                            cshift(bsup1(:,num_u3-2),shift= 1))/2.0
    Av_bsup1(0,Num_u3-1) = 0.0
    Av_bsup1(Num_u1-1,Num_u3-1) = 0.0

    av_bsup3 = (cshift(cshift(bsup3,shift=-1,dim=2),shift=-1,dim=1) + &
                cshift(cshift(bsup3,shift= 1,dim=2),shift=-1,dim=1) + &
                cshift(cshift(bsup3,shift=-1,dim=2),shift= 1,dim=1) + &
                cshift(cshift(bsup3,shift= 1,dim=2),shift= 1,dim=1))/4.0
    Av_bsup3(:,0) = 0.0
    Av_bsup3(:,Num_u3-1) = 0.0
    Av_bsup3(0,:) = 0.0
    Av_bsup3(Num_u1-1,:) = 0.0

    !Rotate to bi
    b1 = g11*bsup1 + g13*av_bsup3
    b2 = g22*bsup2
    b3 = g13*av_bsup1 + g33*bsup3

    !b2(0,:) 		= (4.0*b2(2,*) - b2(4,*))/3.0					          !Bob's fortran
    !b2(num_u1-1,:)	= (4.0*b2(num_u1-3,*) - b2(num_u1-5,*))/3.0 !Bob's fortran
    !b3(0,*) 		= (4.0*b3(2,*) - b3(4,*))/3.0					          !Bob's fortran

    ! Outer L Shell Boundary add driver boundary here
    b3(num_u1-1,:) = Driver*evens_u3

    b1p3 = cshift(b1,shift= 1,dim=1)
    b1m3 = cshift(b1,shift=-1,dim=2)
    b1p3(:,num_u3-1)= 0.0
    b1m3(:,0) = 0.0

    b2p1 = cshift(b2,shift= 1,dim=1)
    b2m1 = cshift(b2,shift=-1,dim=1)
    b2p1(num_u1-1,:)= 0.0
    b2m1(0,:) = 0.0

    b2p3 = cshift(b2,shift=1,dim=2)
    b2m3 = cshift(b2,shift=-1,dim=2)
    b2p3(:,num_u3-1)= 0.0
    b2m3(:,0) = 0.0

    b3p1 = cshift(b3,shift=1,dim=1)
    b3m1 = cshift(b3,shift=-1,dim=1)
    b3p1(num_u1-1,:)= 0.0
    b3m1(0,:) = 0.0
  end do
  write(*,*) 'Iterate Done!!!'

end subroutine iterate
