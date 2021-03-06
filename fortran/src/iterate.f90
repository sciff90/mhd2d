subroutine iterate()
  use global
  implicit none

  !Initialize Arrays
  complex,dimension(0:num_u1-1,0:num_u3-1) :: e1,e2,e3,b1,b2,b3
  complex,dimension(0:num_u1-1,0:num_u3-1) :: esup1,esup2,bsup1,bsup2,bsup3
  real,dimension(0:num_u1-1,0:num_u3-1) :: cfac
  complex,dimension(0:num_u1-1) :: odds_u1,evens_u1
  complex,dimension(0:num_u3-1) :: odds_u3,evens_u3
  real,dimension(0:num_u3-1) :: u3_in,driver
  integer,dimension(0:num_u1/2) :: Eidx
  integer,dimension(0:num_u1/2-1) :: Oidx

  complex,dimension(0:num_u1-1,0:num_u3-1) :: e1p3,e1m3,e2p1,e2m1,e2p3,&
    e2m3,e3p1,e3m1,av_bsup1,av_bsup3,b2m3,b3p1,b3m1,b1p3,b1m3,b2p1,b2m1,b2p3,&
    f1,av_f1,f2,av_f2,av_b3,db2_d3,av_db2_d3,av_esup1,av_esup2,av_e3,&
    bnu_M,bph_M,bmu_M,enu_M,eph_M,emu_M
  complex,dimension(0:num_u1/2) :: bsup3_n,bsup3_s
  complex,dimension(0:num_u1-1) :: psiatm_N,psiatmp,psiatmm,b1atm_N,&
    b2atm_N,b1_N,b2_N,b1dif_N,b2dif_N,b2atm_S,b1_S,b2_S,b1dif_S,b2dif_S,&
    psiatm_S,b1atm_s,psig,psigp,psigm,bxg_N,bxg_S,byg_N,byg_S


  !Counters
  integer*4 :: ii,jj,kk,tpnext
  integer*4 :: iplot,nplots
  real :: tt
  integer :: reclength

  !Strings
  character(len=100) :: filename


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
  u3_in(ii) = dble(ii)
  end do

  cfac = dt*sqrt(va2_arr)*0.5/h3
  nplots = tmax/dtplot
  tpnext = dtplot
  iplot = 0
  tt=0.0

  write(*,*)'Starting Time Loop'
  do while (iplot .lt. nplots)
    tt = tt+dt
    call get_driver(tt,u3_in,driver)

    !shift ei's
    e1p3 = cshift(e1,shift=1,dim=2)
    e1m3 = cshift(e1,shift=-1,dim=2)
    e1p3(:,num_u3-1) = 0.0
    e1m3(:,0) = 0.0

    e2p1 = cshift(e2,shift=1,dim=1)
    e2m1 = cshift(e2,shift=-1,dim=1)
    e2p1(num_u1-1,:) = 0.0
    e2m1(0,:) = 0.0

    e2p3 = cshift(e2,shift=1,dim=2)
    e2m3 = cshift(e2,shift=-1,dim=2)
    e2p3(:,num_u3-1) = 0.0
    e2m3(:,0) = 0.0

    e3p1 = cshift(e3,shift=1,dim=1)
    e3m1 = cshift(e3,shift=-1,dim=1)
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
      !bsup3(num_u1-1,:) = 0.d0
      !bsup3(num_u1-1,0) = 0.d0
      !bsup3(num_u1-1,num_u3-1])= 0.d0
    !else
      !bsup2(num_u1-1,:) = (4.0*bsup2(num_u1-3,*) - bsup2(num_u1-5,*))/3.0
      !bsup3(num_u1-1,:) = (4.0*bsup3(num_u1-3,*) - bsup3(num_u1-5,*))/3.0
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

    b1p3 = cshift(b1,shift= 1,dim=2)
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

    !at e1 location used in e1 equation
    f1 = f1b2*(b2p3-b2m3) + f1b3*b3
    !clean up at ends of field lines
    f1(:,0)= 0.0
    f1(:,Num_u3-1) = 0.0

    !at e2 location used in e2 equation
    Av_f1 = (cshift(f1,shift=-1,dim=1) + cshift(f1,shift=1,dim=1))/2
    Av_f1(0,:) = 0.0
    Av_f1(Num_u1-1,:) = 0.0

    !at e2 location used in e2 equation
    f2 = f2b1*(b1p3-b1m3) + f2b3*(b3p1-b3m1)
    f2(:,0) = 0.0
    f2(:,Num_u3-1) = 0.0

    Av_f2 =(cshift(f2,shift=-1,dim=1) + cshift(f2,shift=1,dim=1))/2
    Av_f2(0,:) = 0.0
    Av_f2(Num_u1-1,:) = 0.0

    !Need to shift e's to all grid lvls with e1 and e2
    !need e1 at e2 locations (e2 is not on boundaries)
    Av_esup1 = (cshift(esup1,shift=-1,dim=1) + cshift(esup1,shift=1,dim=1))/2
    !Inner L shells is not required here (cleaning up from shifts)
    Av_esup1(0,:) = 0.0
    !Outer L shells is not required here (cleaning up from shifts)
    Av_esup1(Num_u1-1,:)= 0.0
    !need e2 at e1 locations
    Av_esup2 = (cshift(esup2,shift=-1,dim=1) + cshift(esup2,shift=1,dim=1))/2
    !Inner L shells PEC BC so Esup2 is zero
    Av_esup2(0,:) = 0.0
    !Outer L shells (equivalent to De2du1 = 0)
    Av_esup2(Num_u1-1,:)= esup2(Num_u1-2,:)

    !advance eperp
    esup1 = (e1e1*   esup1 + e1e2*Av_esup2 + e1f1*   f1 + e1f2*Av_f2)
    esup2 = (e2e1*Av_esup1 + e2e2*   esup2 + e2f1*Av_f1 + e2f2*   f2)


    !esup1[0,*] 			= (4.0*esup1(2,*) - esup1(4,*))/3.0					!Bob's fortran
    !esup1[num_u1-1,*] 	= (4.0*esup1(num_u1-3,*) - esup1(num_u1-5,*))/3.0	!Bob's fortran

    !need b3 at e3 (and b1) locations (not on any boundaries)
    av_b3 = (cshift(cshift(b3,shift=-1,dim=2),shift=-1,dim=1) + &
             cshift(cshift(b3,shift= 1,dim=2),shift=-1,dim=1) + &
             cshift(cshift(b3,shift=-1,dim=2),shift= 1,dim=1) + &
             cshift(cshift(b3,shift= 1,dim=2),shift= 1,dim=1))/4.0
    Av_b3(:,0) = 0.0
    Av_b3(:,Num_u3-1) = 0.0
    Av_b3(0,:) = 0.0
    Av_b3(Num_u1-1,:) = 0.0

    !Need to put the db2du3 derivatives on the right grid point
    !for a yee grid and the fact we have a non orthognal system
    db2_d3 = (b2p3-b2m3)

    !Clean up after shifts (not on any boundaries)
    db2_d3(:,0) = 0.0
    db2_d3(:,Num_u3-1) = 0.0
    db2_d3(0,:) = 0.0
    db2_d3(Num_u1-1,:) = 0.0

    av_db2_d3 = (cshift(cshift(db2_d3,shift=-1,dim=2),shift=-1,dim=1) + &
                 cshift(cshift(db2_d3,shift= 1,dim=2),shift=-1,dim=1) + &
                 cshift(cshift(db2_d3,shift=-1,dim=2),shift= 1,dim=1) + &
                 cshift(cshift(db2_d3,shift= 1,dim=2),shift= 1,dim=1))/4.0

    Av_db2_d3(:,1) =        (cshift(db2_d3(:,2),shift=1) + &
                             cshift(db2_d3(:,2),shift=-1))/2.0

    Av_db2_d3(:,Num_u3-2) = (cshift(db2_d3(:,Num_u3-3),shift=1) + &
                             cshift(db2_d3(:,Num_u3-3),shift=-1))/2.0

    !Clean up after shifts (not on any boundaries)
    Av_db2_d3(:,0) = 0.0
    Av_db2_d3(:,Num_u3-1) = 0.0
    Av_db2_d3(0,:) = 0.0
    Av_db2_d3(Num_u1-1,:) = 0.0

    !advance eparallel
    e3 = e3b21*(b2p1-b2m1) + e3b23*(Av_db2_d3) + e3b1*b1  + e3b3*Av_b3
    !Clean up after shifts (not on any boundaries)
    e3(0,:) = 0.0
    e3(num_u1-1,:) = 0.0
    e3(:,0) = 0.0
    e3(:,Num_u3-1) = 0.0

    !killing eparallel for testing
    !e3[*,*] =0.0d0
    !average ei's
    !Av_esup1 	= (Shift(esup1,1,1) + Shift(esup1,1,-1) + Shift(esup1,-1,1) + Shift(esup1,-1,-1))/ 4.d0
    !Av_esup1(*,0)	= 0.0 		& 	Av_esup1(*,Num_u3-1)	= 0.0			; Clean up after shifts due to odd number of points
    !Av_esup1(0,*)	= 0.0 		& 	Av_esup1(Num_u1-1,*)	= 0.0

    !may need to do ionospheric sheet!
    av_e3 = (cshift(cshift(e3,shift=-1,dim=2),shift=-1,dim=1) + &
             cshift(cshift(e3,shift= 1,dim=2),shift=-1,dim=1) + &
             cshift(cshift(e3,shift=-1,dim=2),shift= 1,dim=1) + &
             cshift(cshift(e3,shift= 1,dim=2),shift= 1,dim=1))/4.0

    Av_e3(:,0) = 0.0
    Av_e3(:,Num_u3-1) = 0.0
    Av_e3(0,:) = 0.0
    Av_e3(Num_u1-1,:) = 0.0

    !rotate to ei
    e1 = e1esup1*esup1 + e1e3*Av_e3
    e2 = g22*esup2

    e1(0,0)=(4*e1(2,0)-e1(4,0))/3.0   !Bob's fortran
    e1(num_u1-1,0)=(4*e1(num_u1-3,0)-e1(num_u1-5,0))/3.0

    esup1(0,0) = (4*esup1(2,0)-esup1(4,0))/3.0  !Bob's fortran
    esup1(num_u1-1,0) = (4*esup1(num_u1-3,0)-esup1(num_u1-5,0))/3.0

    !Thin Sheet Ionospheric Boundaries
    !Northern Ionospheric Sheet
    !Average b3 at thin sheet
    !bsup3_n = (bsup3(:,0) + (cshift(bsup3(:,0),shift=-1)+ &
    !                         cshift(bsup3(:,0),shift=1))/2.0)

    !call Extrap_1D(bsup3_n,u1,0)
    !bsup3_n(0) = bsup3(0,0)
    !bsup3_n(Num_u1-1) = bsup3(Num_u1-1,0)

    !Only using B3 points to do fit
    bsup3_n = bsup3(Eidx,0)
    !Nth hemis
    psiatm_N = matmul(psiib3_N,bsup3_n)

    psiatmp = cshift(psiatm_N,shift=1)
    psiatmm = cshift(psiatm_N,shift=-1)
    b1atm_N = (psiatmp-psiatmm)/(dx2)
    call Extrap_1D(b1atm_N,u1,0)
    call Extrap_1D(b1atm_N,u1,Num_u1-1)
    !b1atm_N(0) = (-3.0*psiatm_N(0) + 4.0*psiatm_N(1) - psiatm_N(2))/(dx2)
    !b1atm_N(Num_u1-1) = ( 3.0*psiatm_N(Num_u1-1) - 4.0*psiatm_N(Num_u1-2) &
    !                      + psiatm_N(Num_u1-3))/(dx2)
    b2atm_N = im*psiatm_N

    b1_N = b1(:,1) +(cshift(b1(:,1),shift=-1)+cshift(b1(:,1),shift=1))/2
    call Extrap_1D(b1_N,u1,0)
    call Extrap_1D(b1_N,u1,Num_u1-1)
    !Bob's extrapolation
    !b1_N(0)			= 3.0*b1(1,1) - 3.5*b1(3,1) + 1.5*b1(5,1)
    !b1_N(Num_u1-1)	= 3.0*b1(Num_u1-2,1) - 3.5*b1(Num_u1-4,1) + 1.5*b1(Num_u1-6,1)

    b2_N = b2(:,1) +(cshift(b2(:,1),shift=-1)+cshift(b2(:,1),shift=1))/2
    !call Extrap_1D(b2_N,u1,0)
    b2_N(0) = b2(0,1)
    b2_N(Num_u1-1) = b2(Num_u1-1,1)
    b1dif_N = b1_N-b1atm_N
    b2dif_N = b2_N-b2atm_N

    e1(:,0) = (e1b1atm_N*b1dif_N + e1b2atm_N*b2dif_N)*evens_u1
    e2(:,0) = (e2b1atm_N*b1dif_N + e2b2atm_N*b2dif_N)*odds_u1

    !Southern Ionospheric Sheet
    !bsup3_s 	= (bsup3(*,Num_u3-1) + (cshift(bsup3a:,Num_u3-1,-1)
    !                                +cshift(bsup3(:,Num_u3-1),1))/2.0)
    !call Extrap_1D(bsup3_s,u1,0)
    !bsup3_s(0)	= bsup3(0,Num_u3-1)
    !bsup3_s(Num_u1-1)	= bsup3(Num_u1-1,Num_u3-1)

    !Only using B3 points to do fit
    bsup3_s = bsup3(Eidx,Num_u3-1)
    psiatm_S = matmul(psiib3_s,bsup3_s)

    psiatmp = cshift(psiatm_S,1)
    psiatmm = cshift(psiatm_S,-1)

    b1atm_S = (psiatmp-psiatmm)/(dx2)
    call Extrap_1D(b1atm_S,u1,0)
    call Extrap_1D(b1atm_S,u1,Num_u1-1)
    !b1atm_S(0) 		  	= (-3.0*psiatm_S(0) + 4.0*psiatm_S(1)&
    !                     - psiatm_S(2))/(dx2)
    !b1atm_S(Num_u1-1) 	= ( 3.0*psiatm_S(Num_u1-1) - 4.0*psiatm_S(Num_u1-2)&
    !                       + psiatm_S(Num_u1-3))/(dx2)
    b2atm_S = im*psiatm_S


    b1_S = b1(:,Num_u3-2) +(cshift(b1(:,Num_u3-2),-1)+&
                            cshift(b1(:,Num_u3-2),1))/2
    call Extrap_1D(b1_S,u1,0)
    call Extrap_1D(b1_S,u1,Num_u1-1)
    !Bob's extrapolation
    !b1_S(0)			= 3.0*b1(1,Num_u3-2) - 3.5*b1(3,Num_u3-2)
    !                   + 1.5*b1(5,Num_u3-2)
    !b1_S(Num_u1-1)	= 3.0*b1(Num_u1-2,Num_u3-2) - 3.5*b1(Num_u1-4,Num_u3-2)&
    !                   + 1.5*b1(Num_u1-6,Num_u3-2)

    b2_S = b2(:,Num_u3-2) +(cshift(b2(:,Num_u3-2),-1)+&
                            cshift(b2(:,Num_u3-2),1))/2
    !call Extrap_1D(b2_S,u1,0)
    b2_S(0) = b2(0,Num_u3-2)
    b2_S(Num_u1-1) = b2(Num_u1-1,Num_u3-2)

    b1dif_S = b1_S - b1atm_S
    b2dif_S = b2_S - b2atm_S

    e1(:,num_u3-1) = (e1b1atm_S*b1dif_S + e1b2atm_S*b2dif_S)*evens_u1
    e2(:,num_u3-1) = (e2b1atm_S*b1dif_S + e2b2atm_S*b2dif_S)*odds_u1

    !temporary ionospheric bc - for testing
    !perfectly reflecting Ionosphere (thin sheet)
    !e1(:,0)=0.0d0
    !e2(:,0)=0.0d0

    !e1(:,Num_u3-1)=0.0d0
    !e2(:,Num_u3-1)=0.0d0


    !This rotation folds the field align E
    !into the Ionospheric boundary (thin sheet) - MDS

    !Av_e3_N = (cshift(e3(:,1),1) + cshift(e3(:,1),-1))/2.0
    !esup1(:,0) = gsup11(:,0)*e1(:,0)+gsup13(:,0)*e3(:,1)

    !esup1(:,0) = gsup11(:,0)*e1(:,0)+gsup13(:,0)*Av_e3_N
    !esup2(:,0) = gsup22(:,0)*e2(:,0)

    !Av_e3_S = (csift(e3(:,num_u3-2),-1) + cshift(e3(:,num_u3-2),1))/2.0
    !esup1(:,num_u3-1) = gsup11(:,num_u3-1)*e1(:,num_u3-1)+&
      !gsup13(:,num_u3-1)*e3(:,num_u3-2)

    !esup1(:,num_u3-1) = gsup11(:,num_u3-1)*e1(:,num_u3-1)+&
      !gsup13(:,num_u3-1)*Av_e3_S
    !esup2(:,num_u3-1) = gsup22(:,num_u3-1)*e2(:,num_u3-1)


    !BC on e1's (Bob's Email) st de1du1 = 0
    !e1(0,0) = (4.0*e1(2,0)-e1(4,0))/3.0
    !e1(num_u1-1,0) = (4.0*e1(num_u1-3,0)-e1(num_u1-5,0))/3.0

    !esup1(0,0) = (4.0*esup1(2,0)-esup1(4,0))/3.0
    !esup1(num_u1-1,0) =(4.0*esup1(num_u1-3,0)-esup1(num_u1-5,0))/3.0

    !e1(0,Num_u3-1) = (4.0*e1(2,Num_u3-1)-e1(4,Num_u3-1))/3.0
    !e1(num_u1-1,Num_u3-1) =(4.0*e1(num_u1-3,Num_u3-1)-&
    !                          e1(num_u1-5,Num_u3-1))/3.0

    !esup1(0,Num_u3-1) =(4.0*esup1(2,Num_u3-1)-esup1(4,Num_u3-1))/3.0
    !esup1(num_u1-1,Num_u3-1) =(4.0*esup1(num_u1-3,Num_u3-1)&
    !                            -esup1(num_u1-5,Num_u3-1))/3.0

    !write output and plot
    bnu_M = 0.0
    bph_M = 0.0
    bmu_M = 0.0
    enu_M = 0.0
    eph_M = 0.0
    emu_M = 0.0

    if(tt .ge. tpnext) then
      iplot=iplot+1
      bnu_M = 1.0e9*bsup1*h1/Re        !B in nT:  mu, nu and phi components
      bph_M = 1.0e9*bsup2*h2/Re
      bmu_M = 1.0e9*b3/h3/Re
      enu_M = 1.0e3*esup1*h1           !E in mV/m
      eph_M = 1.0e3*esup2*h2
      emu_M = 1.0e3*e3/h3

      !b1_Nsc = 1.0e9*b1_N/hthatm_N/Re
      !b2_Nsc = 1.0e9*b2_N/hphiatm_N/Re
      !b1atm_Nsc = 1.0e9*b1atm_N/hthatm_N/Re
      !b2atm_Nsc = 1.0e9*b2atm_N/hphiatm_N/Re

      !b1_Ssc = 1.0e9*b1_S/hthatm_S/Re
      !b2_Ssc = 1.0e9*b2_S/hphiatm_S/Re
      !b1atm_Ssc = 1.0e9*b1atm_S/hthatm_S/Re
      !b2atm_Ssc = 1.0e9*b2atm_S/hphiatm_S/Re

      !Calculate ground signatures
      !Nth Hemisphere
      psig = matmul(psigb3_N,bsup3_N)
      psigp = cshift(psig,1)
      psigm = cshift(psig,-1)
      psigp(num_u1-1) = psig(num_u1-3) + 3.d0*(psig(num_u1-1)-psig(num_u1-2))
      psigm(0) = psig(2) + 3.d0*(psig(0)-psig(1))
      bxg_N = -(psigp-psigm)/dx2/hthg_N
      byg_N = im*psig/hphig_N

      !Sth Hemisphere
      psig = matmul(psigb3_S,bsup3_S)
      psigp = cshift(psig,1)
      psigm = cshift(psig,-1)
      psigp(num_u1-1) = psig(num_u1-3) + 3.d0*(psig(num_u1-1)-psig(num_u1-2))
      psigm(0) = psig(2) + 3.d0*(psig(0)-psig(1))
      bxg_S = (psigp-psigm)/dx2/hthg_S
      byg_S = im*psig/hphig_S

      if(iplot.eq.1)then
        write(filename,*),'./data/plotting_constants.dat'

        inquire(iolength = reclength) num_u1,num_u3,&
                                      va2_arr,&
                                      x_arr,y_arr,&
                                      Eidx,dx2

        open(unit=10,file=filename,form='unformatted',access='direct',&
        recl = reclength)

        write(10,rec=1)&
          num_u1,num_u3,&
          sqrt(Va2_arr),&
          x_arr,y_arr,&
          Eidx,dx2
        close(unit=10)

      end if

      !Output from iterate
      !inquire(iolength = reclength)
      write(filename,'(a,i4.4,a4)'),'./data/plot_data',iplot,'.dat'
      inquire(iolength=reclength) tt,&
       Enu_M, Eph_M, Emu_M, Bnu_M, Bph_M, Bmu_M, &
       b1_N, b1atm_N, b2_N, b2atm_N, &
       b1_S, b1atm_S, b2_S, b2atm_S, &
       psiatm_N, psiatm_S, &
       im, &
       bsup3_n,  bsup3_s, &
       hratm_N, hratm_S, &
       bxg_N, byg_N, bxg_S, byg_S

     open(unit=10, file=filename,&
                    form='unformatted',&
                    access='direct',&
                    recl= reclength)

      write(10,rec=1)&
       tt,&
       Enu_M, Eph_M, Emu_M, Bnu_M, Bph_M, Bmu_M, &
       b1_N, b1atm_N, b2_N, b2atm_N, &
       b1_S, b1atm_S, b2_S, b2atm_S, &
       psiatm_N, psiatm_S, &
       im, &
       bsup3_n,  bsup3_s, &
       hratm_N, hratm_S, &
       bxg_N, byg_N, bxg_S, byg_S
      close(unit=10)
      write(*,*) 'File Written to:',filename

      tpnext = tpnext+dtplot

      !Test for crash
      if((maxval(abs(emu_M)).gt.1.0e12) .or. (maxval(abs(eph_M)).gt. 1.0e12)) then
        !CRASHED BREAKOUT OF LOOP
      end if
    end if


  end do
  write(*,*) 'Iterate Done!!!'

end subroutine iterate

subroutine extrap_1D(array,grid,location)
  use global
  implicit none
  integer :: location
  real,dimension(0:num_u1-1),intent(in) :: grid
  complex,dimension(0:num_u1-1) :: array

  integer :: min_idx_array,max_idx_array,min_idx_grid,max_idx_grid

  real :: hh,hh1
  complex :: fir_d,sec_d
  min_idx_array = 0
  max_idx_array = num_u1-1
  min_idx_grid = 0
  max_idx_grid = num_u1-1

  if( location .eq. 0 )then   !Forwards Extrapolation
    hh = grid(min_idx_grid+1)-grid(min_idx_grid+2)

    fir_d = (array(min_idx_array+3)-4.0*array(min_idx_array+2)&
      + 3.0*array(min_idx_array+1))/(2.0*hh)

    sec_d = (-array(min_idx_array+4) + 4.0*array(min_idx_array+3)&
      - 5.0*array(min_idx_array+2) + 2.0*array(min_idx_array+1))/(hh**2)

    hh1 = grid(min_idx_grid)-grid(min_idx_grid+1)

    array(min_idx_array) = array(min_idx_array+1)&
      +hh1*fir_d+((hh1**2)/2.0)*sec_d
  else
    !step size along ionosphere is uniform on B1_S grid
    hh = grid(max_idx_grid-1)- grid(max_idx_grid-2)
    !First Derivative Backwards difference O(h^2)
    fir_d = (array(max_idx_array-3) - 4.0*array(max_idx_array-2) +&
      3.0*array(max_idx_array-1))/(2.0*hh)
    !Second Derivative Backwards difference O(h^2)
    sec_d = (-array(max_idx_array-4) + 4.0*array(max_idx_array-3) -&
      5.0*array(max_idx_array-2) + 2.0*array(max_idx_array-1))/(hh**2)
    !Step size on B1_S_interp grid
    hh1 = grid(max_idx_grid) - grid(max_idx_grid-1)
    array(max_idx_array) = array(max_idx_array-1) + hh1*fir_d &
      + ((hh1**2)/2.0)*sec_d
  end if
end subroutine extrap_1D
