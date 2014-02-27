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
    e2m3,e3p1,e3m1,av_bsup1,av_bsup3,b2m3,b3p1,b3m1,b1p3,b1m3,b2p1,b2m1,b2p3
  complex,dimension(0:num_u1-1,0:num_u3-1)::  f1,av_f1,f2,av_f2,av_b3,av_esup1,av_esup2,av_e3,&
    bnu_M,bph_M,bmu_M,enu_M,eph_M,emu_M
  complex,allocatable,dimension(:,:) :: db2_d3,av_db2_d3
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
  !shift arrays
  integer,allocatable,dimension(:) :: isort,ipiv,evenx,oddx,evenz,oddz,evenxp,&
    evenxm,oddxp,oddxm,evenzp,evenzm,oddzp,oddzm
  integer,allocatable,dimension(:) :: evenxp3,evenxm3,oddxp3,oddxm3
  integer :: nx,nz,nxhalf,nzhalf

  allocate(db2_d3(0:num_u1-1,0:num_u3-1))
  allocate(av_db2_d3(0:num_u1-1,0:num_u3-1))
  nx=num_u1-1
  nz=num_u3-1
  nxhalf=nx/2
  nzhalf=nz/2

  !Zero Arrays
  e1=0.0
  e2=0.0
  e3=0.0
  b1=0.0
  b2=0.0
  b3=0.0
  f1=0.0
  f2=0.0
  esup1=0.0
  esup2=0.0
  bsup1=0.0
  bsup2=0.0
  bsup3=0.0
  cfac=0.0
  odds_u1=0.0
  evens_u1=0.0
  odds_u3=0.0
  evens_u3=0.0
  u3_in=0.0
  driver=0.0
  Eidx=0
  Oidx=0
  av_bsup1=0.0
  av_bsup3=0.0
  av_f1=0.0
  av_f2=0.0
  av_b3=0.0
  av_esup1=0.0
  av_esup2=0.0
  av_e3=0.0
  db2_d3=0.0
  av_db2_d3=0.0

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

  allocate(evenz(0:nzhalf),oddz(nzhalf),evenx(0:nxhalf),oddx(nxhalf),&
    oddxp(nxhalf),oddxm(nxhalf),oddzp(nzhalf),oddzm(nzhalf))
  allocate(evenxp(0:nxhalf),evenxm(0:nxhalf),evenzp(0:nzhalf),&
    evenzm(0:nzhalf),oddxp3(nxhalf),oddxm3(nxhalf),evenxm3(0:nxhalf),evenxp3(0:nxhalf))
  do ii=0,nzhalf
    evenz(ii)=2*ii
  enddo
  do ii=1,nzhalf
    oddz(ii)=2*ii-1
  enddo
  do ii=0,nxhalf
    evenx(ii)=2*ii
  enddo
  do ii=1,nxhalf
    oddx(ii)=2*ii-1
  enddo  
  oddxp=evenx(1:nxhalf)
  oddxm=evenx(0:nxhalf-1)
  oddzp=evenz(1:nzhalf)
  oddzm=evenz(0:nzhalf-1)
  evenxp(0:nxhalf-1)=oddx
  evenxp(nxhalf)=oddx(nxhalf)
  evenxm(1:nxhalf)=oddx
  evenxm(0)=oddx(1)
  evenzp(0:nzhalf-1)=oddz
  evenzp(nzhalf)=oddz(nzhalf)
  evenzm(1:nzhalf)=oddz
  evenzm(0)=oddz(1)
  oddxp3=min(oddx+3,nx)
  oddxm3=max(oddx-3,0)
  evenxp3=min(evenx+3,nx-1)
  evenxm3=max(evenx-3,1)

  write(*,*)'Starting Time Loop'
  do while (iplot .lt. nplots)
    tt = tt+dt
    call get_driver(tt,u3_in,driver)

    !shift ei's
!    e1p3 = cshift(e1,shift=1,dim=2)
    !e1m3 = cshift(e1,shift=-1,dim=2)
    !e1p3(:,num_u3-1) = 0.0
    !e1m3(:,0) = 0.0

    !e2p1 = cshift(e2,shift=1,dim=1)
    !e2m1 = cshift(e2,shift=-1,dim=1)
    !e2p1(num_u1-1,:) = 0.0
    !e2m1(0,:) = 0.0

    !e2p3 = cshift(e2,shift=1,dim=2)
    !e2m3 = cshift(e2,shift=-1,dim=2)
    !e2p3(:,num_u3-1) = 0.0
    !e2m3(:,0) = 0.0

    !e3p1 = cshift(e3,shift=1,dim=1)
    !e3m1 = cshift(e3,shift=-1,dim=1)
    !e3p1(num_u1-1,:) = 0.0
    !e3m1(0,:) = 0.0

    !Advance bsupi
    !bsup1 = bsup1 + b1e2*(e2p3-e2m3) + b1e3*e3
    !bsup2 = bsup2 + b2e1*(e1p3-e1m3) + b2e3*(e3p1-e3m1)
    !bsup3 = bsup3 + b3e1*e1 + b3e2*(e2p1-e2m1)

    !$OMP PARALLEL
      !$OMP DO PRIVATE(ii)
        do ii=1,nx-1,2
          bsup1(ii,oddz)=bsup1(ii,oddz)+b1e2(ii,oddz)*(e2(ii,oddzp)-&
            e2(ii,oddzm))+b1e3(ii,oddz)*e3(ii,oddz)
        enddo
      !$OMP END DO

      !$OMP DO PRIVATE(kk)
        do kk=1,nz-1,2
          bsup2(evenx,kk)=bsup2(evenx,kk)+b2e1(evenx,kk)*(e1(evenx,kk+1)-&
            e1(evenx,kk-1))+b2e3(evenx,kk)*(e3(evenxp,kk)-e3(evenxm,kk))
        enddo
      !$OMP END DO

      !$OMP DO PRIVATE(kk)
        do kk=0,nz,2
          bsup3(evenx,kk)=bsup3(evenx,kk)+b3e1(evenx,kk)*e1(evenx,kk)+&
            b3e2(evenx,kk)*(e2(evenxp,kk)-e2(evenxm,kk))
        enddo
      !$OMP END DO

      !$OMP SINGLE

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
        !$OMP END SINGLE

        !$OMP DO PRIVATE(kk)
          do kk=1,nz-1,2
            av_bsup3(oddx,kk)=0.25*(bsup3(oddxp,kk+1)+bsup3(oddxm,kk+1)+&
              bsup3(oddxp,kk-1)+bsup3(oddxm,kk-1))
          enddo
        !$OMP END DO
        !  bsup1(evenx,evenz)=0.25*(bsup1(evenxp,evenzp)+bsup1(evenxm,evenzp)+bsup1(evenxp,evenzm)+bsup1(evenxm,evenzm))
        !$OMP DO PRIVATE(kk)
          do kk=2,nz-2,2
            av_bsup1(evenx,kk)=0.25*(bsup1(evenxp,kk+1)+bsup1(evenxm,kk+1)+&
              bsup1(evenxp,kk-1)+bsup1(evenxm,kk-1))
          enddo
        !$OMP END DO

        !Clean up after shifts due to odd number of points
        !$OMP SINGLE
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

          Av_bsup3(:,0) = 0.0
          Av_bsup3(:,Num_u3-1) = 0.0
          Av_bsup3(0,:) = 0.0
          Av_bsup3(Num_u1-1,:) = 0.0
        !$OMP END SINGLE

        !Rotate to bi
        !$OMP DO PRIVATE(kk)
          do kk=1,nz-1,2
            b1(oddx,kk)=g11(oddx,kk)*bsup1(oddx,kk)+g13(oddx,kk)*av_bsup3(oddx,kk)
            b2(evenx,kk)=g22(evenx,kk)*bsup2(evenx,kk)
          enddo
        !$OMP END DO

        !$OMP DO PRIVATE(kk)
          do kk=0,nz,2
            b3(evenx,kk)=g13(evenx,kk)*av_bsup1(evenx,kk)+g33(evenx,kk)*bsup3(evenx,kk)
          enddo
        !$OMP END DO

        !$OMP SINGLE
            b3(num_u1-1,:) = Driver*evens_u3
            !b3(0,:) = 0.0
        !$OMP END SINGLE

        !$OMP DO PRIVATE(ii)
          do ii=0,nx,2
            f1(ii,evenz)=f1b2(ii,evenz)*(b2(ii,evenzp)-b2(ii,evenzm))+f1b3(ii,evenz)*b3(ii,evenz)
          enddo
        !$OMP END DO
        !   f2(oddx,evenz)=f2b1(oddx,evenz)*(b1(oddx,evenzp)-b1(oddx,evenzm))+f2b3(oddx,evenz)*(b3(oddxp,evenz)-b3(oddxm,evenz))
        !$OMP DO PRIVATE(kk)
          do kk=2,nz-2,2
            f2(oddx,kk)=f2b1(oddx,kk)*(b1(oddx,kk+1)-b1(oddx,kk-1))+&
              f2b3(oddx,kk)*(b3(oddxp,kk)-b3(oddxm,kk))
          enddo
        !$OMP END DO

        !$OMP SINGLE
          f1(:,0)= 0.0
          f1(:,Num_u3-1) = 0.0
          f2(:,0) = 0.0
          f2(:,Num_u3-1) = 0.0
        !$OMP END SINGLE

        !$OMP DO PRIVATE(kk)
          do kk=0,nz,2
            av_esup1(oddx,kk)=0.5*(esup1(oddxp,kk)+esup1(oddxm,kk))
            av_esup2(evenx,kk)=0.5*(esup2(evenxp,kk)+esup2(evenxm,kk))
            av_f1(oddx,kk)=0.5*(f1(oddxp,kk)+f1(oddxm,kk))
            av_f2(evenx,kk)=0.5*(f2(evenxp,kk)+f2(evenxm,kk))
          enddo
        !$OMP END DO

        !$OMP SINGLE
          Av_f1(0,:) = 0.0
          Av_f1(Num_u1-1,:) = 0.0
          Av_f2(0,:) = 0.0
          Av_f2(Num_u1-1,:) = 0.0
          av_esup1(0,:) = 0.0
          av_esup1(num_u1-1,:) = 0.0
          av_esup2(0,:) = 0.0
          av_esup2(num_u1-1,:) = esup2(num_u1-2,:)
        !$OMP END SINGLE

        !$OMP DO PRIVATE(kk)
          do kk=0,nz,2
            esup1(evenx,kk)=e1e1(evenx,kk)*esup1(evenx,kk)+e1e2(evenx,kk)*av_esup2(evenx,kk)+&
              e1f1(evenx,kk)*f1(evenx,kk)+e1f2(evenx,kk)*av_f2(evenx,kk)
          enddo
        !$OMP END DO

        !$OMP DO PRIVATE(kk)
          do kk=0,nz,2
            esup2(oddx,kk)=e2e1(oddx,kk)*av_esup1(oddx,kk)+e2e2(oddx,kk)*esup2(oddx,kk)+&
              e2f1(oddx,kk)*av_f1(oddx,kk)+e2f2(oddx,kk)*f2(oddx,kk)
          enddo
        !$OMP END DO

        !$OMP DO PRIVATE(kk)
          do kk=1,nz-1,2
            av_b3(oddx,kk)=0.25*(b3(oddxp,kk+1)+b3(oddxp,kk-1)+b3(oddxm,kk+1)+b3(oddxm,kk-1))
          enddo
        !$OMP END DO

        !$OMP SINGLE
          Av_b3(:,0) = 0.0
          Av_b3(:,Num_u3-1) = 0.0
          Av_b3(0,:) = 0.0
          Av_b3(Num_u1-1,:) = 0.0
        !$OMP END SINGLE

        !$OMP DO PRIVATE(kk)
          do kk = 0,nx,2
            db2_d3(kk,evenz) = b2(kk,evenzp)-b2(kk,evenzm)
          end do
        !$OMP END DO

        !db2_d3 = (cshift(b2,shift=1,dim=2)-cshift(b2,shift=-1,dim=2))


        !$OMP SINGLE
          !!Clean up after shifts (not on any boundaries)
          db2_d3(:,0) = 0.0
          db2_d3(:,Num_u3-1) = 0.0
          db2_d3(0,:) = 0.0
          db2_d3(Num_u1-1,:) = 0.0
        !$OMP END SINGLE

        !$OMP DO PRIVATE(kk)
          do kk=2,nz-2,2
            av_db2_d3(oddx,kk)=0.25*(db2_d3(oddxp,kk+1)+db2_d3(oddxp,kk-1)+&
              db2_d3(oddxm,kk+1)+db2_d3(oddxm,kk-1))
          enddo
        !$OMP END DO

        !$OMP SINGLE
          av_db2_d3(:,1) =        (cshift(db2_d3(:,2),shift=1) + &
                                  cshift(db2_d3(:,2),shift=-1))/2.0

          Av_db2_d3(:,Num_u3-2) = (cshift(db2_d3(:,Num_u3-3),shift=1) + &
                                  cshift(db2_d3(:,Num_u3-3),shift=-1))/2.0

          !!Clean up after shifts (not on any boundaries)
          Av_db2_d3(:,0) = 0.0
          Av_db2_d3(:,Num_u3-1) = 0.0
          Av_db2_d3(0,:) = 0.0
          Av_db2_d3(Num_u1-1,:) = 0.0
        !$OMP END SINGLE

        !advance eparallel
        !e3 = e3b21*(b2p1-b2m1) + e3b23*(Av_db2_d3) + e3b1*b1  + e3b3*Av_b3

        !$OMP DO PRIVATE(kk)
          do kk=1,nz-1,2
            e3(oddx,kk)=e3b1(oddx,kk)*b1(oddx,kk)+e3b21(oddx,kk)*(b2(oddxp,kk)-b2(oddxm,kk))+e3b23(oddx,kk)*av_db2_d3(oddx,kk)+e3b3(oddx,kk)*av_b3(oddx,kk)
          enddo
        !$OMP END DO

        !$OMP SINGLE
          e3(0,:) = 0.0
          e3(num_u1-1,:) = 0.0
          e3(:,0) = 0.0
          e3(:,Num_u3-1) = 0.0
        !$OMP END SINGLE

        !$OMP DO PRIVATE(kk)
          do kk=2,nz-2,2
            av_e3(evenx,kk)=0.25*(e3(evenxp,kk+1)+e3(evenxp,kk-1)+e3(evenxm,kk+1)+e3(evenxm,kk-1))
          enddo
        !$OMP END DO

        !$OMP SINGLE
          av_e3(:,0) = 0.0
          Av_e3(:,Num_u3-1) = 0.0
          Av_e3(0,:) = 0.0
          Av_e3(Num_u1-1,:) = 0.0
        !$OMP END SINGLE

        !$OMP DO PRIVATE(kk)
          do kk=0,nz,2
            e1(evenx,kk)=e1esup1(evenx,kk)*esup1(evenx,kk)+e1e3(evenx,kk)*av_e3(evenx,kk)
            e2(oddx,kk)=g22(oddx,kk)*esup2(oddx,kk)
          enddo
        !$OMP END DO

    !$OMP END PARALLEL

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
