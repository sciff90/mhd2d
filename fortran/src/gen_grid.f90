subroutine gen_grid()
    use global
    use vector_operations
    implicit none

    !parameters
    double precision,parameter :: k0 = 8.02e15      !Magnetic moment in Wb.m
    double precision,parameter :: cosColat0 = sqrt(1.0-(r_iono/LVal))
    double precision,parameter :: r0 = r_iono*Re
    double precision,parameter :: Lat0 = asin(cosColat0)
    double precision,parameter :: ColLat0 = acos(cosColat0)
    double precision,parameter :: DelS = -0.1e3 !approx 100m spacing. Step along field line (- going outwards from northern hemisphere)


    !Variables
    double precision :: X,Z,r,Lat,Bx,Bz,Bx_unit,Bz_unit,d3
    double precision:: Sn,dsi,ccosth0,numin,numax
    integer :: ii,jj,temp
    integer*4 :: Npts
    !Arrays
    double precision,dimension(:),allocatable :: XX,ZZ,CLat,u3,ds,s,ss,y,cosCLat,cr,ccosth
    double precision,dimension(:),allocatable :: mu,rinit,dmudx0,dels_arr
    double precision,dimension(:,:),allocatable :: sinth2,sinth02,dmudx,costh02
    double precision,dimension(:,:),allocatable :: costh0,sinth0,costh2,sinth
    double precision,dimension(:,:),allocatable :: bfac
    double precision,dimension(:,:),allocatable :: gsup22,gsup33
    double precision,dimension(:,:),allocatable :: h30

    !start grid_gen
    X = r0*cos(Lat0)
    Z = r0*sin(Lat0)
    Npts = 0

    !Calculate cononical field line length using analytic dipole field
    do while (Z .gt. 0)
        r = sqrt(X**2+Z**2)
        Lat = atan(Z/X)
        Bx = -3.0*k0/r**3*sin(Lat)*cos(Lat)
        Bz = -2.0*k0/r**3*sin(Lat)**2+k0/r**3*cos(Lat)**2

        Bx_unit = Bx/sqrt(Bx**2+Bz**2)
        Bz_unit = Bz/sqrt(Bx**2+Bz**2)

        X = X + Bx_unit*DelS
        Z = Z + Bz_unit*DelS
        Npts = Npts+1
    end do

    Allocate(XX(0:Npts-1))
    Allocate(ZZ(0:Npts-1))
    Allocate(CLat(0:Npts-1))

    X = r0*cos(Lat0)
    Z = r0*sin(Lat0)
    Npts = 0

    do while (z .gt. 0)
            r = sqrt(X**2+Z**2)
            Lat = atan(Z/X)
            XX(Npts) = X
            ZZ(Npts) = Z
            CLat(Npts) = pi/2 -Lat

            Bx = -3.0*k0/r**3*sin(Lat)*cos(Lat)
            Bz = -2.0*k0/r**3*sin(Lat)**2+k0/r**3*cos(Lat)**2

            Bx_unit = Bx/sqrt(Bx**2+Bz**2)
            Bz_unit = Bz/sqrt(Bx**2+Bz**2)

            X = X + Bx_unit*DelS
            Z = Z + Bz_unit*DelS
            Npts = Npts+1
    end do
    Sn = real(Npts-1)*abs(DelS)
    num_u3_half = nint(2.0*Sn/(dsn))
    num_u3 = 2*num_u3_half + 1  !number of points must be odd

    Allocate(u3(0:num_u3-1))
    do ii = 0, num_u3-1         !u3 ordiniate is just index of array
            u3(ii) = ii
    end do


    dsi = 2.0*Sn/(num_u3_half*(num_u3_half+1)) !space coming from a sum of the arithmetic progression Sn = N/2*(2ds0+(N-1)dsi) and ds0=dsi

    Allocate(ds(0:num_u3_half-1))
    Allocate(s(0:num_u3_half))
    Allocate(ss(0:Npts-1))
    Allocate(cosCLat(0:Npts-1))

    do ii = 0, num_u3_half-1
            ds(ii) = dsi+ii*dsi
    end do

    do ii = 0, num_u3_half
            if (ii.eq.0) then
                    s(ii) = 0
            else
                    s(ii) = sum(ds(0:ii-1))
            end if
    end do

    do ii = 0, Npts-1
            ss(ii) = ii*abs(Dels)
    end do

    !Generate points along Cononical (LVal) field line
    do ii = 0, Npts-1
            cosCLat(ii) = cos(Clat(ii))
    end do

    Allocate(y(0:num_u3-1))
    y(0) = cosCLat(0)
    do ii = 1, num_u3_half-1

            do jj = 0,Npts-1
                if (ss(jj).lt.s(ii)) then
                       temp = jj
                end if
            end do

            y(ii) = cosCLat(temp)+(cosCLat(temp)-cosCLat(temp-1))*(s(ii)-ss(temp))/(ss(temp)-ss(temp-1))
            y(num_u3-1-ii) = -y(ii) !fill other hemisphere

    end do

    y(num_u3-1) = -y(0)
    y(num_u3_half) = 0.0        !equator point

    !c prefix for canonical form

    Allocate(ccosth(0:num_u3-1))
    Allocate(cr(0:num_u3-1))
    Allocate(mu(0:num_u3-1))
    ccosth0 = y(0)

    do ii = 0, num_u3-1
            ccosth(ii) = y(ii)
            cr(ii) = LVal*(1.0-y(ii)**2)        !cononical field line radial distance from center of earth (in Re)
            mu(ii) = r_iono**2*ccosth(ii)/(cr(ii)**2*ccosth0)
    end do

    write(*,*) "Number of cells = ",num_u3," rtop = ", maxval(cr)

    !Setup nu grid
    numin = -r_iono/LMin
    numax = -r_iono/LMax
    del_u1 = (numax-numin)/(num_u1-1)
    Allocate(u1(0:num_u1-1))
    do ii = 0, num_u1-1
            u1(ii) = ii*del_u1+numin
    end do

    Allocate(r_arr(0:num_u1-1,0:num_u3-1))
    Allocate(sinth2(0:num_u1-1,0:num_u3-1))
    Allocate(costh2(0:num_u1-1,0:num_u3-1))
    Allocate(sinth02(0:num_u1-1,0:num_u3-1))
    Allocate(sinth0(0:num_u1-1,0:num_u3-1))
    Allocate(costh(0:num_u1-1,0:num_u3-1))
    Allocate(sinth(0:num_u1-1,0:Num_u3-1))
    Allocate(dmudx(0:num_u1-1,0:num_u3-1))
    Allocate(costh02(0:num_u1-1,0:num_u3-1))
    Allocate(costh0(0:num_u1-1,0:num_u3-1))

    do ii = 0, num_u1-1
           sinth02(ii,:) = -u1(ii)
    end do

    costh02(:,:) = 1.0 - sinth02(:,:)   !Ionospheric Colat of field lines squared
    costh0(:,:) = sqrt(costh02(:,:))
    sinth0(:,:) = sqrt(sinth02(:,:))    !odd fnct

    !Northern Hemisphere
    Allocate(rinit(0:num_u3_half-1))

    do ii = 0, num_u1-1
            rinit = r_iono/sqrt(abs(mu(0:num_u3_half-1)*costh0(ii,0)))
            call New_r(rinit,u1(ii),mu(0:num_u3_half-1),costh02(ii,0),r_iono,r_arr(ii,0:num_u3_half-1),num_u3_half)
            sinth2(ii,0:num_u3_half-1) = -u1(ii)*r_arr(ii,0:num_u3_half-1)/r_iono
    end do

    !Extend to southern hemisphere
    do ii = 0, num_u3_half-1
        r_arr(:,num_u3-ii-1) = r_arr(:,ii)
        sinth2(:,num_u3-ii-1) = sinth2(:,ii)
    end do


    r_arr(:,num_u3_half) = -r_iono/u1(:)
    sinth2(:,num_u3_half) = 1.0
    costh2(:,:) = 1.0-sinth2(:,:)
    costh(:,0:num_u3_half-1) = sqrt(costh2(:,0:num_u3_half-1))
    !Northern Hemisphere
    costh(:,num_u3_half:num_u3-1) = -1.0*sqrt(costh2(:,num_u3_half:num_u3-1))
    !Southern Hemisphere
    sinth(:,:) = sqrt(sinth2(:,:))  !theta only runs 0<x<pi so sin(theta) > 0

    Allocate(x_arr(0:num_u1-1,0:num_u3-1))
    Allocate(y_arr(0:num_u1-1,0:num_u3-1))

    x_arr(:,:) = r_arr*sqrt(sinth2(:,:))  !for theta as co_latitude angle
    y_arr(:,:) = r_arr*costh(:,:)

    !define scale factors from mu coordinate to an "index" coordinate system along field line

    Allocate(dmudx0(0:num_u3-1))
    dmudx0(:) = (cshift(mu,shift = 1,dim = 1)-cshift(mu,shift=-1,dim = 1))/d32
    dmudx0(0) = (-3.0*mu(0)+4.0*mu(1)-mu(2))/d32 !backwards difference
    dmudx0(num_u3-1) = (3.0*mu(num_u3-1)-4.0*mu(num_u3-2)+mu(num_u3-3))/d32    !forward difference

    do ii = 0, num_u3-1
            dmudx(:,ii) = dmudx0(ii)
    end do

    Allocate(bfac(0:num_u1-1,0:num_u3-1))
    Allocate(bsqrt(0:num_u1-1,0:num_u3-1))
    bfac(:,:) = 1.0+3.0*costh2(:,:)
    bsqrt(:,:) = sqrt(bfac(:,:))

    !Defining constants used in the model

    Allocate(z_arr(0:num_u1-1,0:num_u3-1))
    z_arr(:,:) = (r_arr-1.0)*re         !heights from earth surface (m)
    d3 = -1.0

    Allocate(g11(0:num_u1-1,0:num_u3-1))
    Allocate(gsup11(0:num_u1-1,0:num_u3-1))
    Allocate(h1(0:num_u1-1,0:num_u3-1))
    g11(:,:) = (r_arr(:,:)/costh0(:,:))**4/r_iono**2/bfac(:,:)**2*((1.0-r_iono/r_arr(:,:))**2+costh2(:,:)*bfac(0,0)**2*0.25/sinth2(:,:))
    gsup11(:,:) = r_iono**2/r_arr(:,:)**4*sinth2(:,:)*bfac(:,:)
    h1 = 1.0/sqrt(gsup11(:,:))

    Allocate(g22(0:num_u1-1,0:num_u3-1))
    Allocate(gsup22(0:num_u1-1,0:num_u3-1))
    Allocate(h2(0:num_u1-1,0:num_u3-1))
    g22 = r_arr(:,:)**2*sinth2(:,:)
    gsup22 = 1.0/g22(:,:)
    h2 = sqrt(g22(:,:))

    Allocate(g13(0:num_u1-1,0:num_u3-1))
    Allocate(gsup13(0:num_u1-1,0:num_u3-1))
    g13 = dmudx(:,:)*r_arr(:,:)**4*costh(:,:)*0.5/r_iono**2/costh0(:,:)/bsqrt(:,:)**2
    gsup13 = -0.5*r_iono**3/r_arr(:,:)**5*sinth02(:,:)*costh(:,:)/costh0**3*bfac/dmudx(:,:)


    Allocate(g33(0:num_u1-1,0:num_u3-1))
    Allocate(gsup33(0:num_u1-1,0:num_u3-1))
    Allocate(h3(0:num_u1-1,0:num_u3-1))
    Allocate(h30(0:num_u1-1,0:num_u3-1))
    g33 = (dmudx(:,:)*r_arr(:,:)**3/r_iono**2*costh0/bsqrt)**2
    gsup33 = r_iono**4/r_arr(:,:)**6/costh0**6*(costh2(:,:)*bfac(:,:)**2*0.25 + sinth2(:,:)*(1.0-r_iono/r_arr(:,:))**2)/dmudx(:,:)**2
    h3 = dmudx(:,:)*r_arr(:,:)**3/r_iono**2*costh0/bsqrt(:,:)
    h30 = abs(r_arr(:,:)**3/r_iono**2*costh0/bsqrt(:,:))

    Allocate(jac(0:num_u1-1,0:num_u3-1))
    jac = dmudx(:,:)*r_arr(:,:)**6*costh0/r_iono**3/bfac(:,:)

    Allocate(hratm_N(0:num_u1-1))
    Allocate(hthatm_N(0:num_u1-1))
    Allocate(hphiatm_N(0:num_u1-1))
    Allocate(hthg_N(0:num_u1-1))
    Allocate(hphig_N(0:num_u1-1))

    Allocate(hratm_S(0:num_u1-1))
    Allocate(hthatm_S(0:num_u1-1))
    Allocate(hphiatm_S(0:num_u1-1))
    Allocate(hthg_S(0:num_u1-1))
    Allocate(hphig_S(0:num_u1-1))

    !Atmospheric Scale Factors
    hratm_N = -2.0*r_iono*(costh0(:,0))**3/(costh(:,0)*bfac(:,0))*dmudx(:,0)
    hthatm_N = -r_iono*costh(:,0)/(2.0*sinth(:,0)*costh02(:,0))
    hphiatm_N = r_iono*sinth(:,0)

    !scale factors mapped to ground in spherical coords
    hthg_N = hthatm_N/r_iono
    hphig_N = hphiatm_N/r_iono

    hratm_S = -2.0*r_iono*(costh0(:,Num_u3-1))**3/(costh(:,Num_u3-1)*bfac(:,Num_u3-1))*dmudx(:,Num_u3-1)
    hthatm_S = -r_iono*costh(:,Num_u3-1)/(2.0*sinth(:,Num_u3-1)*costh02(:,Num_u3-1))
    hphiatm_S   = r_iono*sinth(:,Num_u3-1)

    !scale factors mapped to ground in spherical coords
    hthg_S = hthatm_S/r_iono
    hphig_S = hphiatm_S/r_iono

    Allocate(dels_arr(0:num_u3-1))
    Allocate(outerlength(0:num_u3-1))
    dels_arr = sqrt((cshift(x_arr(num_u1-1,:),-1)-cshift(x_arr(num_u1-1,:),0))**2 + &
      (cshift(y_arr(num_u1-1,:),-1)-cshift(y_arr(num_u1-1,:),0))**2)
    dels_arr(0) = 0.0
    outerlength = sum_cumulative(dels_arr)


    write(*,*) 'Grid Generation Finished!!!'


end subroutine gen_grid

subroutine New_r(r0,u1,u3,costh02,RI,ans,num_u3_half)
    implicit none
    integer:: num_u3_half
    double precision,dimension(0:num_u3_half-1) :: r0,ans,error,fr,df_dr,u3
    double precision:: u1,costh02,RI
    integer:: MaxN,ii
    double precision :: Tol

    MaxN = 100
    Tol = 1.e-3
    error(:) = 1.0
    ii = 0

    do while ((ii .lt. MaxN) .and. (maxval(error).gt.Tol))
            fr = u3(:)**2*r0(:)**4*costh02/RI**3-u1*r0(:)-RI
            df_dr = 4.0*u3**2*r0(:)**3*costh02/RI**3-u1
            ans(:) = r0(:)-fr(:)/df_dr(:)
            error(:) = abs(ans(:)-r0(:))
            r0(:) = ans(:)
            ii = ii+1;
    end do


end subroutine New_r

