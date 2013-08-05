subroutine gen_grid()
    use global 
    implicit none
    
    !parameters
    real,parameter :: r_iono = 1.0+z0/re   !start altitude from earth centre
    real,parameter :: k0 = 8.02e15      !Magnetic moment in Wb.m
    real,parameter :: cosColat0 = sqrt(1.0-(r_iono/LVal))
    real,parameter :: r0 = r_iono*Re
    real,parameter :: Lat0 = asin(cosColat0)
    real,parameter :: ColLat0 = acos(cosColat0)
    real,parameter :: DelS = -0.1e3 !approx 100m spacing. Step along field line (- going outwards from northern hemisphere)
    real,parameter :: d32 = 2.0 !"spacing" between the k+1 and k-1 points in the u3 direction

    !Variables
    real :: X,Z,r,Lat,Bx,Bz,Bx_unit,Bz_unit
    real:: Sn,dsi,ccosth0 
    integer :: num_u3_half,num_u3,ii,jj,temp
    integer*4 :: Npts
    !Arrays
    real,dimension(:),allocatable :: XX,ZZ,CLat,u3,ds,s,ss,y,cosCLat,cr,ccosth 
    real,dimension(:),allocatable :: mu
    
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
                    s(ii) = sum(ds(0:ii))
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
    write (*, *) "Npts = ",Npts
    write (*, *) "Numu3 = ",num_u3
end subroutine gen_grid
