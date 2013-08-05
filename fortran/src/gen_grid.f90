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
    real:: Sn,dsi 
    integer :: num_u3_half,num_u3,ii   
    integer*4 :: Npts
    !Arrays
    real,dimension(:),allocatable :: XX,ZZ,CLat,u3,ds
    
    !start grid_gen
    X = r0*cos(Lat0)    
    Z = r0*sin(Lat0)
    Npts = 0

    !Calculate cononical field line length using analytic dipole field
    do while (Z .gt. 0)
        r = sqrt(X**2+Z**2)
        Lat = atan(Z,X)
        Bx = -3*k0/r**3*sin(Lat)*cos(Lat)
        Bz = -2*k0/r**3*cos(Lat)

        Bx_unit = Bx/sqrt(Bx**2+Bz**2)
        Bz_unit = Bz/sqrt(Bx**2+Bz**2)

        X = X + Bx_unit + DelS 
        Z = Z + Bz_unit + DelS 
        Npts = Npts+1
    end do

    Allocate(XX(Npts))
    Allocate(ZZ(Npts))
    Allocate(CLat(Npts))

    X = r0*cos(Lat0)
    Z = r0*sin(Lat0)
    Npts = 0

    do while (z .gt. 0)
            r = sqrt(X**2+Z**2)
            Lat = atan(Z,X)
            XX(Npts) = X
            ZZ(Npts) = Z
            CLat(Npts) = pi/2 -Lat 

            Bx = -3*k0/r**3 * sin(Lat)*cos(Lat)
            Bz = -2*k0/r**3*sin(Lat)**2+1*k0/r**3*cos(Lat)**2

            Bx_unit  = Bx/sqrt(Bx**2+Bz**2)
            Bz_unit  = Bx/sqrt(Bx**2+Bz**2)

            X = X + Bx_unit*DelS 
            Z = Z + Bz_unit*DelS 
            Npts = Npts + 1
    end do

    Sn = real(Npts-1)*abs(DelS)
    num_u3_half = nint(2.0*Sn/(dsn))   
    num_u3 = 2*num_u3_half + 1  !number of points must be odd

    do ii = 0, num_u3-1         !u3 ordiniate is just index of array
            u3(ii) = ii
    end do

!    d32 = 2.0                   !the spacing k+1 and k-1 points in u3 direction
!    dsi = 2.0*Sn/(num_u3_half*(num_u3_half+1)) !space coming from a sum of the arithmetic progression Sn = N/2*(2ds0+(N-1)dsi) and ds0=dsi
    
!    do ii = 0, num_u3_half-1
!        num_u3_half(ii) = dsi*(1+ii)   
!    end do

!    s = ??
    
!    do ii = 0, Npts-1
!            ii*abs(Dels)
!    end do

    !Generate points along Cononical (LVal) field line
!    cosCLat = cos(Clat)
    
!    Allocate(y(0:num_u3-1))
!    y(0) = cosCLat(0)
!    do ii = 1, num_u3_half-1
            ! code
!    end do
     

    write (*, *) "Npts = ",Npts
    write (*, *) "Numu3 = ",num_u3
end subroutine gen_grid
