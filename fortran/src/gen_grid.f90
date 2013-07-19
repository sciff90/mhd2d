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
    real:: Sn,num_u3_half,num_u3,dsi 
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


    write (*, *) "Npts = ",Npts
end subroutine gen_grid
