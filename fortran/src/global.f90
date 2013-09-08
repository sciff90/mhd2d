module global
      implicit none

      !Time Parameters
      real,parameter :: tmax = 0.0      !Max run time in Seconds
      real,parameter :: cour = 0.85     !fraction of courant number
      
      real,parameter :: m = 2.0         !azimuthal (wavenumber)
      real,parameter :: drive0 = 10.0e-9!amplitude of driver (T)
      real,parameter :: freq = 0.01     !frequency of driver in (Hz)

      !Directories
      character(*),parameter :: data_dir = "./data"     !data directory
      
      !Grid Parameters
      real,parameter :: Lval = 6.0      !conoical L value
      real,parameter :: LMin = 2.0      !min L
      real,parameter :: LMax = 10.0     !max L
      integer,parameter :: num_u1 = 151    !number of field lines
      real,parameter :: Re = 6378.388e3 !Re in metres
      real,parameter :: z0 = 80.0e3     !height of ionospheric thin sheet current (in) m
      !real,parameter :: ds0 = 5.0e3    !grid spacing (along cononical field line) at z0 (in m)
      real,parameter :: dsn = 500.0e3   !grid spacing (along cononical field line) at rtop (in m)
      real,parameter :: modefrac = 0.2  !fraction (of gridpoints in u1) to use as basis function in atmospheric expansions
      real,parameter :: d32 = 2.0 !"spacing" between the k+1 and k-1 points in the u3 direction

      !constants
      real,parameter :: pi = 3.141592
      !Scaled physical constants for unit length in Re
      
      real,parameter :: c_u = 2.997e8            !speed of light in vacuum
      real,parameter :: u0_u = 4.0*pi*1.0e-7     !magnetic permeability (SI units)
      real,parameter :: e0_u = 1.0/(u0_u*c_u**2) !dielectric constant for a vacuum in (si) units
      real,parameter :: c2 = (c_u/Re)**2         !speed of light squared in re^2/c^2
      real,parameter :: mhocgs = 1.0/(4.0*pi*e0_u) !conversion factor needed with bobs ionospheric file 
      real,parameter :: re2 = re**2              !Re^2
      complex,parameter :: im = (0,m)

      !Global Variables
      integer :: num_u3,num_u3_half
      real :: del_u1

      !Global Arrays
      real,dimension(:,:),allocatable :: x_arr,y_arr,z_arr,r_arr
      real,dimension(:,:),allocatable :: costh
      !Grid Arrays
      real,dimension(:,:),allocatable :: h3,h1,h2

end module
