module global
      implicit none

      !Time Parameters
      real,parameter :: tmax = 1200.0      !Max run time in Seconds
      real,parameter :: cour = 0.85     !fraction of courant number
      real,parameter :: dtplot = 1.0

      real,parameter :: m = 2.0         !azimuthal (wavenumber)
      real,parameter :: phi = 22.5      !azimuthal location for plots
      real,parameter :: drive0 = 10.0e-9!amplitude of driver (T)
      real,parameter :: freq = 0.03     !frequency of driver in (Hz)

      !Directories
      character(*),parameter :: data_dir = "./data"     !data directory

      !Grid Parameters
      real,parameter :: Lval = 6.0      !conoical L value
      real,parameter :: LMin = 2.0      !min L
      real,parameter :: LMax = 10.0     !max L
      integer,parameter :: num_u1 = 151    !number of field lines
      real,parameter :: Re = 6378.388e3 !Re in metres
      real,parameter :: z0 = 80.0e3     !height of ionospheric thin sheet current (in) m
      real,parameter :: r_iono = 1.0+z0/re   !start altitude from earth centre
      real,parameter :: ds0 = 5.0e3    !grid spacing (along cononical field line) at z0 (in m)
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
      real :: del_u1,dt,nt

      !Global Arrays
      real,dimension(:,:),allocatable :: x_arr,y_arr,z_arr,r_arr
      real,dimension(:),allocatable :: hratm_N,hthatm_N,hphiatm_N,hthg_N,hphig_N
      real,dimension(:),allocatable :: hratm_S,hthatm_S,hphiatm_S,hthg_S,hphig_S
      real,dimension(:,:),allocatable :: sigp_arr,sigh_arr,sig0_arr,eta_arr
      real,dimension(:),allocatable :: e1b1atm_N,e1b2atm_N,e2b1atm_N,e2b2atm_N
      real,dimension(:),allocatable :: e1b1atm_S,e1b2atm_S,e2b1atm_S,e2b2atm_S


      real,dimension(:,:),allocatable :: costh,va2_arr,eps_arr,va_arr,rho_a
      complex,dimension(:),allocatable :: sigpatm_N,sighatm_N,sig0atm_N
      complex,dimension(:),allocatable :: sigpatm_S,sighatm_S,sig0atm_S
      !Grid Arrays
      real,dimension(:,:),allocatable :: h3,h1,h2,bsqrt,g22,g11,g13,g21,g12,jac,g33,h30
      real,dimension(:,:),allocatable :: gsup11,gsup13,gsup22,gsup33,dmudx
      real,dimension(:),allocatable :: u1,u3,outerlength
      !Factor Arrays
      real,dimension(:,:),allocatable:: gg12,gg21,e1e1,e1e2,e1f1,&
        e1f2,f1b2,e2e1,e2e2,e2f1,e2f2,f2b1,f2b3,e3b23,b1e2,b2e1,b2e3,b3e2,&
        e1esup1,e1e3,e3b21
      complex,dimension(:,:),allocatable :: f1b3,e3b1,e3b3,b1e3,b3e1,&
        psiib3_N,psiib3_S,psigb3_S,psigb3_N
      real :: dx2,dxsq


end module
