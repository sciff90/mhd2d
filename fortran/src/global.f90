module global
      implicit none

      !Time Parameters
      double precision,parameter :: tmax = 300.0      !Max run time in Seconds
      double precision,parameter :: cour = 0.85     !fraction of courant number
      double precision,parameter :: dtplot = 1.0

      double precision,parameter :: m = 2.0         !azimuthal (wavenumber)
      double precision,parameter :: drive0 = 10.0e-9!amplitude of driver (T)
      double precision,parameter :: freq = 0.01     !frequency of driver in (Hz)

      !Directories
      character(*),parameter :: data_dir = "./data"     !data directory

      !Grid Parameters
      double precision,parameter :: Lval = 6.0      !conoical L value
      double precision,parameter :: LMin = 2.0      !min L
      double precision,parameter :: LMax = 10.0     !max L
      integer,parameter :: num_u1 = 151    !number of field lines
      double precision,parameter :: Re = 6378.388e3 !Re in metres
      double precision,parameter :: z0 = 80.0e3     !height of ionospheric thin sheet current (in) m
      double precision,parameter :: r_iono = 1.0+z0/re   !start altitude from earth centre
      !double precision,parameter :: ds0 = 5.0e3    !grid spacing (along cononical field line) at z0 (in m)
      double precision,parameter :: dsn = 500.0e3   !grid spacing (along cononical field line) at rtop (in m)
      double precision,parameter :: modefrac = 0.2  !fraction (of gridpoints in u1) to use as basis function in atmospheric expansions
      double precision,parameter :: d32 = 2.0 !"spacing" between the k+1 and k-1 points in the u3 direction

      !constants
      double precision,parameter :: pi = 3.141592
      !Scaled physical constants for unit length in Re

      double precision,parameter :: c_u = 2.997e8            !speed of light in vacuum
      double precision,parameter :: u0_u = 4.0*pi*1.0e-7     !magnetic permeability (SI units)
      double precision,parameter :: e0_u = 1.0/(u0_u*c_u**2) !dielectric constant for a vacuum in (si) units
      double precision,parameter :: c2 = (c_u/Re)**2         !speed of light squared in re^2/c^2
      double precision,parameter :: mhocgs = 1.0/(4.0*pi*e0_u) !conversion factor needed with bobs ionospheric file
      double precision,parameter :: re2 = re**2              !Re^2
      complex*16,parameter :: im = (0,m)

      !Global Variables
      integer :: num_u3,num_u3_half
      double precision :: del_u1,dt

      !Global Arrays
      double precision,dimension(:,:),allocatable :: x_arr,y_arr,z_arr,r_arr
      double precision,dimension(:),allocatable :: hratm_N,hthatm_N,hphiatm_N,hthg_N,hphig_N
      double precision,dimension(:),allocatable :: hratm_S,hthatm_S,hphiatm_S,hthg_S,hphig_S
      double precision,dimension(:,:),allocatable :: sigp_arr,sigh_arr,sig0_arr,eta_arr
      double precision,dimension(:),allocatable :: e1b1atm_N,e1b2atm_N,e2b1atm_N,e2b2atm_N
      double precision,dimension(:),allocatable :: e1b1atm_S,e1b2atm_S,e2b1atm_S,e2b2atm_S


      double precision,dimension(:,:),allocatable :: costh,va2_arr,eps_arr
      !Grid Arrays
      double precision,dimension(:,:),allocatable :: h3,h1,h2,bsqrt,g22,g11,g13,g21,g12,jac,g33
      double precision,dimension(:,:),allocatable :: gsup11,gsup13
      double precision,dimension(:),allocatable :: u1,outerlength
      !Factor Arrays
      double precision,dimension(:,:),allocatable:: gg12,gg21,e1e1,e1e2,e1f1,&
        e1f2,f1b2,e2e1,e2e2,e2f1,e2f2,f2b1,f2b3,e3b23,b1e2,b2e1,b2e3,b3e2,&
        e1esup1,e1e3,e3b21
      complex*16,dimension(:,:),allocatable :: f1b3,e3b1,e3b3,b1e3,b3e1,&
        psiib3_N,psiib3_S,psigb3_S,psigb3_N
      double precision :: dx2,dxsq


end module
