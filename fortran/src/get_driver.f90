subroutine get_driver(tt,u3_in,driver)
  use global
  implicit none

  real :: tt
  real,dimension(:),intent(in) :: u3_in
  real,dimension(0:num_u3-1),intent(out) :: driver

  !Local Definitions
  real :: amp,half_width,omega,width,carr,Env
  real :: dnorm,decay,tdep
  real,dimension(0:num_u3-1) :: z,spatial

  !Amplitude
  !Convert tesla into model scaled units
  Amp = drive0*(h3(num_u1-1,num_u3_half))*Re  !/ABS(h3(num_u1-1,num_u3_half))
  !spatial Distribution
  half_width = 2.0  !in Re's
  Z = (outerlength - outerlength(Num_u3/2))/half_width
  Spatial = exp(-Z**2)
  !Spatial 	= ABS(reform(h3(num_u1-1,*))) * exp(-((u3-num_u3_half)/half_width)^2)   ; Spatial distribution of Driver

  !Time Dependance
  !Monocromatic Driver
  omega = 2*pi*freq
  !tdep 		= sin(omega*t)

  !Bandlimited Driver
  omega = 2*pi*freq
  Width = 0.5*freq  !width of band (relative to central freq)
  carr = sin(omega*tt)
  Env = exp(-((tt-2.0/Width)*Width)**2.)
  !tdep		= Carr*Env

  !Pulse
  Dnorm = 5.0
  Decay = 5.0
  tdep = tt/DNorm*(1.0-tt/DNorm)*Exp(-Decay*tt/Dnorm)*2*Dnorm

  !The Driver
  Driver = Amp*spatial*tdep
  !test driver
  !driver = 1.0
  Driver(0) = 0.0
  Driver(Num_u3-1) = 0.0

end subroutine get_driver
