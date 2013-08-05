
Pro Get_Driver, drive0, freq, t, h3, h30, u3, num_u3_half, num_u1, num_u3, Re, Driver, OuterLength
;	run driver.pro to get desired spectral shape of time dependent part
;	Driver can be defined on ALL grid points on outer boundary


;	Amplitude
Amp			= drive0*(h3(num_u1-1,num_u3_half))*Re;/ABS(h3(num_u1-1,num_u3_half))				; Convert tesla into model scaled units

;	Spatial Distribution
half_width 	= 2.0										; in Re's
Z		= (outerlength - outerlength(Num_u3/2))/half_width
Spatial = exp(-Z^2)
;Spatial 	= ABS(reform(h3(num_u1-1,*))) * exp(-((u3-num_u3_half)/half_width)^2)   ; Spatial distribution of Driver

;	Time Dependance
;	Monocromatic Driver
omega		= 2*!dpi*freq
;tdep 		= sin(omega*t)

;	Bandlimited Driver
omega		= 2*!dpi*freq
Width		= 0.5*freq									; width of band (relative to central freq)
carr		= sin(omega*t)
Env			= exp(-((t-2.0/Width)*Width)^2.)
;tdep		= Carr*Env

;	Pulse
Dnorm		= 5.0
Decay		= 5.0
tdep		= t/DNorm*(1.0-t/DNorm)*Exp(-Decay*t/Dnorm)*2*Dnorm


;	The Driver
Driver				= Amp*spatial*tdep
Driver(0)			= 0.0
Driver(Num_u3-1)	= 0.0


end			; End of Driver

