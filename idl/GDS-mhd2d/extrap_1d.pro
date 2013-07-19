pro Extrap_1D, array,grid,location
;	Extrapolation to end point by second order Taylor expansion

;Extrap_1D,Bsup_n,u1,0

Noelm 	= N_elements(grid)

If location Eq 0 then $				; forwards extrap
begin
 hh				= grid(1)-grid(2) 																							; step size along ionosphere is uniform
 FirD 			= ( array(3) - 4.0*array(2) + 3.0*array(1))/(2.0*hh)													; First Derivative Forwards difference O(h^2)
 SecD 			= (-array(4) + 4.0*array(3) - 5.0*array(2) + 2.0*array(1))/(hh^2) 										; Second Derivative Forwards difference O(h^2)
 hh1			= grid(0) - grid(1)
 array(0)		= array(1) + hh1*FirD + ((hh1^2)/2.0)*SecD
end

If location Ne 0 then $				; backwardsward extrap
begin
 hh				= grid(Noelm-2)- grid(Noelm-3) 																				; step size along ionosphere is uniform on B1_S grid
 FirD 			= ( array(Noelm-4) - 4.0*array(Noelm-3) + 3.0*array(Noelm-2))/(2.0*hh)									; First Derivative Backwards difference O(h^2)
 SecD 			= (-array(Noelm-5) + 4.0*array(Noelm-4) - 5.0*array(Noelm-3) + 2.0*array(Noelm-2))/(hh^2) 				; Second Derivative Backwards difference O(h^2)
 hh1			= grid(Noelm-1) - grid(Noelm-2)																				; Step size on B1_S_interp grid
 array(Noelm-1)	= array(Noelm-2) + hh1*FirD + ((hh1^2)/2.0)*SecD
end

end

