
;	This is a rewrite of the "alfdip2dv3_mds_Yee.pro" code (from Minnesota visit).
;	Full meridional plane
;	Using a Yee grid and with Height distributed Ionosphere and tilt dipole coordinates
;	MDS - 14th May 2012
;	MDS - 5th June 2012 - Beta version code given to CLW
;	MDS - 14th June 2012 - Code for Anlytic Ionosphere added (I hope)
;	v2 MDS - 3rd July 2012 - Wave speed modified by collisions included
;	v3 MDS - 6th July 2012 - Mod Odds grid points in Num_u1
;	v4 MDS - 30th Oct 2012 - Only using B3 pts to calculate Psi's (not average then fit)
;	v5 MDS - 7th Dec 2012 - Added two bases function set to fit (not average then fit)
;	v6 MDS - 21st Feb 2013 - New (improved) Grid generation (fix for cononical field line length)
Resolve_all
;	Time Parameters
tmax		= 10.0	             	; sets max run time [in sec]
cour		= 0.85					; Fraction of Courant number to use in time step
dtplot		= 1.0               	; plot every...X seconds
m			= 2.0               	; azimuthal variation (wave number)
Phi			= 22.5					; azimuthal location for plots
drive0		= 10.0e-9             	; amplitude of driver (T)
Freq  		= 0.010					; Frequency of driver (in Hz)

;	Directories
Pth			= '/media/Raid_Data/Data/mhd2d/idl/data/'					; Directory for data file
out_pth		= '/media/Raid_Data/Data/mhd2d/idl/plots/'					; Directiory for images
inp_pth		= '../iono_data/'				; Directory for Ionopsheric data (eg. Conductances etc)
Neutral_file= inp_pth+'NeutralAtm_Min_Mod.txt'					; Neutral Atmosphere File
plot_png=0												; 1= plot png files
Plot_fields =0												; 1 = plot e,b fields each dtplot time step

;	Grid Parameters
LVal		= 6.0             		; Cononical L Value
LMin		= 2.0             		; min L
LMax		= 10.0            		; max L
Num_u1		= 151           		; number of field lines (should be even)
Re			= 6378.388e3          	; Re in m

z0			= 80.0e3             	; height of ionospheric thin sheet current  (in m)
;ds0			= 5.0e3              	; grid spacing (along cononical field line) at z0 (in m)
dsn			= 500.0e3            	; grid spacing (along cononical field line) at rtop (in m)

modefrac	= 0.20					; Fraction (of Grid points in u1) to use as basis function in atmospheric expansions

;  cpu, tpool_nthreads=1, tpool_min_elts=1000				; Multi threading

;	Color Table with white back ground
Device,decomposed=0
loadct,4,/silent
TVLCT,r,g,b,/get
r(0)=255 & g(0)=255 & b(0)=255
r(255)=0 & g(255)=0 & b(255)=0
TVLCT,r,g,b

;	Scaled physical constants where unit for length is in Re
c_u 	= 2.997d8               	; Speed of light in a Vacuum (in SI units)
u0_u 	= 4.0d0*!dpi*1.0d-7	   		; Magnetic Permeability (in SI units)
e0_u 	= 1.0d0/(u0_u*c_u^2)   		; Dielectric Constant for a vacuum (in SI units)

c2 		= (c_u/(Re))^2				; Speed of light squared in Re^2/common s^2
mhocgs	= 1.0d0/(4.0d0*!dpi*e0_u)	; Conversion factor need with Bob's ionospheric file (will remove once we have an SI version)
re2		= re*re						; Re^2 [m^2]
im		= Dcomplex(0.0d0,m)			;

;	Generate Grid and Metrics
Print,'Generating Grid...'
Gen_grid, 		LMin, LVal, LMax, num_u1, num_u3, num_u3_half, del_u1, d32, $
            	ds0, dsn, z0, x_arr, y_arr, r_arr, u1, dmudx, h1, h2, h3, h30, $
            	g11, g13, g22, g33, gsup11, gsup13, gsup22, gsup33, jac, $
            	sinth0, costh0, costh02, bfac, cr, costh, plot_png, out_pth, $
            	hthatm_N, hphiatm_N, hth_g_N, hphi_g_N, hratm_N, $
            	hthatm_S, hphiatm_S, hth_g_S, hphi_g_S, hratm_S, Re, $
            	hthg_N, hthg_S, hphig_N, hphig_S, z_arr, r_iono, u3, d3, bsqrt, OuterLength


;	Get/creating Va, Conductivities./Conductances
Print,'Generating Va and Conductance profiles...'
Get_Va,  		c_u, e0_u, u0_u, x_arr, y_arr, r_arr, z_arr, costh, Re, Num_u1, num_u3, plot_png, out_pth, Va_arr, eps_arr, rho_a, $
				cour, d32,  del_u1, h1, h2, h3, m, nt, tmax, dt, va2_arr, No_tot, $
		 		Neutral_file, num_u3_half, z0, mhocgs, hphiatm_N, hphiatm_S, hthatm_N, hthatm_S, bsqrt, $
				eta_arr, sig0_arr, sigp_arr, sigh_arr, sigpatm_N, sighatm_N, sig0atm_N, sigpatm_S, sighatm_S, sig0atm_S, $
				e1b1atm_N, e1b2atm_N, e2b1atm_N, e2b2atm_N, e1b1atm_S, e1b2atm_S, e2b1atm_S, e2b2atm_S

Print,'Generating Basis Functions...'
Get_BasisFn,	Num_u1, m, u1, del_u1, modefrac, r_iono, psiib3_N, psiib3_S, plot_png, out_pth, $
			 	hratm_N, hratm_S, psigb3_N, psigb3_S, ixbc_0, ixbc_n

Print,'Generating Numerical Factors...'
Get_Facts,		sigp_arr, sigh_arr, eps_arr, eta_arr, g11, g22, g33, g13, gsup11, gsup13, gsup33, h1, h2, h3, Jac, dt, $
				va_arr, va2_arr, im, d32, del_u1, $
				e1e1, e1e2, e1f1, e1f2, f1b2, f1b3, e2e1, e2e2, e2f1, e2f2, f2b1, f2b3, e3b1, e3b21, e3b23, e3b3, $
				b1e2, b1e3, b2e1, b2e3, b3e1, b3e2, e1esup1, e1e3, dx2, dxsq, z_arr, r_arr, Plot_PNG, out_pth

Print,'Saving Parameters'
filenme = 'Faligned2D_Parameters.sav'
save,filename=out_pth+filenme, $
				sigp_arr, sigh_arr, eps_arr, eta_arr, g11, g22, g33, g13, gsup13, gsup33, h1, h2, h3, h30, Jac, dt, tmax, nt,$
				va_arr,	nt, LMin, LVal, LMax, num_u1, num_u3, x_arr, y_arr, z_arr, r_arr, u1, u3, m,  $
				r_iono, rho_a, cour, dtplot, drive0, Freq, eta_arr, sigp_arr, sigh_arr, $
				sigpatm_N, sighatm_N, sig0atm_N, sigpatm_S, sighatm_S, sig0atm_S, z0, ds0, dsn, modefrac
;				cmax,dtx, dtz

Print,'Starting Time Iterations...'
Iterate,   		num_u1,num_u3, num_u3_half, dt, va_arr, h1, h2, h3, h30, va2_arr, d32, del_u1, $
				g11,g22, g33, g13, u1, u3,gsup11,gsup22,gsup33, gsup13, out_pth, im, Phi, $
				e1e1, e1e2, e1f1, e1f2, f1b2, f1b3, e2e1, e2e2, e2f1, e2f2, f2b1, f2b3, e3b1, e3b21, e3b23, e3b3, $
				b1e2, b1e3, b2e1, b2e3, b3e1, b3e2, e1esup1, e1e3, dx2, dxsq, tmax, dtplot, $
				e1b2atm_N,e1b1atm_N, e2b1atm_N, e2b2atm_N, e1b2atm_S,e1b1atm_S, e2b1atm_S, e2b2atm_S, psiib3_N, psiib3_S, $
				psigb3_N, psigb3_S, drive0, Freq, Plot_fields, x_arr, y_arr, hratm_N, hratm_S ,$
				hphiatm_N, hphiatm_S, hthatm_N, hthatm_S,hthg_N, hthg_S, hphig_N, hphig_S, Re, ixbc_0, ixbc_n, OuterLength

stop
END				; End Main program

