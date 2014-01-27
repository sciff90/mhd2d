Pro mhd2d_idl
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

;	Time Parameters
tmax		= 0.0	             	; sets max run time [in sec]
cour		= 0.85					; Fraction of Courant number to use in time step
dtplot		= 1.0               	; plot every...X seconds
m			= 2.0               	; azimuthal variation (wave number)
Phi			= 22.5					; azimuthal location for plots
drive0		= 10.0e-9             	; amplitude of driver (T)
Freq  		= 0.010					; Frequency of driver (in Hz)

;	Directories
Pth			= '/home/gareths/Data/mhd2d/idl/data/'					; Directory for data file
out_pth		= '/home/gareths/Data/mhd2d/idl/plots/'					; Directiory for images
inp_pth		= '/home/gareths/Code/mhd2d/idl/iono_data/'				; Directory for Ionopsheric data (eg. Conductances etc)
Neutral_file= inp_pth+'NeutralAtm_Min_Mod.txt'					; Neutral Atmosphere File
plot_png 	= 1												; 1= plot png files
Plot_fields = 1												; 1 = plot e,b fields each dtplot time step

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

c2 		= (c_u/(Re))^2				; Speed of light squared in Re^2/s^2
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

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Pro Gen_grid,	LMin, LVal, LMax, num_u1, num_u3, num_u3_half, del_u1, d32, $
            	ds0, dsn, z0, x_arr, y_arr, r_arr, u1, dmudx, h1, h2, h3, h30, $
            	g11, g13, g22, g33, gsup11, gsup13, gsup22, gsup33, jac, $
            	sinth0, costh0, costh02, bfac, cr, costh, plot_png, out_pth, $
            	hthatm_N, hphiatm_N, hth_g_N, hphi_g_N, hratm_N, $
            	hthatm_S, hphiatm_S, hth_g_S, hphi_g_S, hratm_S, Re, $
            	hthg_N, hthg_S, hphig_N, hphig_S, z_arr, r_iono, u3, d3, bsqrt, OuterLength

; d1 	-> del_u1
; to generate grid for 2D - with E parallel
; LVal             	; cononical L Value
; LMin             	; min L
; LMax             	; max L
; num_u1           	; num of field lines
; ds0				; dist between grid points at ionos
; dsn 				; dist between grid points at equator plane

; num_u1	= fix(num_u1/2)*2					; Must be even (number of field lines)
 num_u1	= fix(num_u1/2)*2+1					; Must be even (number of field lines)


; Define fixed parameters
;re		= 6378.388    			; Re in km
r_iono		= 1.0d0+z0/re          	; start altitude in Re from Earth centre -> 1.015678


k0			= 8.02e15						; Magnetic moment? in Wb.m
cosColat0	= sqrt(1.0-(r_iono/LVal))
r0 			= r_iono*Re
Lat0		= Asin(cosColat0)
CoLat0		= Acos(cosColat0)
DelS		= -0.1d3						; ~ 100m spacing space! Step along field line (- going outwards from northern hemisphere)

X 			= r0*cos(Lat0)
Z 			= r0*sin(Lat0)
Npts 		= Long(0)

;	Calculate cononical field line length using Analytic dipole field

Repeat Begin

r 	= sqrt(X^2+Z^2)
Lat	= atan(Z,X)

Bx 	= -3.0*k0/r^3 * sin(Lat) * cos(Lat)
Bz	= -2.0*k0/r^3 * sin(Lat)^2 + 1.0*k0/r^3 * cos(Lat)^2

Bxunit	= Bx/sqrt(Bx^2 + Bz^2)
Bzunit	= Bz/sqrt(Bx^2 + Bz^2)

X = X + Bxunit*DelS
Z = Z + Bzunit*DelS

Npts++

ENDREP UNTIL Z le 0.0


XX = Dblarr(Npts)
ZZ = Dblarr(Npts)


CLat = Dblarr(Npts)

X = r0*cos(Lat0)
Z = r0*sin(Lat0)

Npts = Long(0)

Repeat Begin

r 	= sqrt(X^2+Z^2)
Lat	= atan(Z,X)

xx(Npts) = X
zz(Npts) = Z
CLat(Npts)	= !dpi/2.0-Lat


Bx 	= -3.0*k0/r^3 * sin(Lat) * cos(Lat)
Bz	= -2.0*k0/r^3 * sin(Lat)^2 + 1.0*k0/r^3 * cos(Lat)^2

Bxunit	= Bx/sqrt(Bx^2 + Bz^2)
Bzunit	= Bz/sqrt(Bx^2 + Bz^2)

X = X + Bxunit*DelS
Z = Z + Bzunit*DelS

Npts++

ENDREP UNTIL Z  le 0.0 ;OR sqrt( X^2  + Z^2 ) gt 100.0 ;

Sn 			= Double(Npts-1)*ABS(DelS)

num_u3_half	= (floor(2.0*sn/(dsn)))								;
num_u3 		= 2*num_u3_half + 1			; number points along field line (must be odd)
u3			= findgen(num_u3)			; u3 ordinate is just the index of the array!
d32 		= 2.0               		; the "spacing" between the k+1 and k-1 (2) points in the u3 direction

dsi			= 2.0*Sn/(Num_u3_half*(Num_u3_half+1.0))		; Space coming from a sum of the arithmetic progression Sn = N/2(2ds0+(N-1)dsi) and ds0 = dsi
ds			= dsi + Dindgen(num_u3_half)*dsi

s			= [0,total(ds,/cumulative,/double)]
ss			= Findgen(Npts)*ABS(Dels)

; generate points along Cononical (LVal) field line
CosCLat 	= cos(Clat)
y			= dblarr(num_u3)
y[0]		= CosCLat(0)
for i=1,num_u3_half-1 do begin
 ii				= max(where(ss lt s[i]))
 print,ii
 y[i]			= CosCLat[ii]+(CosCLat[ii]-CosCLat[ii-1])*(s[i]-ss[ii])/(ss[ii]-ss[ii-1])
 y[num_u3-1-i]	= -y[i]       					;	fill other hemis
endfor
 y[num_u3-1] 	= -y(0)
 y[num_u3_half]	= 0.0          						; equator point

; 'c' prefix for canonical field line
ccosth	= y											; Colats of points along cononical field line
ccosth0	= y(0)
cr		= LVal*(1.0d0-y*y)							; Cononical field line Radial Distance from centre of earth (in Re)
mu		= r_iono*r_iono*ccosth/cr/cr/ccosth0      	; Dipole mu coordinate for cononical field line
print,num_u3,' cells, rtop=',max(cr)         		; rmax

; set up Nu grid
numin	= -r_iono/LMin								; Starting Colatitude
numax	= -r_iono/LMax
del_u1	= (numax-numin)/double(num_u1-1)
u1		= findgen(num_u1)*del_u1 + numin

r_arr 	= dblarr(num_u1,num_u3)
sinth2	= dblarr(num_u1,num_u3)
sinth02	= dblarr(num_u1,num_u3)
costh 	= dblarr(num_u1,num_u3)
dmudx	= dblarr(num_u1,num_u3)

for i=0,num_u1-1 do sinth02[i,*]=replicate(-u1[i],num_u3)

costh02	= 1.0 - sinth02								; Ionospheric Colat of field lines squared
costh0	= sqrt(costh02)
sinth0	= sqrt(sinth02)                          	; odd ftn

; Northern hemisphere
for i=0,num_u1-1 do begin
 rinit=r_iono/sqrt(abs(mu[0:num_u3_half-1]*costh0[i,0]))
 New_r,rinit,u1[i],mu[0:num_u3_half-1],costh02[i,0],r_iono,ans
 r_arr[i,0:num_u3_half-1] = ans
 sinth2[i,0:num_u3_half-1] = -u1[i]*r_arr[i,0:num_u3_half-1]/r_iono
endfor

; Extend to Southern hemisphere
For i=0,num_u3_half-1 do begin
 r_arr[*,num_u3-i-1] = r_arr[*,i]
 sinth2[*,num_u3-i-1] = sinth2[*,i]
endfor

r_arr[*,num_u3_half] 	= -r_iono/u1     			; equator point
sinth2[*,num_u3_half] 	= 1.0

costh2 	= 1.0d0 - sinth2
costh[*,0:num_u3_half-1] 		=          sqrt(costh2[*,0:num_u3_half-1])                 	; Northern Hemis
costh[*,num_u3_half:num_u3-1] 	= -1.0d0 * sqrt(costh2[*,num_u3_half:num_u3-1])   			; Southern Hemis
sinth 	= sqrt(sinth2)						; Theta only runs 0<x<pi so sin(theta) > 0

x_arr	= r_arr*sqrt(sinth2)                ; for theta as co_Latitude angle
y_arr	= r_arr*costh

;define scale factors from mu coordinate to an "index" coordinate system along field line
dmudx0			= (shift(mu,-1)-shift(mu,1))/d32                 			; diff -> dmudx0(0) and dmudx0(num_u3-1) need to be calc from mu
dmudx0[0]		=(-3.0*mu[0] + 4.0*mu[1]-mu[2])/d32          				; back diff
dmudx0[num_u3-1]=(3.0*mu[num_u3-1] - 4.0*mu[num_u3-2] + mu[num_u3-3])/d32  	; forward diff
for k=0,num_u3-1 do dmudx[*,k]=replicate(dmudx0[k],num_u1)

bfac 	= 1.0d0 + 3.0d0*costh2
bsqrt 	= sqrt(bfac)

;	Defining constants etc used in model
z_arr 	= (r_arr-1.0d0)*re          ; height from Earth surface (m)
r_iono	= 1.0d0+z0/re            	; start altitude in Re from Earth centre
u3 		= Findgen(Num_u3)
d3 		= -1.d0

g11 	= (r_arr/costh0)^4/r_iono^2/bfac^2*((1.0d0-r_iono/r_arr)^2+costh2*bfac[0]^2*0.25/sinth2)
gsup11 	= r_iono^2/r_arr^4*sinth2*bfac
h1 		= 1.0/sqrt(gsup11)

g22 	= r_arr*r_arr*sinth2
gsup22 	= 1.0d0/g22
h2 		= sqrt(g22)

g13 	= dmudx*r_arr^4*costh*0.5/r_iono^2/costh0/bsqrt^2
gsup13 	= -0.5*r_iono^3/r_arr^5*sinth02*costh/costh0^3*bfac/dmudx

;g33 	= h3*h3
g33 	= (dmudx*r_arr^3/r_iono^2*costh0/bsqrt)^2
gsup33 	= r_iono^4/r_arr^6/costh0^6*(costh2*bfac^2*0.25  +sinth2*(1.0-r_iono/r_arr)^2)/dmudx^2
h3 		= dmudx*r_arr^3/r_iono^2*costh0/bsqrt
h30		=   abs(r_arr^3/r_iono^2*costh0/bsqrt)   									; h_mu, used to rotate bmu

jac		= dmudx*r_arr^6*costh0/r_iono^3/bfac

;	Atmospheric Scale factors
hratm_N 	= -2.0*r_iono*(costh0[*,0])^3/(costh[*,0]*bfac[*,0])*dmudx[*,0]
hthatm_N  	= -r_iono*costh[*,0]/(2.0*sinth[*,0]*costh02[*,0])
hphiatm_N 	= r_iono*sinth[*,0]
;scale factors mapped to ground in spherical coords
hthg_N 		= hthatm_N/r_iono
hphig_N 	= hphiatm_N/r_iono


hratm_S		= -2.0*r_iono*(costh0[*,Num_u3-1])^3/(costh[*,Num_u3-1]*bfac[*,Num_u3-1])*dmudx[*,Num_u3-1]
hthatm_S  	= -r_iono*costh[*,Num_u3-1]/(2.0*sinth[*,Num_u3-1]*costh02[*,Num_u3-1])
hphiatm_S 	= r_iono*sinth[*,Num_u3-1]
;scale factors mapped to ground in spherical coords
hthg_S 		= hthatm_S/r_iono
hphig_S		= hphiatm_S/r_iono

device,decomposed=0
;loadct,0
!p.multi=0
window,0,xsize=800,ysize=700,title='Grid'
!P.charsize=1.6
!P.Charthick=1.2
!p.symsize=0.2

plot,x_arr,y_arr,psym=2,title='2D Faligned Grid',xtitle='[R!dE!n]',ytitle='[R!dE!n]',$
	 xstyle=1,ystyle=1,/isotropic

If Plot_png EQ 1 then $
begin
 Filename 	= out_pth+'2D_Faligned_Grid.png'
 image 		= TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
end

device,decompose=0
window,1,Title='Title Grid',Xsize=1200,Ysize=700
!P.multi=[0,4,3]
!P.charsize=2.5
!P.Charthick=1.2
Lidx = min(where (r_arr(*,num_u3_half) GE Lval))

plot,u3,gsup11(Lidx,*),Title='g11_con',Xtitle='u3',xstyle=1
plot,u3,gsup22(Lidx,*),Title='g22_con',Xtitle='u3',xstyle=1
plot,u3,gsup33(Lidx,*)*dmudx(Lidx,*)^2,Title='g33_con',Xtitle='u3',xstyle=1
plot,u3,gsup13(Lidx,*)*dmudx(Lidx,*),Title='g13_con',Xtitle='u3',xstyle=1

plot,u3,g11(Lidx,*),Title='g11_cov',Xtitle='u3',/ylog,xstyle=1
plot,u3,g22(Lidx,*),Title='g22_cov',Xtitle='u3',/ylog,xstyle=1
plot,u3,g33(Lidx,*)/dmudx(Lidx,*)^2,Title='g33_cov',Xtitle='u3',/ylog,xstyle=1
plot,u3,g13(Lidx,*)/dmudx(Lidx,*),Title='g13_cov',Xtitle='u3',xstyle=1


plot,u3,Jac(Lidx,*)/dmudx(Lidx,*),Title='Jacobian',Xtitle='u3',/ylog,xstyle=1,ystyle=1
plot,u3,r_arr(Lidx,*),Title='R array',Xtitle='u3',xstyle=1
plot,u3,Acos(costh(Lidx,*))*180/!pi,Title='theta array',Xtitle='u3',xstyle=1

plot,x_arr,y_arr,Title='Grid',psym=3,xstyle=1
oplot,x_arr(Lidx,*),y_arr(Lidx,*),color=50,thick=2

!p.multi=0

If Plot_PNG EQ 1 then $
Begin
 Filename = out_pth+'2D_Faligned_Metric at L '+String(Lval,format='(F4.1)')+'.PNG'
 image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
end

DelS 		= Reform(sqrt((Shift(x_arr(Num_u1-1,*),1)-Shift(x_arr(Num_u1-1,*),0))^2 + (Shift(y_arr(Num_u1-1,*),1)-Shift(Y_arr(Num_u1-1,*),0))^2))
DelS(0)			= 0.0
OuterLength 	= Total(DelS,/cumulative)

Print,'Grid Spacing dsi = ',dsi/1.0e3,' km'
print,'min spacing in u3', min(abs(d32*h3))*Re/1.0e3,' km'
;stop

end			; end of Gen_grid

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

; Newtons Method Procedure
; Solve for R using Newton's method

pro New_r,r0,u1,u3,costh02,RI,ans
 MaxN=100
 Tol=1.e-3
 err=1.0
 ii=0
 While (ii lt MaxN) AND (max(err) gt Tol) do begin
  fr=u3*u3*r0^4*costh02/RI^3 - u1*r0 - RI
  df_dr=4.0*u3*u3*r0^3*costh02/RI^3 - u1
  ans=r0-fr/df_dr
  err=abs(ans-r0)
  r0=ans
 end
end

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


Pro Get_BasisFn, Num_u1, m, u1, del_u1, modefrac, r_iono, psiib3_N, psiib3_S, plot_png, out_pth, $
				 hratm_N, hratm_S, psigb3_N, psigb3_S, ixbc_0, ixbc_n

; set up basis functions for Atmospheric Solution
Scaleup	= 4						; Solve eigen problem on dense grid "scaleup" bigger that model grid
NBsets 	= 1						; Number of basis Sets to use

;  nmodes	= floor(size(evals,/n_elements)*modefrac)   ; select max spatial scale here
nmodes	= floor(Num_u1*modefrac)   ; select max spatial scale here
nm		= nmodes-1

evecs	= dcomplexarr(Num_u1,Num_u1-2)
ev_tmp = dcomplexarr(Num_u1,Nmodes,NBsets)
evl_tmp = dcomplexarr(Nmodes,NBsets)

For sets = 0,NBsets-1 do begin

If sets EQ 0 then ixbc_0 	= 1						; Low L Boundary Condition switch for basis functions ixbc eq 1 means derivative =0 BC else zero
If sets EQ 0 then ixbc_n 	= 1						; High L Boundary Condition switch for basis functions ixbc eq 1 means derivative =0 BC else zero

If sets EQ 1 then ixbc_0 	= 0						; Low L Boundary Condition switch for basis functions ixbc eq 1 means derivative =0 BC else zero
If sets EQ 1 then ixbc_n 	= 0						; High L Boundary Condition switch for basis functions ixbc eq 1 means derivative =0 BC else zero

Npts 	= Scaleup*(Num_u1-1)+1
IPts	= Npts-2

du1_sc 	= (u1(Num_u1-1)-u1(0))/(Npts-1)
u1_sc	= findgen(Npts)*du1_sc + u1(0)
du1_2 	= 2.0d0*du1_sc
du1_sq 	= du1_sc^2
indx 	= indgen(Npts/Scaleup+1)*Scaleup

a 		= replicate(0.0d0,Ipts,Ipts)  ; finite diff matrix
d2co 	= 4.0d0*u1_sc*(1+u1_sc)/du1_sq
d1co	= (2.0d0+3.0d0*u1_sc)/du1_sc
d0co	= m*m/u1_sc

;populate inner coeffs into matrix a
for i=2,Ipts-1 do begin
  a[i-2,i-1] 	=  d2co[i]-d1co[i]
  a[i-1,i-1] 	= -(d0co[i]+2*d2co[i])
  a[i,i-1] 		= d2co[i]+d1co[i]
end

if (ixbc_0 eq 1) then begin                ; ixbc eq 1 means derivative =0 BC else zero
  p0=d2co[1]-d1co[1]
endif else begin
  p0=0.0
endelse

if (ixbc_n eq 1) then begin                ; ixbc eq 1 means derivative =0 BC else zero
  pn=d2co[Npts-2]+d1co[Npts-2]
endif else begin
  pn=0.0
endelse

a[0,0] 				= -2.0d0*d2co[1]-d0co[1]+4.d0*p0/3.d0
a[1,0] 				= d2co[1]+d1co[1]-p0/3.d0
a[Npts-3,Npts-3] 	= -2.d0*d2co[Npts-2]-d0co[Npts-2]+4.d0*pn/3.d0
a[Npts-4,Npts-3]	= d2co[Npts-2]-d1co[Npts-2]-pn/3.d0

eval_sc				= La_eigenproblem(a,eigenvectors=evec1_sc)
evec_sc				= dcomplexarr(Npts,Npts-2)
evec_sc[1:Npts-2,*]	= evec1_sc

if(ixbc_0 eq 1) then begin
  evec_sc[0,*]		= (4.d0*evec_sc[1,*]     -evec_sc[2,*]     )/3.d0
endif else begin
  evec_sc[0,*]		= 0.0
endelse

if(ixbc_n eq 1) then begin
  evec_sc[Npts-1,*]	= (4.d0*evec_sc[Npts-2,*]-evec_sc[Npts-3,*])/3.d0
endif else begin
  evec_sc[Npts-1,*]	= 0.0
endelse

; sort eigenvectors
isort 	= sort(eval_sc)
evals 	= eval_sc[isort]

for ll 	= 0,Num_u1-3 do evecs(*,ll) = evec_sc[indx,isort(ll)]			; Sort and resample to smaller model grid

ev_tmp[*,0:nm,sets]		= evecs[*,0:nm]
evl_tmp[0:nm,sets]		= evals[0:nm]

evec_sc2       = complexarr(Npts,Npts-2)
openr,1,'~/Code/mhd2d/evecs.dat',/F77_UNFORMATTED
readu,1,evec_sc2
close,1

end				; end basis set loop.

If NBsets GT 1 then Ev  = [[ev_tmp(*,*,0)],[ev_tmp(*,*,1)]] else Ev = reform(ev_tmp(*,*,0))				; Concatinate Basis sets
If NBsets GT 1 then Evl =  [evl_tmp(*,0),evl_tmp(*,1)] else Evl = reform(evl_tmp(*,0))						; Concatinate Eigenvalues sets


;normalize
dufac 			= del_u1/(2.d0*sqrt(1.d0+u1))
dufac[0] 		= dufac[0]*0.5
dufac[num_u1-1] = dufac[num_u1-1]*0.5

for n=0,N_elements(ev[0,*])-1 do begin
  ev2 	= ev[*,n]*conj(ev[*,n])*dufac
  normfac = sqrt(total(ev2))
  ev[*,n] = ev[*,n]/normfac
endfor


nulm		= (sqrt(1.d0 + 4.d0*evl)-1.d0)*0.5
evt			= transpose(ev)
eve			= evt#ev
eveinv		= La_invert(eve)

;	For fitting only at Br locations
Eidx = Indgen(num_u1/2+1)*2					; Array of index location of even grid points
Oidx = Indgen(num_u1/2)*2+1					; Array of index location of odd grid points

b3coef		= eveinv # evt(*,Eidx)
;b3coef		= eveinv # evt(*,*)
b3b3		= ev # b3coef

;psiicoef 	= dcomplexarr(NBsets*nmodes,num_u1)
;psigcoef	= dcomplexarr(NBsets*nmodes,num_u1)

nufac0 		= nulm*r_iono*(1.d0 - 1.d0/r_iono^(2*nulm+1))
nufaci 		= (1.d0 + nulm/(nulm+1)/r_iono^(2*nulm+1))/nufac0
nufacg		= (2*nulm+1)/(nulm+1)/r_iono^nulm/nufac0

;	For Southern Ionosphere/Atmopshere
psifaci 	= nufaci#hratm_S(Eidx)
psifacg 	= nufacg#hratm_S(Eidx)
psiicoef	= b3coef*psifaci
psigcoef	= b3coef*psifacg
psiib3_S 	= ev # psiicoef   			;	now b3b3 should be nxp1 x nxp1 unit matrix
psigb3_S 	= ev # psigcoef

;	For Northern Ionosphere/Atmopshere
psifaci		= nufaci#hratm_N(Eidx)
psifacg		= nufacg#hratm_N(Eidx)
psiicoef	= b3coef*psifaci
psigcoef	= b3coef*psifacg
psiib3_N	= ev # psiicoef   			;	now b3b3 should be nxp1 x nxp1 unit matrix
psigb3_N	= ev # psigcoef

;set_plot,'win'
device,decomposed=0
;loadct,4
!p.multi=[0,1,2]
window,4,xsize=800,ysize=700,title='Basis Functions'
!P.charsize=1.6
!P.Charthick=1.2
!p.symsize=0.2

plot,u1,ev(*,0),title='Basis Functions set 0',xtitle='u1',ytitle='Amplitude',$
	 xstyle=1,ystyle=1,Yrange=[-max(ABS(ev)),max(ABS(ev))]
Oplot,u1,ev(*,0),color=1,psym=2
For ii = 0,9 do Oplot,u1,ev(*,ii),color=ii*20+1

If NBsets GT 1 then begin
plot,u1,ev(*,Nmodes),title='Basis Functions set 1',xtitle='u1',ytitle='Amplitude',$
	 xstyle=1,ystyle=1,Yrange=[-max(ABS(ev)),max(ABS(ev))]
Oplot,u1,ev(*,Nmodes),color=1,psym=2
For ii = 0,9 do Oplot,u1,ev(*,Nmodes+ii),color=ii*20+1
end

!p.multi=0

If Plot_png EQ 1 then $
begin
 Filename 	= out_pth+'Basis Functions.png'
 image 		= TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
end


; stop

End 		; End of  Basis Function Routine

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


Pro Get_Va,  	c_u, e0_u, u0_u, x_arr, y_arr, r_arr, z_arr, costh, Re, Num_u1, num_u3, plot_png, out_pth, Va_arr, eps_arr, rho_a, $
				cour, d32,  del_u1, h1, h2, h3, m, nt, tmax, dt, va2_arr, No_tot, $
		 		Neutral_file, num_u3_half, z0, mhocgs, hphiatm_N, hphiatm_S, hthatm_N, hthatm_S, bsqrt, $
				eta_arr, sig0_arr, sigp_arr, sigh_arr, sigpatm_N, sighatm_N, sig0atm_N, sigpatm_S, sighatm_S, sig0atm_S, $
				e1b1atm_N, e1b2atm_N, e2b1atm_N, e2b2atm_N, e1b1atm_S, e1b2atm_S, e2b1atm_S, e2b2atm_S


rho_a	= Dblarr(Num_u1,Num_u3)
flr_a	= Dblarr(Num_u1,8)        ; calc 7 harmonics
Va	 	= Dblarr(Num_u1,Num_u3)

No_E 	= dblarr(num_u1,num_u3)
No_F1 	= dblarr(num_u1,num_u3)
No_F2 	= dblarr(num_u1,num_u3)
No_Mag	= dblarr(num_u1,num_u3)



Eregion = 100e3								; 	E region Height in m
Mol_H	= 1.0e-3							;	Molecular Mass of Hydrogen (in SI)
Av_n	= 6.023e23							;	Avogardros Number
Mp		= 1.67e-27							;	Mass of Proton (kg)
Me		= 9.109e-31							;	Mass of electron (kg)
q_i		= 1.602e-19							;	Charge (Coulomb)
q_e		=-1.602e-19							;	Charge (Coulomb)


;	Magnetic field Strength
LVal	= Dblarr(Num_u1,Num_u3)
for k=0,num_u3-1 do LVal(*,k) = x_arr(*,num_u3/2)
K_0		= 8.0e15              ; Earth magnetic moment
Lat		= !pi/2 - Acos(costh)
B0		= sqrt(4.0-3.0*cos(Lat)^2)/(cos(Lat)^6) * K_0/(LVal*(Re))^3

;	using Colins FLR Va functions
xd 		= LVal(*,Num_u3/2)
;tanhm	= 2.0*(xd-4.5)									; Plasmapause at L = 4.5 (Used for 1st 2D Model Paper)
;tterm	= (55.0/xd^2.6)*(tanh(tanhm)+1.0)/2.0

tanhm=4.0*(xd-5.6)								;	Plasmapause at L = 5.6
tterm=(157./xd^2.6)*(tanh(tanhm)+1.0)/1.6

va_sp	= (2.09/xd^1.1+tterm)*990.0-50.0				; Va speed at equator
va_sp	= va_sp*1000.0


B_eq	= B0(*,num_u3/2)
nH		= B_eq^2/(u0_u*va_sp^2)						; Equatorial plasma density
nH_eq	= Dblarr(Num_u1,Num_u3)
for k=0,num_u3-1 do nH_eq(*,k) = nH

;nO		= 0.0e7										; Oxygen density (O2 Kg/Re^3)
;nO_Ht	= 250e3/Re									; Scale Height of Oxygen

Plaw	= 4.0										; Power law of density profile
rho_mag	= nH_eq*(Lval/r_arr)^Plaw ;+ nO/(Re)^3*Exp(-(r_arr - 1.0)/nO_Ht)

No_Mag	= rho_mag/mp								; Number of electrons (and protons) in magnetosphere from Power Law FLR profile

; Find_flr,xd,Num_u3,rho_a(ii,*),har
; flr_a(ii,0:7)=har(0:7)             				; Toroidal FLR frequency for this Va

;	Analytic profile near Ionosphere using Chapman layers (see Lysak '99)

For ii = 0, Num_u1-1 do begin
For jj = 0, Num_u3-1 do begin

;	E region
Ne_E 	= 1.5e11			;	Number of electron at E region peek (per m^3)
Alt_E	= Eregion			;	Altitude of e region max (in m)
Esht_t	= 300.0e3			;	Scale Height above E region (in m)
Esht_b	= 20.0e3			;	Scale Height below E region (in m)

If z_arr(ii,jj) GT Alt_E then X = (z_arr(ii,jj)-Alt_E)/Esht_t
If z_arr(ii,jj) LE Alt_E then X = (z_arr(ii,jj)-Alt_E)/Esht_b
No_E(ii,jj)	= Ne_E*Exp(1.0-X-Exp(-X))				; Using a Chapman profile

;	F1 region
Ne_F1 	= 2.5e11			;	Number of electron at F1 peek (per m^3)
Alt_F1	= 200.0e3			;	Altitude of F1 region (in m)
F1sht_t	= 200.0e3			;	Scale Height above F1 region (in m)
F1sht_b	= 40.0e3			;	Scale Height below F1 region (in m)

If z_arr(ii,jj) GT Alt_F1 then X = (z_arr(ii,jj)-Alt_F1)/F1sht_t
If z_arr(ii,jj) LE Alt_F1 then X = (z_arr(ii,jj)-Alt_F1)/F1sht_b
No_F1(ii,jj)= Ne_F1*Exp(1.0-X-Exp(-X))				; Using a Chapman profile

;	F2 region
Ne_F2 	= 2.0e12			;	Number of electron at F2 peek (per m^3)
Alt_F2	= 350.0e3			;	Altitude of F2 region (in km)
F2sht_t	= 175.0e3			;	Scale Height above F2 region (in m)
F2sht_b	= 75.0e3			;	Scale Height below F2 region (in m)

If z_arr(ii,jj) GT Alt_F2 then X = (z_arr(ii,jj)-Alt_F2)/F2sht_t
If z_arr(ii,jj) LE Alt_F2 then X = (z_arr(ii,jj)-Alt_F2)/F2sht_b
No_F2(ii,jj)= Ne_F2*Exp(1.0-X-Exp(-X))				; Using a Chapman profile

end
end

No_tot 	= No_E + No_F1 + No_F2 + No_Mag			; Total number of electrons (and Protons)

; Total Mass Density (assuming a Hydrogen Plasma)
Rho_a	= (No_tot/Av_N)*Mol_H				;   Mass Density in Kg/m^3


; Read ionospheric data from file

; Define output arrays
;eps_arr		=dblarr(num_u1,num_u3)
sigp_arr	=dblarr(num_u1,num_u3)
sigh_arr	=dblarr(num_u1,num_u3)
sig0_arr	=dblarr(num_u1,num_u3)
eta_arr		=dblarr(num_u1,num_u3)

sigpatm_N 	=Dcomplexarr(num_u1)
sighatm_N 	=Dcomplexarr(num_u1)
sig0atm_N 	=Dcomplexarr(num_u1)

sigpatm_S 	=Dcomplexarr(num_u1)
sighatm_S 	=Dcomplexarr(num_u1)
sig0atm_S 	=Dcomplexarr(num_u1)

No_N 		= Dblarr(Num_u1,Num_u3) 	; Number of Neutrals
Mav 		= Dblarr(Num_u1,Num_u3) 	; Average Molecular mass of Neutrals
Temp		= Dblarr(Num_u1,Num_u3)		; Temperature of Neutrals


;	Magnetic field Strength
LVal	= Dblarr(Num_u1,Num_u3)
for k=0,num_u3-1 do LVal(*,k) = x_arr(*,num_u3/2)
K_0		= 8.0e15              ; Earth magnetic moment
Lat		= !pi/2 - Acos(costh)
B0		= sqrt(4.0-3.0*cos(Lat)^2)/(cos(Lat)^6) * K_0/(LVal*(Re))^3

;	Gyrofrequency
Omega_p		=  q_i*B0/(Mp)						;	Gyrofrequency of Protons (Hydrogen in SI)
Omega_e		=  q_e*B0/(Me)						;	Gyrofrequency of Electron (in SI)

;	Neutral Atmosphere for Collison Frequency (from Kelley)
st1	= ''
y1 	= 0.0 & y2 	= 0.0 & y3 	= 0.0 & y4 	= 0.0e0

Pt			= 29
Neutral		= Dblarr(Pt+6,4)

OpenR,u,Neutral_file,/Get_Lun
Print,'Reading Data from ',Neutral_file
ReadF,u,st1								;	Header
ReadF,u,st1								;	Header
ReadF,u,st1								;	Header

;	Order in Neutral Atmopshere File:  Height[km], Temperature [k], Average Mass [g], No of Neutrals [cc]
;	SCale Height above reported region in Kelley [in km]
ScleHt = 2637.

For rr=1,Pt do begin $
 ReadF,u,y1,y2,y3,y4
Neutral(rr,0) = y1 & Neutral(rr,1)= y2 & Neutral(rr,2) = y3 & Neutral(rr,3) = y4
end
Neutral(0  ,0) = 0.     & Neutral(0,1) = Neutral(1  ,1)  & Neutral(0,2) = Neutral(1,2) & Neutral(0,3) = Neutral(1,3)
Neutral(rr  ,0) = 3000.  & Neutral(rr  ,1)= y2 & Neutral(rr  ,2) = y3 & Neutral(rr  ,3) = y4*exp(-(Neutral(rr  ,0)-Neutral(Pt,0))/ScleHt)
Neutral(rr+1,0) = 4000.	 & Neutral(rr+1,1)= y2 & Neutral(rr+1,2) = y3 & Neutral(rr+1,3) = y4*exp(-(Neutral(rr+1,0)-Neutral(Pt,0))/ScleHt)
Neutral(rr+2,0) = 5000.	 & Neutral(rr+2,1)= y2 & Neutral(rr+2,2) = y3 & Neutral(rr+2,3) = y4*exp(-(Neutral(rr+2,0)-Neutral(Pt,0))/ScleHt)
Neutral(rr+3,0) = 6000.	 & Neutral(rr+3,1)= y2 & Neutral(rr+3,2) = y3 & Neutral(rr+3,3) = y4*exp(-(Neutral(rr+3,0)-Neutral(Pt,0))/ScleHt)
Neutral(rr+4,0) = 100000  & Neutral(rr+4,1)= y2 & Neutral(rr+4,2) = y3 & Neutral(rr+4,3) = y4*exp(-(Neutral(rr+4,0)-Neutral(Pt,0))/ScleHt)
Free_Lun,u

;	Place Neutral Values on Model Grid Points (Assumes Radial Profile)

For ii = 0,Num_u1-1 do begin

Fline_No_N 	= INTERPOL( Neutral(*,3), Neutral(*,0)*1.0e3, z_arr(ii,*),/QUADRATIC)
Fline_Mav 	= INTERPOL( Neutral(*,2), Neutral(*,0)*1.0e3, z_arr(ii,*),/QUADRATIC)
Fline_Temp 	= INTERPOL( Neutral(*,1), Neutral(*,0)*1.0e3, z_arr(ii,*),/QUADRATIC)

No_N(ii,*)	= Fline_No_N
Mav(ii,*)	= Fline_Mav
Temp(ii,*)	= Fline_Temp

end

;	Collision Frequencies

Ne_cc	= No_tot/1.0e6							; number of electrons (per cc)

;;	from Kelley (possible mistake)
;Nu_ei	= [34.0+4.18*ALog((Temp^3)/Ne_cc)]*Ne_cc*Temp^(-1.5)
;Nu_en	= 5.4e-10*No_N*sqrt(Temp)
;Nu_in	= 2.6e-9*(No_N+Ne_cc)/sqrt(Mav)

;	From my 1D code (Zhang and Cole)
Nu_ei	= [34.0+4.18*ALog((Temp^3))/Ne_cc]*Ne_cc*Temp^(-1.5)
Nu_en	= 5.4e-10*No_N*sqrt(Temp)
Nu_in	= 2.6e-9*(No_N+Ne_cc)*sqrt(Mav)

Nu_e	= Nu_en + Nu_ei
Nu_i	= Nu_in ;+ Nu_ei											; Leaving Nu_ei in Nu_i leads to small negatives in Hall but gives reasonable profile

Lambda	= sqrt(Me/(4.0*!dpi*1.0e-7*No_tot*q_e^2))					; from Lysak99 p10,019

;	Conductivities (from 1d codes)
;sig0_arr	=	 1.0/(4.0*!dpi*1.0e-7*(Nu_e)*Lambda^2)				; from Lysak99 p10,019
sig0_arr		=	  No_tot*(q_e^2/(Me*Nu_e) + q_i^2/(Mp*Nu_i))	; from 1D code

sigp_arr	=  (No_tot*q_i^2/(Mp))* (Nu_i)/((Nu_i)^2+omega_p^2) $
           	  +(No_tot*q_e^2/(Me))* (Nu_e)/((Nu_e)^2+omega_e^2)

sigh_arr	= -( (No_tot*q_i^2/(Mp))* (omega_p)/((Nu_i)^2+omega_p^2) $
		        +(No_tot*q_e^2/(Me))* (omega_e)/((Nu_e)^2+omega_e^2) )

eta_arr		= (1.0/(u0_u*sig0_arr))/Re^2			; Re^2 comes from Scaling? (it works!)

;	Conductances in Thin Sheet (assume expodential decay from last "known" value, integrated analytically)
Thick_Sd = z_arr(*,0)
Scale_Sd = 1.0/0.2e3 															; 1/ Direct scale Height (in m)
sig0atm_N = sig0_arr(*,0)       *((1.0-Exp(-Scale_sd*Thick_sd))/Scale_sd)
sig0atm_S = sig0_arr(*,Num_u3-1)*((1.0-Exp(-Scale_sd*Thick_sd))/Scale_sd)

Thick_Sp = z_arr(*,0)
Scale_Sp = 1.0/0.2e3 															; 1/ Pederson scale Height (in m)
sigpatm_N = sigp_arr(*,0)       *((1.0-Exp(-Scale_sp*Thick_sp))/Scale_sp)
sigpatm_S = sigp_arr(*,Num_u3-1)*((1.0-Exp(-Scale_sp*Thick_sp))/Scale_sp)

Thick_Sh = z_arr(*,0)
Scale_Sh = 1.0/0.2e3 															; 1/ Hall scale Height (in m)
sighatm_N = sigh_arr(*,0)       *((1.0-Exp(-Scale_sh*Thick_sh))/Scale_sh)
sighatm_S = sigh_arr(*,Num_u3-1)*((1.0-Exp(-Scale_sh*Thick_sh))/Scale_sh)

Tol = 0.0e-0
Zero = where (sigh_arr LE Tol)
sigh_arr(Zero) = Tol

Tol = min(eta_arr)
;zero = where (eta_arr(0.3*Num_u1:0.7*Num_u1,*) GE Tol)+0.3*(Num_u1-1)
;eta_arr(*,0.3*Num_u3:0.7*Num_u3) = Tol

; stop


;;;	Atmospheric thin sheet conductances
;sigpatm_N(*) = 1.05
;sigpatm_S(*) = 1.05
;sighatm_N(*) = 1.05
;sighatm_S(*) = 1.05
;sig0atm_N(*) = 500.0
;sig0atm_S(*) = 500.0



;;	Height Distributed Conductivity Arrays (set to zero if thin sheet is only required)
;sigp_arr[*,*]= 0.0
;sigh_arr[*,*]= 0.0

;sigh_arr[0,*]		= 0.0		; Sets inner L shells Sigma_h to zero -> no mode conversion
;sigh_arr[Num_u1-1,*]= 0.0		; Sets inner L shells Sigma_h to zero -> no mode conversion


;;	Set Corner Conductance in thin sheet to ?
;sigpatm_N(0) = 0.0
;sighatm_N(0) = 0.0

;sigpatm_S(0) = 0.0
;sighatm_S(0) = 0.0


; Calculate Va and Wave Speed
; ~~~~~~~~~~~~~~~~~~~~~~~~~~~

V_arr		= B0/sqrt(u0_u*rho_a)						; Alfven Speed in m/s ()
Coll_Mod	= omega_p^2/(omega_p^2 + Nu_i^2)			; see Lysak '99 eq 7: Note ignore electron term since me << Mav
Eperp		= e0_u*(1+(c_u^2/V_arr^2)*Coll_mod)			; E perp in SI (including collision modification)
eps_arr		= Eperp;/Re
Wave_sp		= c_u/sqrt(Eperp/e0_u)						; The actually wave speed including collisions modification

Va_arr		= Wave_sp/Re
;Va_arr		= V_arr/Re									; No Modification by collisions to wave speed

va2_arr 	= va_arr*va_arr



; plot,z_arr(100,0:50),Nu_i(100,0:50), /ylog,/xlog
;oplot,z_arr(100,0:50),Omega_p(100,0:50),color=150
 ; plot,z_arr(100,0:50),Va_arr(100,0:50)*Re,/ylog,/xlog
;oplot,z_arr(100,0:50),V_arr(100,0:50),color=150


; determine time step
dtz=cour*min(abs(d32*h3)/Va_arr)
dtx=cour*min(abs(2.0*del_u1*h1)/Va_arr)
if(m ne 0) then begin                        ; we have azi variation
  dty=cour*min(h2/m/Va_arr)
  dt=1.0/sqrt(3.0)/(sqrt(1.0/dtx^2+1.0/dty^2+1.0/dtz^2))
  print,'dtx,y,z=',dtx,dty,dtz
;    printf,u_o,'dtx,y,z=',dtx,dty,dtz
endif else begin
  dt=1.0/sqrt(2.0)/sqrt(1.0/dtx^2+1.0/dtz^2)
  print,'dtx,z=',dtx,dtz
;    printf,u_o,'dtx,z=',dtx,dtz
endelse

nt=tmax/dt

dtzmax=cour*max(abs(d32*h3)/Va_arr)
dtxmax=cour*max(abs(2.0*del_u1*h1)/Va_arr)


print,'min spacing in u1', min(abs(2.0*del_u1*h1))*Re,' m'
print,'min spacing in u3', min(abs(d32*h3))*Re,' m'
print,'max Va', max(abs(Va_arr))*Re,' m/s'

print,'* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
print,'Total simulation time, tmax: ',tmax
print,'Time step: ',dt
print,'so number of time loops = ', nt

courx	= Va_arr*dt/abs(del_u1*h1)
cmax	= max(courx,idim)
i1z	= floor(idim/(num_u1))
i1x	= idim-i1z*(num_u1)
print,'Max Cour_len u1: ',cmax,' at i1x,i1z = ',i1x,i1z

courz	= Va_arr*dt/abs(h3)
cmaz	= max(courz,idim)
i3z	= floor(idim/(num_u1))
i3x	= idim-i3z*(num_u1)
print,'Max Cour_len u3: ',cmaz,' at i3x,i3z = ',i3x,i3z

print,'num_u1 x nzp1 x nt=',(num_u1)*(num_u3)*nt
print,'Max travel time between grid points in u1, u3 directions = ',dtxmax,dtzmax
print,'* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'

flr_a=dblarr(Num_u1,8)        ; calc 6 harmonics

For ii=0,Num_u1-1 do $        ; chooses field line
Begin
xd=X_arr(ii,Num_u3/2)
Find_flr,xd,Num_u3,rho_a(ii,*),har
flr_a(ii,0:7)=har(0:7)             ; Fundamental FLR
end



;	Factors involving conductances in the numerics of thin sheet current sheet
;	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

; 	musigfac gives muSig in s/Re for 1 mho
musigfac 	= u0_u*Re

;	For Northern Hemipshere
cosalpha_N 	= -2.0*costh[*,0]/bsqrt[*,0]				; set as postive for now, costh0 only goes 0< x < 90 deg
sin2alpha_N	= 1.0d - cosalpha_N*cosalpha_N

musigp 	= sigpatm_N*musigfac
musigh 	= sighatm_N*musigfac
musig0 	= sig0atm_N*musigfac

musigzz	= musig0 + (musigp-musig0)*sin2alpha_N
sigtt 	= musig0*musigp/musigzz
sigtp 	=-musig0*musigh*cosalpha_N/musigzz
sigpt 	= musig0*musigh*cosalpha_N/musigzz
sigpp	= musigp + musigh*musigh*sin2alpha_N/musigzz

sig11 	= sigtt*hphiatm_N/hthatm_N
sig12 	= sigtp
sig21 	= sigpt
sig22 	= sigpp*hthatm_N/hphiatm_N
sigden 	= sig11*sig22-sig12*sig21

e1b1atm_N =-sig12/sigden
e1b2atm_N =-sig22/sigden
e2b1atm_N = sig11/sigden
e2b2atm_N = sig21/sigden

;	For southern Hemipshere
cosalpha_S	= -2.0*costh[*,Num_u3-1]/bsqrt[*,Num_u3-1]
sin2alpha_S	= 1.0d - cosalpha_S*cosalpha_S

musigp 	= sigpatm_S*musigfac
musigh 	= sighatm_S*musigfac
musig0 	= sig0atm_S*musigfac

musigzz	= musig0 + (musigp-musig0)*sin2alpha_S
sigtt 	= musig0*musigp/musigzz
sigtp 	=-musig0*musigh*cosalpha_S/musigzz
sigpt	= musig0*musigh*cosalpha_S/musigzz
sigpp	= musigp + musigh*musigh*sin2alpha_S/musigzz

sig11	= sigtt*hphiatm_S/hthatm_S
sig12 	= sigtp
sig21 	= sigpt
sig22 	= sigpp*hthatm_S/hphiatm_S
sigden 	= sig11*sig22-sig12*sig21

e1b1atm_S =-sig12/sigden
e1b2atm_S =-sig22/sigden
e2b1atm_S = sig11/sigden
e2b2atm_S = sig21/sigden



device,decomposed=0
!p.multi=0
window,5,xsize=800,ysize=800,title='Va Profile'
!P.charsize=1.6
!P.Charthick=1.2
!p.symsize=0.2

!p.multi=[0,2,3]

Lsh1 	= 9.0
Lidx1 	= min(where (r_arr(*,num_u3/2) GE Lsh1))
Lsh2 	= 6.0
Lidx2 	= min(where (r_arr(*,num_u3/2) GE Lsh2))
Lsh3 	= 3.0
Lidx3 	= min(where (r_arr(*,num_u3/2) GE Lsh3))

contour,Alog10(va_arr),/fill,nlevels=255, xstyle=1,ystyle=1,$
 title='Log!d10!n (V) (modified by colliisons)',xtitle='u1 index',ytitle='u3 index'
Plots,i1x,i1z,color=0,thick=2,psym=2,symsize=3
Plots,i3x,i3z,color=250,thick=2,psym=2,symsize=3

contour,Alog10(va_arr),x_arr,y_arr,/fill,nlevels=255,/isotropic, xstyle=1,ystyle=1,$
 title='Log!d10!n (V) (modified by colliisons)',xtitle='Re',ytitle='Re'
Oplot,x_arr(Lidx1,*),y_arr(Lidx1,*),color=0,thick=1
Oplot,x_arr(Lidx2,*),y_arr(Lidx2,*),color=150,thick=1
Oplot,x_arr(Lidx3,*),y_arr(Lidx3,*),color=200,thick=1
Plots,x_arr(i1x,i1z),y_arr(i1x,i1z),color=0,thick=2,psym=2,symsize=3
Plots,x_arr(i3x,i3z),y_arr(i3x,i3z),color=250,thick=2,psym=2,symsize=3

contour,Alog10(rho_a),/fill,nlevels=255, xstyle=1,ystyle=1, $
 title='Log rho',xtitle='u1 index',ytitle='u3 index'
contour,Alog10(rho_a),x_arr,y_arr,/fill,nlevels=255,/isotropic, xstyle=1,ystyle=1, $
 title='Log rho',xtitle='Re',ytitle='Re'
Oplot,x_arr(Lidx1,*),y_arr(Lidx1,*),color=0,thick=1
Oplot,x_arr(Lidx2,*),y_arr(Lidx2,*),color=150,thick=1
Oplot,x_arr(Lidx3,*),y_arr(Lidx3,*),color=200,thick=1

contour,Eperp,/fill,nlevels=255, xstyle=1,ystyle=1, $
 title='E perp',xtitle='u1 index',ytitle='u3 index'
;contour,Eperp,x_arr,y_arr,/fill,nlevels=255,/isotropic, xstyle=1,ystyle=1, $
; title='E perp',xtitle='Re',ytitle='Re'
;Oplot,x_arr(Lidx1,*),y_arr(Lidx1,*),color=0,thick=1
;Oplot,x_arr(Lidx2,*),y_arr(Lidx2,*),color=150,thick=1
;Oplot,x_arr(Lidx3,*),y_arr(Lidx3,*),color=200,thick=1

Plot,X_arr(*,Num_u3/2),Flr_a(*,0),Xrange=[0,10], Yrange=[0.1,1000],Ytype=1,$
 Title = 'Harmonic Structure (mHz) ',Xtitle='L Shell (in Re)',Ytitle='Frequency (mHz)'
Oplot,X_arr(*,Num_u3/2),Flr_a(*,1),linestyle=1
Oplot,X_arr(*,Num_u3/2),Flr_a(*,2),linestyle=2
Oplot,X_arr(*,Num_u3/2),Flr_a(*,3),linestyle=1
Oplot,X_arr(*,Num_u3/2),Flr_a(*,4),linestyle=2
Oplot,X_arr(*,Num_u3/2),Flr_a(*,5),linestyle=1
Oplot,X_arr(*,Num_u3/2),Flr_a(*,6),linestyle=2


!p.multi=0

If Plot_png EQ 1 then $
begin
 Filename 	= out_pth+'Va Profile.png'
 image 		= TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
end


device,decomposed=0
!p.multi=0
window,15,xsize=1200,ysize=600,title='Electron Profile'
!P.charsize=1.2
!P.Charthick=1.2
!p.symsize=0.2

!p.multi=[0,2,1]

;Lval = 6
;Lpt  = max(where(r_arr(*,Num_u3/2) LE Lval))
;Zval = 10000e3
;Zpt  = max(where(z_arr(Lpt,0:Num_u3/2) LE Zval))

 plot,No_tot(Lidx1,*),z_arr(Lidx1,*)/Re,/xlog,/ylog,$
 Title='Electron Profile (NH) along L='+String(Lsh1,format='(F3.1)')+','+String(Lsh2,format='(F3.1)')+','+String(Lsh3,format='(F3.1)')+ ' field line',$
 Xtitle='N!de!n [m!u-3!n]', Ytitle='Radial Distance (Height) [R!dE!n]'
oplot,No_E(Lidx1,*),z_arr(Lidx1,*)/Re,linestyle=1
oplot,No_F1(Lidx1,*),z_arr(Lidx1,*)/Re,linestyle=2
oplot,No_F2(Lidx1,*),z_arr(Lidx1,*)/Re,linestyle=3
oplot,No_Mag(Lidx1,*),z_arr(Lidx1,*)/Re,linestyle=4

oplot,No_tot(Lidx2,*),z_arr(Lidx2,*)/Re,color=150
oplot,No_E(Lidx2,*),z_arr(Lidx2,*)/Re,linestyle=1,color=150
oplot,No_F1(Lidx2,*),z_arr(Lidx2,*)/Re,linestyle=2,color=150
oplot,No_F2(Lidx2,*),z_arr(Lidx2,*)/Re,linestyle=3,color=150
oplot,No_Mag(Lidx2,*),z_arr(Lidx2,*)/Re,linestyle=4,color=150

oplot,No_tot(Lidx3,*),z_arr(Lidx3,*)/Re,color=200
oplot,No_E(Lidx3,*),z_arr(Lidx3,*)/Re,linestyle=1,color=200
oplot,No_F1(Lidx3,*),z_arr(Lidx3,*)/Re,linestyle=2,color=200
oplot,No_F2(Lidx3,*),z_arr(Lidx3,*)/Re,linestyle=3,color=200
oplot,No_Mag(Lidx3,*),z_arr(Lidx3,*)/Re,linestyle=4,color=200

 plot,Va_arr(Lidx1,*)*Re/1.0e3,z_arr(Lidx1,*)/Re,/xlog,/ylog,$
 Title='Wave Speed along L='+String(Lsh1,format='(F3.1)')+','+String(Lsh2,format='(F3.1)')+','+String(Lsh3,format='(F3.1)')+ ' field line',$
  Xtitle='Wave Speed [kms!u-1!n]', Ytitle='Radial Distance (Height) [R!dE!n]'
oplot,Va_arr(Lidx2,*)*Re/1.0e3,z_arr(Lidx2,*)/Re,color=150
oplot,Va_arr(Lidx3,*)*Re/1.0e3,z_arr(Lidx3,*)/Re,color=200

!p.multi=0

If Plot_png EQ 1 then $
begin
 Filename 	= out_pth+'Electron Profile.png'
 image 		= TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
end



device,decompose=0
window,2,Title='Conductances',Xsize=1200,Ysize=900
!P.multi=[0,3,4]
!P.charsize=2.5
!P.Charthick=1.2

Lsh1 	= 9.0
Lidx1 	= min(where (r_arr(*,num_u3_half) GE Lsh1))
Lsh2 	= 6.0
Lidx2 	= min(where (r_arr(*,num_u3_half) GE Lsh2))
Lsh3 	= 3.0
Lidx3 	= min(where (r_arr(*,num_u3_half) GE Lsh3))

Contour,Alog10(sig0_arr),/fill,nlevels=255,xstyle=1,ystyle=1, $
 Title='log!d10!n !7r!3!d0!n',Xtitle='u1 max = '+String(max(sig0_arr),format='(e10.3)')+' min = '+String(min(sig0_arr),format='(e10.3)'), Ytitle='u3'

Contour,Alog10(sig0_arr),x_arr,y_arr,/fill,nlevels=255,/isotropic,xstyle=1,ystyle=1
Oplot,x_arr(Lidx1,*),y_arr(Lidx1,*),color=0,thick=1
Oplot,x_arr(Lidx2,*),y_arr(Lidx2,*),color=150,thick=1
Oplot,x_arr(Lidx3,*),y_arr(Lidx3,*),color=200,thick=1

 Plot,sig0atm_N, Title='!7R!3!d0!n',Xtitle='u1'
OPlot,sig0atm_S,color=150

Contour,Alog10(sigp_arr),/fill,nlevels=255,xstyle=1,ystyle=1, $
 Title='log!d10!n !7r!3!dP!n',Xtitle='u1 max = '+String(max(sigp_arr),format='(e10.3)')+' min = '+String(min(sigp_arr),format='(e10.3)'), Ytitle='u3'

Contour,Alog10(sigp_arr),x_arr,y_arr,/fill,nlevels=255,/isotropic,xstyle=1,ystyle=1
Oplot,x_arr(Lidx1,*),y_arr(Lidx1,*),color=0,thick=1
Oplot,x_arr(Lidx2,*),y_arr(Lidx2,*),color=150,thick=1
Oplot,x_arr(Lidx3,*),y_arr(Lidx3,*),color=200,thick=1

 Plot,sigpatm_N, Title='!7R!3!dP!n',Xtitle='u1'
OPlot,sigpatm_S,color=150

Contour,Alog10(ABS(sigh_arr)),/fill,nlevels=255,xstyle=1,ystyle=1, $
 Title='log!d10!n !7r!3!dH!n',Xtitle='u1 max = '+String(max(sigh_arr),format='(e10.3)')+' min = '+String(min(sigh_arr),format='(e10.3)'), Ytitle='u3'

Contour,Alog10(sigh_arr),x_arr,y_arr,/fill,nlevels=255,/isotropic,xstyle=1,ystyle=1
Oplot,x_arr(Lidx1,*),y_arr(Lidx1,*),color=0,thick=1
Oplot,x_arr(Lidx2,*),y_arr(Lidx2,*),color=150,thick=1
Oplot,x_arr(Lidx3,*),y_arr(Lidx3,*),color=200,thick=1

Plot,sighatm_N, Title='!7R!3!dH!n',Xtitle='u1'
OPlot,sighatm_N,color=150

Contour,Alog10(eta_arr),/fill,nlevels=255,xstyle=1,ystyle=1, $
 Title='log!d10!n !7g!3!d0!n',Xtitle='u1 max = '+String(max(eta_arr),format='(e10.3)')+' min = '+String(min(eta_arr),format='(e10.3)'), Ytitle='u3'

Contour,Alog10(eta_arr),x_arr,y_arr,/fill,nlevels=255,/isotropic,xstyle=1,ystyle=1
Oplot,x_arr(Lidx1,*),y_arr(Lidx1,*),color=0,thick=1
Oplot,x_arr(Lidx2,*),y_arr(Lidx2,*),color=150,thick=1
Oplot,x_arr(Lidx3,*),y_arr(Lidx3,*),color=200,thick=1

 Plot,sig0_arr(Lidx1,*),z_arr(Lidx1,*)/Re,/xlog,/ylog,$
 Xrange=[min(sigh_arr(Lidx1,*)),max(sig0_arr(Lidx1,*))],$
 Title='Along L='+String(Lsh1,format='(F4.1)')+','+String(Lsh2,format='(F4.1)')+','+String(Lsh3,format='(F4.1)')+' field Lines (NH)',Xtitle='Mho',Ytitle='R!dE!n'
oPlot,sigp_arr(Lidx1,*),z_arr(Lidx1,*)/Re,linestyle=1
oPlot,sigh_arr(Lidx1,*),z_arr(Lidx1,*)/Re,linestyle=2

oPlot,sig0_arr(Lidx2,*),z_arr(Lidx2,*)/Re,color=150
oPlot,sigp_arr(Lidx2,*),z_arr(Lidx2,*)/Re,linestyle=1,color=150
oPlot,sigh_arr(Lidx2,*),z_arr(Lidx2,*)/Re,linestyle=2,color=150

oPlot,sig0_arr(Lidx3,*),z_arr(Lidx3,*)/Re,color=200
oPlot,sigp_arr(Lidx3,*),z_arr(Lidx3,*)/Re,linestyle=1,color=200
oPlot,sigh_arr(Lidx3,*),z_arr(Lidx3,*)/Re,linestyle=2,color=200


;plot,No_N(Lidx1,*),z_arr(Lidx1,*)/Re,/xlog,/ylog


;Plot,Sig0atm_N, Title='!7R!3!d0!n',Xtitle='u1'
;OPlot,Sig0atm_S,color=150

If Plot_PNG EQ 1 then $
Begin
 Filename = out_pth+'2D_Conductances.PNG'
 image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
end

;stop

End						; End of Va routine

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

; Routine to generate the numerical factors used in the model coming from the discretisation etc

Pro Get_Facts,	sigp_arr, sigh_arr, eps_arr, eta_arr, g11, g22, g33, g13, gsup11, gsup13, gsup33, h1, h2, h3, Jac, dt, $
				va_arr, va2_arr, im, d32, del_u1,$
				e1e1, e1e2, e1f1, e1f2, f1b2, f1b3, e2e1, e2e2, e2f1, e2f2, f2b1, f2b3, e3b1, e3b21, e3b23, e3b3, $
				b1e2, b1e3, b2e1, b2e3, b3e1, b3e2, e1esup1, e1e3, dx2, dxsq, z_arr, r_arr, Plot_PNG, out_pth


; factors for integration
  sigp_eps_arr 	= sigp_arr / eps_arr
  sigh_eps_arr 	= sigh_arr / eps_arr

  esig_arr 	= exp(-sigp_eps_arr*dt)
  esig2_arr = exp(-sigp_eps_arr*dt*0.5)
  csig_arr 	= cos(sigh_eps_arr*dt)
  csig2_arr = cos(sigh_eps_arr*dt*0.5)
  ssig_arr 	= sin(sigh_eps_arr*dt)
  ssig2_arr = sin(sigh_eps_arr*dt*0.5)

; factors for generator
  gg12	=g22*g33/h3/jac
  gg21	=(g11*g33-g13*g13)/h3/jac

  e1e1 	= esig_arr*csig_arr
  e1e2 	= gg12*esig_arr*ssig_arr

  e1f1 	= dt*esig2_arr*csig2_arr
  e1f2 	= dt*gg12*esig2_arr*ssig2_arr

  f1b2 	= -va2_arr/jac/d32
  f1b3 	=  va2_arr*im/jac

  e2e1 	= -gg21*esig_arr*ssig_arr
  e2e2 	= e1e1

  e2f1 	= -dt*gg21*esig2_arr*ssig2_arr
  e2f2 	=  dt*esig2_arr*csig2_arr

  f2b1 	= va2_arr/jac/d32
  f2b3 	=-va2_arr/jac/del_u1*0.5

  e3b1 	= -im*eta_arr*g33/jac
  e3b21	= eta_arr/jac*g33/del_u1*0.5
  e3b23	=-eta_arr/jac*g13/d32
  e3b3 	= im*eta_arr*g13/jac

  b1e2 	= dt/d32/jac
  b1e3 	=-dt*im/jac

  b2e1 	=-dt/d32/jac
  b2e3 	= dt/del_u1/jac*0.5

  b3e1 		= dt*im/jac
  b3e2 		=-dt/del_u1/jac*0.5

;  e1esup1 	= g11/(1-g13*gsup13)
;  e1e3 		= g13*gsup33/(1-g13*gsup13)

  e1esup1	=1.0/gsup11						; From Bob's fortran code (and word document)
  e1e3		=-gsup13/gsup11					; From Bob's fortran code (and word document)

  dx2 	= 2.0d0*del_u1
  dxsq 	= del_u1*del_u1






device,decompose=0
window,24,Title='Factors',Xsize=1500,Ysize=500
!P.multi=[0,3,1]
!P.charsize=2.5
!P.Charthick=1.2

u3_half = N_elements(r_arr(1,*))/2

Lsh1 	= 9.0
Lidx1 	= min(where (r_arr(*,u3_half) GE Lsh1))
Lsh2 	= 6.0
Lidx2 	= min(where (r_arr(*,u3_half) GE Lsh2))
Lsh3 	= 3.0
Lidx3 	= min(where (r_arr(*,u3_half) GE Lsh3))

H1idx 	= where(r_arr(Lidx1,0:u3_half)*6370 LE (1100+6370) )
H2idx 	= where(r_arr(Lidx2,0:u3_half)*6370 LE (1100+6370) )
H3idx 	= where(r_arr(Lidx3,0:u3_half)*6370 LE (1100+6370) )

 Plot,(sigh_eps_arr(Lsh1,H1idx)*dt),(r_arr(Lsh1,H1idx)-1.0)*6370,color=255, $
 xtitle='sigh_eps_arr * dt',ytitle='Altitude [km]',Title='',/ylog,Ystyle=1,Yrange=[50,1100]
OPlot,(sigh_eps_arr(Lsh1,H1idx)*dt),(r_arr(Lsh1,H1idx)-1.0)*6370,color=255,psym=2

OPlot,(sigh_eps_arr(Lsh2,H2idx)*dt),(r_arr(Lsh2,H2idx)-1.0)*6370,color=150
OPlot,(sigh_eps_arr(Lsh2,H2idx)*dt),(r_arr(Lsh2,H2idx)-1.0)*6370,color=150,psym=2

OPlot,(sigh_eps_arr(Lsh3,H3idx)*dt),(r_arr(Lsh3,H3idx)-1.0)*6370,color=200
OPlot,(sigh_eps_arr(Lsh3,H3idx)*dt),(r_arr(Lsh3,H3idx)-1.0)*6370,color=200,psym=2

 Plot,(sigh_arr(Lsh1,H1idx)),(r_arr(Lsh1,H1idx)-1.0)*6370,color=255, $
 xtitle='sigh_arr',ytitle='Altitude [km]',Title='',/ylog,Ystyle=1,Yrange=[50,1100]
oPlot,(sigh_arr(Lsh1,H1idx)),(r_arr(Lsh1,H1idx)-1.0)*6370,color=255,Psym=2

oPlot,(sigh_arr(Lsh2,H2idx)),(r_arr(Lsh2,H2idx)-1.0)*6370,color=150
oPlot,(sigh_arr(Lsh2,H2idx)),(r_arr(Lsh2,H2idx)-1.0)*6370,color=150,Psym=2

oPlot,(sigh_arr(Lsh3,H3idx)),(r_arr(Lsh3,H3idx)-1.0)*6370,color=200
oPlot,(sigh_arr(Lsh3,H3idx)),(r_arr(Lsh3,H3idx)-1.0)*6370,color=200,Psym=2

 Plot,(eps_arr(Lsh1,H1idx)),(r_arr(Lsh1,H1idx)-1.0)*6370,color=255, $
 xtitle='eps_arr',ytitle='Altitude [km]',Title='',/ylog,Ystyle=1,Yrange=[50,1100]
oPlot,(eps_arr(Lsh1,H1idx)),(r_arr(Lsh1,H1idx)-1.0)*6370,color=255,Psym=2

oPlot,(eps_arr(Lsh2,H2idx)),(r_arr(Lsh2,H2idx)-1.0)*6370,color=150
oPlot,(eps_arr(Lsh2,H2idx)),(r_arr(Lsh2,H2idx)-1.0)*6370,color=150,Psym=2

oPlot,(eps_arr(Lsh3,H3idx)),(r_arr(Lsh3,H3idx)-1.0)*6370,color=200
oPlot,(eps_arr(Lsh3,H3idx)),(r_arr(Lsh3,H3idx)-1.0)*6370,color=200,Psym=2



If Plot_PNG EQ 1 then $
Begin
 Filename = out_pth+'Factors.PNG'
 image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
end



stop
end 				; End of Factors routine

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 pro Iterate,  	num_u1,num_u3, num_u3_half, dt, va_arr, h1, h2, h3, h30, va2_arr, d32, del_u1, $
				g11,g22, g33, g13, u1, u3,gsup11,gsup22,gsup33, gsup13, out_pth, im, Phi, $
				e1e1, e1e2, e1f1, e1f2, f1b2, f1b3, e2e1, e2e2, e2f1, e2f2, f2b1, f2b3, e3b1, e3b21, e3b23, e3b3, $
				b1e2, b1e3, b2e1, b2e3, b3e1, b3e2, e1esup1, e1e3, dx2, dxsq, tmax, dtplot, $
				e1b2atm_N,e1b1atm_N, e2b1atm_N, e2b2atm_N, e1b2atm_S,e1b1atm_S, e2b1atm_S, e2b2atm_S, psiib3_N, psiib3_S, $
				psigb3_N, psigb3_S, drive0, Freq, Plot_fields, x_arr, y_arr, hratm_N, hratm_S ,$
				hphiatm_N, hphiatm_S, hthatm_N, hthatm_S,hthg_N, hthg_S, hphig_N, hphig_S, Re, ixbc_0, ixbc_n, OuterLength

;initialize arrays
e1 		= dcomplexarr(num_u1,num_u3)
e2 		= dcomplexarr(num_u1,num_u3)
e3 		= dcomplexarr(num_u1,num_u3)

b1 		= dcomplexarr(num_u1,num_u3)
b2 		= dcomplexarr(num_u1,num_u3)
b3 		= dcomplexarr(num_u1,num_u3)

esup1 	= dcomplexarr(num_u1,num_u3)
esup2 	= dcomplexarr(num_u1,num_u3)
bsup1 	= dcomplexarr(num_u1,num_u3)
bsup2 	= dcomplexarr(num_u1,num_u3)
bsup3 	= dcomplexarr(num_u1,num_u3)

cfac	=dt*va_arr*0.5/h3
nplots	=tmax/dtplot

Odds_u3 = Dcomplexarr(Num_u3)
Evens_u3= Dcomplexarr(Num_u3)
Odds_u1 = Dcomplexarr(Num_u1)
Evens_u1= Dcomplexarr(Num_u1)

For ii = 0,Num_u3/2-1 do $
begin
  Odds_u3(2*ii+1) 		= 1.0
  Evens_u3(2*ii) 		= 1.0
end

Evens_u3(Num_u3-1) = 1.0

For ii = 0,Num_u1/2-1 do $
begin
  Odds_u1(2*ii+1) 		= 1.0
  Evens_u1(2*ii) 		= 1.0
end

Evens_u1(0) = 0.0							; For odd number of grid points
Evens_u1(Num_u1-1) = 0.0					; For odd number of grid points


Eidx = Indgen(num_u1/2+1)*2					; Array of index location of even grid points
Oidx = Indgen(num_u1/2)*2+1					; Array of index location of odd grid points


i		= long(0)
t 		= 0.0d0
tpnext 	= dtplot
iplot 	= 0;
u3		= findgen(num_u3)


; ################################################
; Start time loop
; ################################################

Print,'Start of time Loop...'
while !d.window GT -1 do Wdelete 				; kill ALL open windows
Starttime = SYSTIME(/SECONDS)

while (iplot lt nplots) do begin
  t 	= t+dt									; Time

Get_Driver,   drive0, freq, t, h3, h30, u3, num_u3_half, num_u1, num_u3, Re, Driver, OuterLength

; shift ei's
    e1p3 	= shift(e1,0,-1)
    e1m3 	= shift(e1,0,1)
    e1p3[*,num_u3-1]= 0.0
    e1m3[*,0] 		= 0.0

    e2p1 	= shift(e2,-1,0)
    e2m1 	= shift(e2,1,0)
    e2p1[num_u1-1,*]= 0.0
    e2m1[0,*] 		= 0.0

    e2p3 	= shift(e2,0,-1)
    e2m3 	= shift(e2,0,1)
    e2p3[*,num_u3-1]= 0.0
    e2m3[*,0] 		= 0.0

    e3p1 	= shift(e3,-1,0)
    e3m1 	= shift(e3,1,0)
    e3p1[num_u1-1,*]= 0.0
    e3m1[0,*] 		= 0.0

; 	advance bsupi
; -----------------------------------------------------------
    bsup1 = bsup1 + b1e2*(e2p3-e2m3) + b1e3*e3
    bsup2 = bsup2 + b2e1*(e1p3-e1m3) + b2e3*(e3p1-e3m1)
    bsup3 = bsup3 + b3e1*e1 + b3e2*(e2p1-e2m1)

;if (ixbc_0 eq 0) then begin
; 	bsup3[0,*] 				= 0.d0
; 	bsup3[0,0] 				= 0.d0
;    bsup3[0,num_u3-1] 		= 0.d0
;endif else begin
;    bsup2[0,*] 				= (4.0*bsup2(2,*) - bsup2(4,*))/3.0    									; Inner L shell Perfectly reflecting
;    bsup3[0,*] 				= (4.0*bsup3(2,*) - bsup3(4,*))/3.0

    bsup3[0,0] 				= (4.0*bsup3(2,0) - bsup3(4,0))/3.0										; Corners
    bsup3[0,num_u3-1] 		= (4.0*bsup3(2,num_u3-1) - bsup3(4,num_u3-1))/3.0
;endelse
;
;if (ixbc_n eq 0) then begin
;    bsup3[num_u1-1,*] 		= 0.d0
;    bsup3[num_u1-1,0] 		= 0.d0
;    bsup3[num_u1-1,num_u3-1]= 0.d0
;endif else begin
;    bsup2[num_u1-1,*] 				= (4.0*bsup2(num_u1-3,*) - bsup2(num_u1-5,*))/3.0    			; Outer L shell Perfectly reflecting???????
;    bsup3[num_u1-1,*] 		= (4.0*bsup3(num_u1-3,*) - bsup3(num_u1-5,*))/3.0

    bsup3[num_u1-1,0] 		= (4.0*bsup3(num_u1-3,0) - bsup3(num_u1-5,0))/3.0						; Corners
    bsup3[num_u1-1,num_u3-1]= (4.0*bsup3(num_u1-3,num_u3-1) - bsup3(num_u1-5,num_u3-1))/3.0
;endelse


; average bsups
    Av_bsup1 		= (Shift(bsup1,1,1) + Shift(bsup1,1,-1) + Shift(bsup1,-1,1) + Shift(bsup1,-1,-1))/ 4.d0
    Av_bsup1(*,0)	= 0.0 		& 	Av_bsup1(*,Num_u3-1)	= 0.0			; Clean up after shifts due to odd number of points
    Av_bsup1(0,*)	= 0.0 		& 	Av_bsup1(Num_u1-1,*)	= 0.0

    Av_bsup1(*,0)	= (Shift(bsup1(*,1),1) + Shift(bsup1(*,1),-1))/2									; for b3 field aligned at North Ionosphere
    Av_bsup1(0,0)	= 0.0		& 	    Av_bsup1(Num_u1-1,0)	= 0.0

    Av_bsup1(*,Num_u3-1)	= (Shift(bsup1(*,Num_u3-2),1) + Shift(bsup1(*,Num_u3-2),-1))/2 				; for b3 field aligned at Southern Ionosphere
    Av_bsup1(0,Num_u3-1)	= 0.0		& 	    Av_bsup1(Num_u1-1,Num_u3-1)	= 0.0


    Av_bsup3 		= (Shift(bsup3,1,1) + Shift(bsup3,1,-1) + Shift(bsup3,-1,1) + Shift(bsup3,-1,-1))/ 4.d0
    Av_bsup3(*,0)	= 0.0 		& 	Av_bsup3(*,Num_u3-1)	= 0.0
    Av_bsup3(0,*)	= 0.0 		& 	Av_bsup3(Num_u1-1,*)	= 0.0

; rotate to bi
    b1 = g11*bsup1 + g13*Av_bsup3
    b2 = g22*bsup2
    b3 = g13*Av_bsup1 + g33*bsup3


;    b2[0,*] 		= (4.0*b2(2,*) - b2(4,*))/3.0					; Bob's fortran
;	b2[num_u1-1,*]	= (4.0*b2(num_u1-3,*) - b2(num_u1-5,*))/3.0    	; Bob's fortran
;    b3[0,*] 		= (4.0*b3(2,*) - b3(4,*))/3.0					; Bob's fortran
    b3[num_u1-1,*] 	= Driver*evens_u3     		 					; Outer L Shell Boundary add driver boundary here


; shift bi's
    b1p3 	= shift(b1,0,-1)
    b1m3 	= shift(b1,0,1)
    b1p3[*,num_u3-1]= 0.0
    b1m3[*,0] 		= 0.0

    b2p1 	= shift(b2,-1,0)
    b2m1 	= shift(b2,1,0)
    b2p1[num_u1-1,*]= 0.0
    b2m1[0,*]		= 0.0

    b2p3 	= shift(b2,0,-1)
    b2m3 	= shift(b2,0,1)
    b2p3[*,num_u3-1]= 0.0
    b2m3[*,0] 		= 0.0

    b3p1 	= shift(b3,-1,0)
    b3m1 	= shift(b3,1,0)
    b3p1[num_u1-1,*]= 0.0
    b3m1[0,*]		= 0.0

    f1 		= f1b2*(b2p3-b2m3) + f1b3*b3            						; 	at e1 location used in e1 equation
    f1(*,0)	= 0.0															; clean up at ends of field lines
    f1(*,Num_u3-1)	= 0.0

    Av_f1	=(shift(f1,1,0) + shift(f1,-1,0))/2								;  	at e2 location used in e2 equation
    Av_f1(0,*)			= 0.0
    Av_f1(Num_u1-1,*)	= 0.0

    f2 		= f2b1*(b1p3-b1m3) + f2b3*(b3p1-b3m1)							; 	at e2 location used in e2 equation
    f2(*,0)	= 0.0
    f2(*,Num_u3-1)	= 0.0

    Av_f2	=(shift(f2,1,0) + shift(f2,-1,0))/2
    Av_f2(0,*)			= 0.0
    Av_f2(Num_u1-1,*)	= 0.0

;	Need to shift e's to all grid lvls with e1 and e2
	Av_esup1 			= (shift(esup1,1,0) + shift(esup1,-1,0))/2 			; 	need e1 at e2 locations (e2 is not on boundaries)
    Av_esup1(0,*) 		= 0.0												;	Inner L shells is not required here (cleaning up from shifts)
    Av_esup1(Num_u1-1,*)= 0.0												;	Outer L shells is not required here (cleaning up from shifts)

	Av_esup2 			= (shift(esup2,1,0) + shift(esup2,-1,0))/2			; 	need e2 at e1 locations
    Av_esup2(0,*)		= 0.0												;	Inner L shells PEC BC so Esup2 is zero
    Av_esup2(Num_u1-1,*)= esup2(Num_u1-2,*)									;	Outer L shells (equivalent to De2du1 = 0)

; 	advance eperp
; -----------------------------------------------------------
    esup1 = (e1e1*   esup1 + e1e2*Av_esup2 + e1f1*   f1 + e1f2*Av_f2)
    esup2 = (e2e1*Av_esup1 + e2e2*   esup2 + e2f1*Av_f1 + e2f2*   f2)


;	esup1[0,*] 			= (4.0*esup1(2,*) - esup1(4,*))/3.0					; Bob's fortran
;    esup1[num_u1-1,*] 	= (4.0*esup1(num_u1-3,*) - esup1(num_u1-5,*))/3.0	; Bob's fortran


    Av_b3 		= (Shift(b3,1,1) + Shift(b3,1,-1) + Shift(b3,-1,1) + Shift(b3,-1,-1))/ 4.d0		; 	need b3 at e3 (and b1) locations (not on any boundaries)
    Av_b3(*,0)	= 0.0 		& 	Av_b3(*,Num_u3-1)	= 0.0
    Av_b3(0,*)	= 0.0 		& 	Av_b3(Num_u1-1,*)	= 0.0

;	Need to put the db2du3 derivatives on the right grid point for a yee grid and the fact we have a non orthognal system
	db2_d3 		= (b2p3-b2m3)
    db2_d3(*,0)	= 0.0 		& 	db2_d3(*,Num_u3-1)	= 0.0			; 	Clean up after shifts (not on any boundaries)
;    db2_d3(0,*)	= 0.0 		& 	db2_d3(Num_u1-1,*)	= 0.0

    Av_db2_d3	= (Shift(db2_d3,1,1) + Shift(db2_d3,1,-1) + Shift(db2_d3,-1,1) + Shift(db2_d3,-1,-1))/ 4.d0

	Av_db2_d3(*,1)			= (Shift(db2_d3(*,2),-1) + Shift(db2_d3(*,2),1))/2.0
	Av_db2_d3(*,Num_u3-2)	= (Shift(db2_d3(*,Num_u3-3),-1) + Shift(db2_d3(*,Num_u3-3),1))/2.0

	Av_db2_d3(*,0)	= 0.0 		& 	Av_db2_d3(*,Num_u3-1)	= 0.0 	; 	Clean up after shifts (not on any boundaries)
    Av_db2_d3(0,*)	= 0.0 		& 	Av_db2_d3(Num_u1-1,*)	= 0.0

; 	advance eparallel
; 	-----------------------------------------------------------
    e3 		= e3b21*(b2p1-b2m1) + e3b23*(Av_db2_d3) + e3b1*b1  + e3b3*Av_b3

    e3[0,*] = 0.0	&     e3[num_u1-1,*] = 0.0 						; 	Clean up after shifts (not on any boundaries)
 	e3[*,0] = 0.0	&     e3[*,Num_u3-1] = 0.0

;;;	killing eparallel
;    e3[*,*] =0.0d0         ; for testing

;; 	average ei's
;    Av_esup1 	= (Shift(esup1,1,1) + Shift(esup1,1,-1) + Shift(esup1,-1,1) + Shift(esup1,-1,-1))/ 4.d0
;    Av_esup1(*,0)	= 0.0 		& 	Av_esup1(*,Num_u3-1)	= 0.0			; Clean up after shifts due to odd number of points
;    Av_esup1(0,*)	= 0.0 		& 	Av_esup1(Num_u1-1,*)	= 0.0

    Av_e3 		= (Shift(e3,1,1) + Shift(e3,1,-1) + Shift(e3,-1,1) + Shift(e3,-1,-1))/ 4.d0						; ################# may need to do ionospheric sheet!
    Av_e3(*,0)	= 0.0 		& 	Av_e3(*,Num_u3-1)	= 0.0
    Av_e3(0,*)	= 0.0 		& 	Av_e3(Num_u1-1,*)	= 0.0

; 	rotate to ei
    e1 = e1esup1*esup1 + e1e3*Av_e3
    e2 = g22*esup2

    e1(0,0)=(4*e1(2,0)-e1(4,0))/3.0											; Bob's fortran
    e1(num_u1-1,0)=(4*e1(num_u1-3,0)-e1(num_u1-5,0))/3.0

    esup1(0,0)			=(4*esup1(2,0)-esup1(4,0))/3.0						; Bob's fortran
    esup1(num_u1-1,0)	=(4*esup1(num_u1-3,0)-esup1(num_u1-5,0))/3.0




;; 	Thin Sheet Ionospheric Boundaries
;;	Northern Ionospheric Sheet
;	bsup3_n 	= (bsup3[*,0] + (shift(bsup3[*,0],1)+shift(bsup3[*,0],-1))/2.0) 		; Average b3 at thin sheet
;;	Extrap_1D,bsup3_n,u1,0
;	bsup3_n(0)	= bsup3[0,0]
;	bsup3_n(Num_u1-1)	= bsup3[Num_u1-1,0]

	bsup3_n 	= bsup3[Eidx,0]					; Only using B3 points to do fit
    psiatm_N 	= psiib3_N # bsup3_n  			; Nth hemis

    psiatmp 	= shift(psiatm_N,-1)
    psiatmm 	= shift(psiatm_N,1)
    b1atm_N 			= (psiatmp-psiatmm)/(dx2)
	Extrap_1D,b1atm_N,u1,0
	Extrap_1D,b1atm_N,u1,Num_u1-1
;    b1atm_N(0) 		  	= (-3.0*psiatm_N(0)        + 4.0*psiatm_N(1)        - psiatm_N(2))/(dx2)
;    b1atm_N(Num_u1-1) 	= ( 3.0*psiatm_N(Num_u1-1) - 4.0*psiatm_N(Num_u1-2) + psiatm_N(Num_u1-3))/(dx2)
    b2atm_N 			= im*psiatm_N

	b1_N    	= b1[*,1] +(shift(b1[*,1],1)+shift(b1[*,1],-1))/2
  	Extrap_1D,b1_N,u1,0
	Extrap_1D,b1_N,u1,Num_u1-1
;  	;	Bob's extrapolation
;	b1_N(0)			= 3.0*b1(1,1) - 3.5*b1(3,1) + 1.5*b1(5,1)
;	b1_N(Num_u1-1)	= 3.0*b1(Num_u1-2,1) - 3.5*b1(Num_u1-4,1) + 1.5*b1(Num_u1-6,1)

	b2_N    	= b2[*,1] +(shift(b2[*,1],1)+shift(b2[*,1],-1))/2
;	Extrap_1D,b2_N,u1,0
	b2_N(0)			=b2(0,1)
	b2_N(Num_u1-1)	=b2(Num_u1-1,1)

    b1dif_N		= b1_N-b1atm_N
    b2dif_N		= b2_N-b2atm_N

    e1[*,0] 	= (e1b1atm_N*b1dif_N + e1b2atm_N*b2dif_N)*evens_u1
    e2[*,0] 	= (e2b1atm_N*b1dif_N + e2b2atm_N*b2dif_N)*odds_u1

;;	Southern Ionospheric Sheet
;	bsup3_s 	= (bsup3[*,Num_u3-1] + (shift(bsup3[*,Num_u3-1],1)+shift(bsup3[*,Num_u3-1],-1))/2.0)
;;	Extrap_1D,bsup3_s,u1,0
;	bsup3_s(0)	= bsup3[0,Num_u3-1]
;	bsup3_s(Num_u1-1)	= bsup3[Num_u1-1,Num_u3-1]

	bsup3_s 	= bsup3[Eidx,Num_u3-1]				; Only using B3 points to do fit
    psiatm_S 	= psiib3_s # bsup3_s

    psiatmp 	= shift(psiatm_S,-1)
    psiatmm 	= shift(psiatm_S,1)

    b1atm_S 			= (psiatmp-psiatmm)/(dx2)
	Extrap_1D,b1atm_S,u1,0
	Extrap_1D,b1atm_S,u1,Num_u1-1
;    b1atm_S(0) 		  	= (-3.0*psiatm_S(0)        + 4.0*psiatm_S(1)        - psiatm_S(2))/(dx2)
;    b1atm_S(Num_u1-1) 	= ( 3.0*psiatm_S(Num_u1-1) - 4.0*psiatm_S(Num_u1-2) + psiatm_S(Num_u1-3))/(dx2)
    b2atm_S 			= im*psiatm_S


	b1_S    	= b1[*,Num_u3-2] +(shift(b1[*,Num_u3-2],1)+shift(b1[*,Num_u3-2],-1))/2
    Extrap_1D,b1_S,u1,0
    Extrap_1D,b1_S,u1,Num_u1-1
;;	Bob's extrapolation
;	b1_S(0)			= 3.0*b1(1,Num_u3-2) - 3.5*b1(3,Num_u3-2) + 1.5*b1(5,Num_u3-2)
;	b1_S(Num_u1-1)	= 3.0*b1(Num_u1-2,Num_u3-2) - 3.5*b1(Num_u1-4,Num_u3-2) + 1.5*b1(Num_u1-6,Num_u3-2)


	b2_S    	= b2[*,Num_u3-2] +(shift(b2[*,Num_u3-2],1)+shift(b2[*,Num_u3-2],-1))/2
;    Extrap_1D,b2_S,u1,0
	b2_S(0)			= b2(0,Num_u3-2)
	b2_S(Num_u1-1)	= b2(Num_u1-1,Num_u3-2)

    b1dif_S		= b1_S - b1atm_S
    b2dif_S		= b2_S - b2atm_S

    e1[*,num_u3-1] 	= (e1b1atm_S*b1dif_S + e1b2atm_S*b2dif_S)*evens_u1
    e2[*,num_u3-1] 	= (e2b1atm_S*b1dif_S + e2b2atm_S*b2dif_S)*odds_u1


;; temporary ionospheric bc - for testing perfectly reflecting Ionosphere (thin sheet)
;  e1[*,0]=0.0d0
;  e2[*,0]=0.0d0
;
;  e1[*,Num_u3-1]=0.0d0
;  e2[*,Num_u3-1]=0.0d0
;
;; -----------------------------------------

; This rotation folds the field align E into the Ionospheric boundary (thin sheet) - MDS
;#######


;	Av_e3_N = (Shift(e3[*,1],-1) + Shift(e3[*,1],1))/2.0
;;    esup1[*,0] = gsup11[*,0]*e1[*,0]+gsup13[*,0]*e3[*,1]
;
;    esup1[*,0] = gsup11[*,0]*e1[*,0]+gsup13[*,0]*Av_e3_N
;    esup2[*,0] = gsup22[*,0]*e2[*,0]
;
;	Av_e3_S = (Shift(e3[*,num_u3-2],-1) + Shift(e3[*,num_u3-2],1))/2.0
;;    esup1[*,num_u3-1] = gsup11[*,num_u3-1]*e1[*,num_u3-1]+gsup13[*,num_u3-1]*e3[*,num_u3-2]
;
;    esup1[*,num_u3-1] = gsup11[*,num_u3-1]*e1[*,num_u3-1]+gsup13[*,num_u3-1]*Av_e3_S
;    esup2[*,num_u3-1] = gsup22[*,num_u3-1]*e2[*,num_u3-1]
;

;;	BC on e1's (Bob's Email) st de1du1 = 0
;;	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ???
;e1(0,0)			=(4.0*e1(2,0)-e1(4,0))/3.0
;e1(num_u1-1,0)	=(4.0*e1(num_u1-3,0)-e1(num_u1-5,0))/3.0
;
;esup1(0,0)			=(4.0*esup1(2,0)-esup1(4,0))/3.0
;esup1(num_u1-1,0)	=(4.0*esup1(num_u1-3,0)-esup1(num_u1-5,0))/3.0
;
;e1(0,Num_u3-1)			=(4.0*e1(2,Num_u3-1)-e1(4,Num_u3-1))/3.0
;e1(num_u1-1,Num_u3-1)	=(4.0*e1(num_u1-3,Num_u3-1)-e1(num_u1-5,Num_u3-1))/3.0
;
;esup1(0,Num_u3-1)			=(4.0*esup1(2,Num_u3-1)-esup1(4,Num_u3-1))/3.0
;esup1(num_u1-1,Num_u3-1)	=(4.0*esup1(num_u1-3,Num_u3-1)-esup1(num_u1-5,Num_u3-1))/3.0


; write output and plot
bnu_M 	= 0.0 & bph_M 	= 0.0 & bmu_M 	= 0.0
enu_M 	= 0.0 & eph_M 	= 0.0 & emu_M 	= 0.0
if (t ge tpnext) then begin
iplot=iplot+1
bnu_M 	= 1.0e9*bsup1*h1/Re        ; B in nT:  mu, nu and phi components
bph_M 	= 1.0e9*bsup2*h2/Re
bmu_M 	= 1.0e9*b3/h3/Re
enu_M 	= 1.0e3*esup1*h1           ; E in mV/m
eph_M 	= 1.0e3*esup2*h2
emu_M 	= 1.0e3*e3/h3

;b1_Nsc 		= 1.0e9*b1_N/hthatm_N/Re
;b2_Nsc 		= 1.0e9*b2_N/hphiatm_N/Re
;b1atm_Nsc 	= 1.0e9*b1atm_N/hthatm_N/Re
;b2atm_Nsc 	= 1.0e9*b2atm_N/hphiatm_N/Re
;
;b1_Ssc 		= 1.0e9*b1_S/hthatm_S/Re
;b2_Ssc 		= 1.0e9*b2_S/hphiatm_S/Re
;b1atm_Ssc 	= 1.0e9*b1atm_S/hthatm_S/Re
;b2atm_Ssc 	= 1.0e9*b2atm_S/hphiatm_S/Re


Iter_output,   	Enu_M, Eph_M, Emu_M, Bnu_M, Bph_M, Bmu_M, Va_arr, Plot_fields, t, iplot, tpnext, x_arr, y_arr, out_pth,  im, Phi,$
				b1_N, b1atm_N, b2_N, b2atm_N,  b1_S, b1atm_S, b2_S, b2atm_S, psiatm_N, psiatm_S, $
				Num_u1, Num_u3,  bsup3_n,  bsup3_s, hratm_N, hratm_S, psigb3_N, psigb3_S, $
				hthg_N, hthg_S, hphig_N, hphig_S, dx2, Eidx

tpnext=tpnext+dtplot                ; dtplot is interval for plots


;test for crash
if(max(abs(emu_M)) gt 1.0e12 or max(abs(eph_M)) gt 1.0e12) then stop
nanx=where(finite(emu_M) eq 0,cx)
nany=where(finite(eph_M) eq 0,cy)
if(cx gt 0 or cy gt 0) then stop



end


endwhile			; End of time loop

Finishtime = SYSTIME(/SECONDS)
;Print,'Average Time per time iteration = ',(Finishtime - Intstarttime)/Nt
Print,'Total time taken = ',Finishtime - Starttime,' seconds'

close,/all
!p.multi=0

stop
End			; End of Iterate
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


Pro Iter_output,Enu, Eph, Emu, Bnu, Bph, Bmu, Va_arr, Plot_fields, t, iplot, tpnext, x_arr, y_arr, out_pth,  im, Phi,$
				b1_N, b1atm_N, b2_N, b2atm_N,  b1_S, b1atm_S, b2_S, b2atm_S, psiatm_N, psiatm_S, $
				Num_u1, Num_u3,  bsup3_n,  bsup3_s, hratm_N, hratm_S, psigb3_N, psigb3_S, $
				hthg_N, hthg_S, hphig_N, hphig_S, dx2, Eidx


;	Plot Results
if(Plot_fields Eq 1) then begin

;	Calculate ground signatures
psig = psigb3_N # bsup3_N		           ; Nth hemis
psigp = shift(psig,-1)
psigm = shift(psig,1)
psigp[num_u1-1] = psig[num_u1-3] + 3.d0*(psig[num_u1-1]-psig[num_u1-2])
psigm[0] = psig[2] + 3.d0*(psig[0]-psig[1])
bxg_N = -(psigp-psigm)/dx2/hthg_N
byg_N = im*psig/hphig_N

psig = psigb3_S # bsup3_S		           ; Sth hemis
psigp = shift(psig,-1)
psigm = shift(psig,1)
psigp[num_u1-1] = psig[num_u1-3] + 3.d0*(psig[num_u1-1]-psig[num_u1-2])
psigm[0] = psig[2] + 3.d0*(psig[0]-psig[1])
bxg_S = (psigp-psigm)/dx2/hthg_S
byg_S = im*psig/hphig_S


Label = StrTrim(string(format='(I4)',iplot),1)
If iplot LT 1 Then Label=Label
If iplot LT 10 Then Label='0'+Label
If iplot LT 100 Then Label='0'+Label
If iplot LT 1000 Then Label='0'+Label

filenme = 'Faligned2D_'+Label+'.sav'
save,filename=out_pth+filenme,	t, Num_u1, Num_u3, Enu, Eph, Emu, Bnu, Bph, Bmu, b1_N, b1atm_N, b2_N, b2atm_N, $
						b1_S, b1atm_S, b2_S, b2atm_S, psiatm_N, psiatm_S, im, bsup3_n,  bsup3_s, hratm_N, hratm_S, $
						bxg_N, byg_N, bxg_S, byg_S

bnu 	= bnu + (Shift(bnu,1,0)+Shift(bnu,-1,0))/2
bnu		= bnu + (Shift(bnu,0,1)+Shift(bnu,0,-1))/2
bnu(0,*)= 0.0 & bnu(Num_u1-1,*) = 0.0 & bnu(*,0) = 0.0 & bnu(*,Num_u3-1) = 0.0

bph 	= bph + (Shift(bph,1,0)+Shift(bph,-1,0))/2
bph	= bph + (Shift(bph,0,1)+Shift(bph,0,-1))/2
bph(0,*) = 0.0 & bph(Num_u1-1,*) = 0.0 & bph(*,0) = 0.0 & bph(*,Num_u3-1) = 0.0

bmu 	= bmu + (Shift(bmu,1,0)+Shift(bmu,-1,0))/2
bmu	= bmu + (Shift(bmu,0,1)+Shift(bmu,0,-1))/2
bmu(0,*) = 0.0 & bmu(Num_u1-1,*) = 0.0
bmu(*,0) = 0.0 & bmu(*,Num_u3-1) = 0.0

enu 	= enu + (Shift(enu,1,0)+Shift(enu,-1,0))/2
enu	= enu + (Shift(enu,0,1)+Shift(enu,0,-1))/2
enu(0,*) = 0.0 & enu(Num_u1-1,*) = 0.0 & enu(*,0) = 0.0 & enu(*,Num_u3-1) = 0.0

eph 	= eph + (Shift(eph,1,0)+Shift(eph,-1,0))/2
eph	= eph + (Shift(eph,0,1)+Shift(eph,0,-1))/2
eph(0,*) = 0.0 & eph(Num_u1-1,*) = 0.0 & eph(*,0) = 0.0 & eph(*,Num_u3-1) = 0.0

emu 	= emu + (Shift(emu,1,0)+Shift(emu,-1,0))/2
emu	= emu + (Shift(emu,0,1)+Shift(emu,0,-1))/2
emu(0,*) = 0.0 & emu(Num_u1-1,*) = 0.0 & emu(*,0) = 0.0 & emu(*,Num_u3-1) = 0.0






device,decomposed=0
;loadct,4,/silent
window,6,xsize=1200,ysize=800,title='Fields time = '+StrTrim(t,2)+' '+StrTrim(iplot)
!P.charsize=1.6
!P.Charthick=1.2
!p.symsize=0.2

!p.multi=[0,3,2]



;contour,ABS(Enu),/fill,nlevels=255,xstyle=1,ystyle=1,title='Enu'
;contour,ABS(Eph),/fill,nlevels=255,xstyle=1,ystyle=1,title='Eph'
;contour,ABS(Emu),/fill,nlevels=255,xstyle=1,ystyle=1,title='Emu'
;
;contour,ABS(Bnu),/fill,nlevels=255,xstyle=1,ystyle=1,title='Bnu'
;contour,ABS(Bph),/fill,nlevels=255,xstyle=1,ystyle=1,title='Bph'
;contour,ABS(Bmu),/fill,nlevels=255,xstyle=1,ystyle=1,title='Bmu'



Np=50.5			; Number of Contour levels (/2)


vamax=max(abs(enu))
tlab=string(format="('Enu, t=',f7.1)",t)
print,'t=',t
print,"max Enu=",vamax
if(vamax lt 1.0e-15) then vamax=1
dcont=(vamax)/Np
clev=float( (findgen(2*Np+1)-Np)*dcont)
contour,enu*exp(im*Phi*!pi/180),x_arr,y_arr,levels=clev,/fill,xtitle='x [Re]',ytitle='z [Re]',title='Enu    max='+string(max(abs(enu))),/xstyle,/ystyle,/ISOTROPIC

vamax=max(abs(eph))
print,"max Eph=",vamax
if(vamax lt 1.0e-15) then vamax=1
dcont=(vamax)/Np
clev=float( (findgen(2*Np+1)-Np)*dcont)
contour,eph*exp(im*Phi*!pi/180),x_arr,y_arr,levels=clev,/fill,xtitle='x [Re]',ytitle='z [Re]',title='Eph    max='+string(max(abs(eph))),/xstyle,/ystyle,/ISOTROPIC


tlab=string(format="('Emu, t=',f7.1)",t)
vamax=max(abs(emu))
print,"max Emu=",vamax
if(vamax lt 1.0e-15) then vamax=1
dcont=(vamax)/Np
clev=float( (findgen(2*Np+1)-Np)*dcont)
contour,emu*exp(im*Phi*!pi/180),x_arr,y_arr,levels=clev,/fill,xtitle='x [Re]',ytitle='z [Re]',title=tlab+'  max='+string(max(abs(emu))),/xstyle,/ystyle,/ISOTROPIC


 vamax=max(abs(bnu))
print,"max Bnu=",vamax
if(vamax lt 1.0e-15) then vamax=1
dcont=(vamax)/Np
clev=float( (findgen(2*Np+1)-Np)*dcont)
contour,bnu*exp(im*Phi*!pi/180),x_arr,y_arr,levels=clev,/fill,xtitle='x [Re]',ytitle='z [Re]',title='Bnu    max='+string(max(abs(bnu))),/xstyle,/ystyle,/ISOTROPIC

vamax=max(abs(bph))
print,"max Bph=",vamax
if(vamax lt 1.0e-15) then vamax=1
dcont=(vamax)/Np
clev=float( (findgen(2*Np+1)-Np)*dcont)
contour,bph*exp(im*Phi*!pi/180),x_arr,y_arr,levels=clev,/fill,xtitle='x [Re]',ytitle='z [Re]',title='Bph    max='+string(max(abs(bph))),/xstyle,/ystyle,/ISOTROPIC

vamax=max(abs(bmu))
print,"max Bmu=",vamax
if(vamax lt 1.0e-15) then vamax=1
dcont=(vamax)/Np
clev=float( (findgen(2*Np+1)-Np)*dcont)
contour,bmu*exp(im*Phi*!pi/180),x_arr,y_arr,levels=clev,/fill,xtitle='x [Re]',ytitle='z [Re]',title='Bmu    max='+string(max(abs(bmu))),/xstyle,/ystyle,/ISOTROPIC


Filename 	= out_pth+'Fields'+Label+'.png'
image 		= TVRD(0,0,!d.X_size,!d.Y_size,true=1)
Write_PNG,filename,image,r,g,b



window,7,xsize=1200,ysize=800,title='Ionopsheric Fields time = '+StrTrim(t,2)+' '+StrTrim(iplot)

!P.charsize=2.6
!P.Charthick=1.2
!p.symsize=0.2

!p.multi=[0,4,2]

Label = StrTrim(string(format='(I4)',iplot),1)
If iplot LT 1 Then Label=Label
If iplot LT 10 Then Label='0'+Label
If iplot LT 100 Then Label='0'+Label
If iplot LT 1000 Then Label='0'+Label


Bmax = max( [max(ABS(b1_N)), max(ABS(b2_N))] )
 Plot, Double(b1_N),Title='Magnetospheric Bs (NH)'+string(format="('t=',f7.1)",t),Yrange=[-bmax,bmax]
oPlot, Imaginary(b1_N),linestyle=1
oPlot, Double(b2_N),Color=150
oPlot, Imaginary(b2_N),Color=150,linestyle=1
oPlot, Eidx,Double(bsup3_N*hratm_N),Color=50
oPlot, Eidx,Imaginary(bsup3_n*hratm_N),Color=50,linestyle=1


Bmax = max( [max(ABS(b1atm_N)), max(ABS(b2atm_N))] )
 Plot,Double(b1atm_N),Title='Atmospheric Bs (NH)',Yrange=[-bmax,bmax]
oPlot,Imaginary(b1atm_N),linestyle=1
oPlot,Double(b2atm_N),Color=150
oPlot,Imaginary(b2atm_N),Color=150,linestyle=1

Bmax = max( [max(ABS(bxg_N)), max(ABS(byg_N))] )
 Plot,Double(bxg_N),Title='Ground Bs (NH)',Yrange=[-bmax,bmax]
oPlot,Imaginary(bxg_N),linestyle=1
oPlot,Double(byg_N),Color=150
oPlot,Imaginary(byg_N),Color=150,linestyle=1

Bmax = max(max(ABS(Psiatm_N)))
 Plot,Double(psiatm_N),Title='Psi (NH)',Yrange=[-bmax,bmax]
oPlot,Imaginary(psiatm_N),linestyle=1

Bmax = max( [max(ABS(b1_S)), max(ABS(b2_S))] )
 Plot, Double(b1_S),Title='Magnetospheric Bs (SH)',Yrange=[-bmax,bmax]
oPlot, Imaginary(b1_S),linestyle=1
oPlot, b2_S,Color=150
oPlot, Imaginary(b2_S),Color=150,linestyle=1
oPlot, Eidx,Double(bsup3_s*hratm_S),Color=50
oPlot, Eidx,Imaginary(bsup3_s*hratm_S),Color=50,linestyle=1

Bmax = max( [max(ABS(b1atm_S)), max(ABS(b2atm_S))] )
 Plot,Double(b1atm_S),Title='Atmospheric Bs (SH)',Yrange=[-bmax,bmax]
oPlot,Imaginary(b1atm_S),linestyle=1
oPlot,Double(b2atm_S),Color=150
oPlot,Imaginary(b2atm_S),Color=150,linestyle=1

Bmax = max( [max(ABS(bxg_S)), max(ABS(byg_S))] )
 Plot,Double(bxg_S),Title='Ground Bs (SH)',Yrange=[-bmax,bmax]
oPlot,Imaginary(bxg_S),linestyle=1
oPlot,Double(byg_S),Color=150
oPlot,Imaginary(byg_S),Color=150,linestyle=1

Bmax = max(max(ABS(Psiatm_S)))
 Plot,Double(psiatm_S),Title='Psi (SH)',Yrange=[-bmax,bmax]
oPlot,Imaginary(psiatm_S),linestyle=1

Filename 	= out_pth+'IonosGnd'+Label+'.png'
image 		= TVRD(0,0,!d.X_size,!d.Y_size,true=1)
Write_PNG,filename,image,r,g,b

endif



;;test for crash
;if(max(abs(emu)) gt 1.0e12 or max(abs(eph)) gt 1.0e12) then stop
;nanx=where(finite(emu) eq 0,cx)
;nany=where(finite(eph) eq 0,cy)
;if(cx gt 0 or cy gt 0) then stop


end 			; End of Iterate Ouput

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Pro Find_flr,LVal,sizen,DEN,har
; sizen=200
 harm = 8                     ;	no. of harmonics to calc.
 D = FLTARR(23)               ;	constants and reference parameters to pass to functions
; Lmbd = DBLARR(sizen)         ; Density values - or lamda to load DENSITY
 HAR = FLTARR(harm)           ; harmonics
 Ref = 6378.e3                ; Earth radius
 Anchor = 130.e3              ;	Altitude of Ionosphere anchor point
 Lat = acos(sqrt((Ref + anchor)/(Ref*LVal)))
 h = 2*Lat/(sizen + 1)        ; step size of lambda
 mu = 1.256637061E-6          ; Permeability of free space
 D = [LVal*Ref , 7378E3 , 9.81 , 7.2921E-5 , 1.38E-23 , 8.02E15 , 2903 , $
	  12.0407E10 , 6.02035E10 , 5.418315E10 , 6.02035E9 , 0 , 0 , 0 , $
	  2.676E-26 , 1.6725E-27 , 6.69E-27 , 0 , 0 , 0 , LVal, Ref, mu]
;
; Details of array D are:
; D = [L*Re , Refhght , g at refhght , earth rot , kb , ko , ref iontemp ,
;     elec den ref, ox den ref , hyd den ref, hel den ref , 0,0,0,
;     mass ox , mass hyd , mass hel, 0,0,0, L, Re, mu ]

; FOR i = 1,sizen DO Lmbd(i-1) = -Lat+h*i     ; load lambda
; Load r^-m model
; rho_0=DEN(sizen/2)                          ; Equatorial density value
; f = cos(Lmbd)*cos(Lmbd)
;stop
; DEN(*) =  rho_0*(1.d0/f)^3

;  Solve for FLRs
 A = DBLARR(sizen,sizen)          ; Finite difference matrix
 B = DBLARR(sizen+2,sizen)
 C = DBLARR(3)
 EVAL = COMPLEXARR(sizen)      ; eigenvalues
 ORDER = INTARR(sizen)         ; harmonic order
 ltd = fltarr(sizen)
 con = double(((D(5)/(D(0)^4))^2)/D(22))
 For i = 1, sizen  do $
  Begin
  lda = -Lat + h*float(i)
  ltd(i-1) = lda
  C = [1.d0/h - tan(lda)/2,-2/h,1.d0/h + tan(lda)/2.d0]
  B(i-1,i-1) =  C*con/(h*DEN(i-1)*(cos(lda)^14.))
 ENDFOR
 A = B(1:sizen,0:sizen-1)             ;Truncation to square matrix
 EVAL = HQR(ELMHES(A),/double)  ;Calculates eigenvalues
;EVEC = EIGENVEC(A, EVAL, /DOUBLE);          Calculates eigenvectors
 ORDER = SORT(EVAL)             ;Sort the eigen values small to large
 For i = 0,harm-1 do HAR(i) = sqrt(-1.d0*FLOAT(EVAL(ORDER(i))))*500.d0/!dpi  ; freqs in mHz
 ;Print,'Freqs:',LVal,HAR

END
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




