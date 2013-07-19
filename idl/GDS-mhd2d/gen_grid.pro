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

;set_plot,'ps'
;device,decomposed=0
;;loadct,0
;!p.multi=0
;window,0,xsize=800,ysize=700,title='Grid'
;!P.charsize=1.6
;!P.Charthick=1.2
;!p.symsize=0.2

;plot,x_arr,y_arr,psym=2,title='2D Faligned Grid',xtitle='[R!dE!n]',ytitle='[R!dE!n]',$
;	 xstyle=1,ystyle=1,/isotropic

If Plot_png EQ 1 then $
begin
 Filename 	= out_pth+'2D_Faligned_Grid.png'
 image 		= TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
end

;set_plot,'win'
;device,decompose=0
;window,1,Title='Title Grid',Xsize=1200,Ysize=700
;!P.multi=[0,4,3]
;!P.charsize=2.5
;!P.Charthick=1.2
;Lidx = min(where (r_arr(*,num_u3_half) GE Lval))

;plot,u3,gsup11(Lidx,*),Title='g11_con',Xtitle='u3',xstyle=1
;plot,u3,gsup22(Lidx,*),Title='g22_con',Xtitle='u3',xstyle=1
;plot,u3,gsup33(Lidx,*)*dmudx(Lidx,*)^2,Title='g33_con',Xtitle='u3',xstyle=1
;plot,u3,gsup13(Lidx,*)*dmudx(Lidx,*),Title='g13_con',Xtitle='u3',xstyle=1

;plot,u3,g11(Lidx,*),Title='g11_cov',Xtitle='u3',/ylog,xstyle=1
;plot,u3,g22(Lidx,*),Title='g22_cov',Xtitle='u3',/ylog,xstyle=1
;plot,u3,g33(Lidx,*)/dmudx(Lidx,*)^2,Title='g33_cov',Xtitle='u3',/ylog,xstyle=1
;plot,u3,g13(Lidx,*)/dmudx(Lidx,*),Title='g13_cov',Xtitle='u3',xstyle=1


;plot,u3,Jac(Lidx,*)/dmudx(Lidx,*),Title='Jacobian',Xtitle='u3',/ylog,xstyle=1,ystyle=1
;plot,u3,r_arr(Lidx,*),Title='R array',Xtitle='u3',xstyle=1
;plot,u3,Acos(costh(Lidx,*))*180/!pi,Title='theta array',Xtitle='u3',xstyle=1

;plot,x_arr,y_arr,Title='Grid',psym=3,xstyle=1
;oplot,x_arr(Lidx,*),y_arr(Lidx,*),color=50,thick=2

;!p.multi=0

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

