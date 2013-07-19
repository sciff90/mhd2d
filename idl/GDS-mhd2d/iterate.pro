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





set_plot,'win'
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


