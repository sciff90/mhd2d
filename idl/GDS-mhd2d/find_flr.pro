
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
