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
Neutral(rr  ,0) = 0.     & Neutral(0,1) = Neutral(1  ,1)  & Neutral(0,2) = Neutral(1,2) & Neutral(0,3) = Neutral(1,3)
Neutral(rr  ,0) = 3000.  & Neutral(rr  ,1)= y2 & Neutral(rr  ,2) = y3 & Neutral(rr  ,3) = y4*exp(-(Neutral(rr  ,0)-Neutral(Pt,0))/ScleHt)
Neutral(rr+1,0) = 4000.	 & Neutral(rr+1,1)= y2 & Neutral(rr+1,2) = y3 & Neutral(rr+1,3) = y4*exp(-(Neutral(rr+1,0)-Neutral(Pt,0))/ScleHt)
Neutral(rr+2,0) = 5000.	 & Neutral(rr+2,1)= y2 & Neutral(rr+2,2) = y3 & Neutral(rr+2,3) = y4*exp(-(Neutral(rr+2,0)-Neutral(Pt,0))/ScleHt)
Neutral(rr+3,0) = 6000.	 & Neutral(rr+3,1)= y2 & Neutral(rr+3,2) = y3 & Neutral(rr+3,3) = y4*exp(-(Neutral(rr+3,0)-Neutral(Pt,0))/ScleHt)
Neutral(rr+4,0) = 100000  & Neutral(rr+4,1)= y2 & Neutral(rr+4,2) = y3 & Neutral(rr+4,3) = y4*exp(-(Neutral(rr+4,0)-Neutral(Pt,0))/ScleHt)
Free_Lun,u

;	Place Neutral Values on Model Grid Points (Assumes Radial Profile)

For ii = 0,Num_u1-1 do begin

Fline_No_N 	= INTERPOL( Neutral(*,3), Neutral(*,0)*1.0e3, z_arr(ii,*))
Fline_Mav 	= INTERPOL( Neutral(*,2), Neutral(*,0)*1.0e3, z_arr(ii,*))
Fline_Temp 	= INTERPOL( Neutral(*,1), Neutral(*,0)*1.0e3, z_arr(ii,*))

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



set_plot,'win'
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

set_plot,'win'
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


set_plot,'win'
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


