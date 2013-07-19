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

eval_sc				= La_eigenproblem((a),eigenvectors=evec1_sc,/double)
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

set_plot,'win'
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
