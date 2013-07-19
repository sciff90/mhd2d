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





set_plot,'win'
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

