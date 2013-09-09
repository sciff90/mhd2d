subroutine get_facts()
    use global
    implicit none

    !arrays
    real,dimension(0:num_u1-1,0:Num_u3-1) :: sigp_eps_arr,sigh_eps_arr
    real,dimension(0:num_u1-1,0:Num_u3-1) :: esig_arr,esig2_arr,csig_arr,csig2_arr
    real,dimension(0:num_u1-1,0:Num_u3-1) :: ssig_arr,ssig2_arr
    real,dimension(0:num_u1-1,0:Num_u3-1) :: e1e1,e1e2,e1f1,e1f2,f1b2,e2e1,e2e2
    real,dimension(0:num_u1-1,0:Num_u3-1) :: e2f1,e2f2,f2b1,f2b3,e3b21,e3b23
    real,dimension(0:num_u1-1,0:Num_u3-1) :: b1e2,b2e1,b2e3,b3e2
    complex,dimension(0:num_u1-1,0:Num_u3-1) :: f1b3,e3b1,e3b3,b1e3,b3e1
    real,dimension(0:num_u1-1,0:Num_u3-1) :: gg12,gg21,e1esup1,e1e3

    !variables
    real :: dx2,dxsq
    !factors for integration
    sigp_eps_arr      = sigp_arr / eps_arr
    sigh_eps_arr    = sigh_arr / eps_arr

    esig_arr      = exp(-sigp_eps_arr*dt)
    esig2_arr = exp(-sigp_eps_arr*dt*0.5)
    csig_arr  = cos(sigh_eps_arr*dt)
    csig2_arr = cos(sigh_eps_arr*dt*0.5)
    ssig_arr      = sin(sigh_eps_arr*dt)
    ssig2_arr = sin(sigh_eps_arr*dt*0.5)
    
    !factors for generator
    gg12 =g22*g33/h3/jac
    gg21 =(g11*g33-g13*g13)/h3/jac

    e1e1 = esig_arr*csig_arr
    e1e2 = gg12*esig_arr*ssig_arr

    e1f1 = dt*esig2_arr*csig2_arr
    e1f2 = dt*gg12*esig2_arr*ssig2_arr

    f1b2 = -va2_arr/jac/d32
    f1b3 =  va2_arr*im/jac

    e2e1 = -gg21*esig_arr*ssig_arr
    e2e2 = e1e1

    e2f1 = -dt*gg21*esig2_arr*ssig2_arr
    e2f2 = dt*esig2_arr*csig2_arr

    f2b1      = va2_arr/jac/d32
    f2b3 = -va2_arr/jac/del_u1*0.5

    e3b1 = -im*eta_arr*g33/jac
    e3b21 = eta_arr/jac*g33/del_u1*0.5
    e3b23 =-eta_arr/jac*g13/d32
    e3b3 = im*eta_arr*g13/jac

    b1e2 = dt/d32/jac
    b1e3 =-dt*im/jac
    b2e1 =-dt/d32/jac
    b2e3 = dt/del_u1/jac*0.5

    b3e1 = dt*im/jac
    b3e2 =-dt/del_u1/jac*0.5
    
    !e1esup1 = g11/(1-g13*gsup13)
    !e1e3 =  g13*gsup33/(1-g13*gsup13)

    e1esup1 =1.0/gsup11
    !From Bob's fortran code (and word document)
    e1e3 =-gsup13/gsup11
    
    !From Bob's fortran code(and word document)
    dx2 = 2.0d0*del_u1
    dxsq = del_u1*del_u1

    
end subroutine get_facts
