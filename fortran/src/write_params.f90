subroutine write_params
  use global
  implicit none
  !Output Run Parameters

  character(len=100) :: filename
  integer :: reclength

  filename='./data/run_parameters.dat'
  inquire(iolength=reclength) &
    num_u1,num_u3,&
    x_arr,y_arr,z_arr,r_arr,u1,u3,&
    g11,g22,g33,g13,gsup13,gsup11,gsup22,gsup33,costh,dmudx,&
    h1,h2,h3,h30,Jac,dt,nt,&
    va_arr,rho_a,&
    m,r_iono,cour,dtplot,drive0,phi,Freq,z0,ds0,dsn,modefrac,&
    eta_arr,sigp_arr,sigh_arr,&
    sigpatm_N,sighatm_N,sig0atm_N,&
    sigpatm_S,sighatm_S,sig0atm_S
    


  open(unit=10, file=filename,&
                form='unformatted',&
                access='direct',&
                recl= reclength)

  write(10,rec=1)&
    num_u1,num_u3,&
    x_arr,y_arr,z_arr,r_arr,u1,u3,&
    g11,g22,g33,g13,gsup13,gsup11,gsup22,gsup33,costh,dmudx,&
    h1,h2,h3,h30,Jac,dt,nt,&
    va_arr,rho_a,&
    m,r_iono,cour,dtplot,drive0,phi,Freq,z0,ds0,dsn,modefrac,&
    eta_arr,sigp_arr,sigh_arr,&
    sigpatm_N,sighatm_N,sig0atm_N,&
    sigpatm_S,sighatm_S,sig0atm_S
  close(unit=10)
  write(*,*) 'File Written to:',filename

end subroutine write_params
