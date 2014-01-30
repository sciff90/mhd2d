pro idl_plots
  num_u1 = long(51)
  num_u3 = long(371)

  dx2 = 0.0d

  data_dir = './data/'

  
  va_arr = dblarr(num_u1,num_u3)
  x_arr = dblarr(num_u1,num_u3)
  y_arr = dblarr(num_u1,num_u3)
  Eidx = lonarr(num_u1/2+1)
  print,'arrays assigned!!'
  openr,unit1, data_dir +'plotting_constants.dat',/get_lun,/F77_UNFORMATTED
  readu,unit1,num_u1,num_u3,va_arr,x_arr,y_arr,Eidx,dx2
  close,unit1
  free_lun,unit1

  print,'Plotting Constants Assigned'

  enu = dcomplexarr(num_u1,num_u3)
  eph = dcomplexarr(num_u1,num_u3)
  emu = dcomplexarr(num_u1,num_u3)
  bnu = dcomplexarr(num_u1,num_u3)
  bph = dcomplexarr(num_u1,num_u3)
  bmu = dcomplexarr(num_u1,num_u3)

  b1_n = dcomplexarr(num_u1)
  b1atm_n = dcomplexarr(num_u1)
  b2_n = dcomplexarr(num_u1)
  b2atm_n = dcomplexarr(num_u1)
  b1_S = dcomplexarr(num_u1)
  b1atm_s = dcomplexarr(num_u1)
  b2_s = dcomplexarr(num_u1)
  b2atm_s = dcomplexarr(num_u1)

  psiatm_n = dcomplexarr(num_u1)
  psiatm_s = dcomplexarr(num_u1)

  hratm_n = dblarr(num_u1)
  hratm_s = dblarr(num_u1)

  bxg_n = dcomplexarr(num_u1)
  bxg_s = dcomplexarr(num_u1)
  byg_n = dcomplexarr(num_u1)
  byg_s = dcomplexarr(num_u1)

  im = dcomplex(0.0)
  bsup3_n = dcomplexarr(num_u1/2+1)
  bsup3_s = dcomplexarr(num_u1/2+1)
  tt = 0.0d

  openr,unit1, data_dir +'plot_data0010.dat',/get_lun,/F77_UNFORMATTED
  readu,unit1,tt,$
       Enu, Eph, Emu, Bnu, Bph, Bmu, $
       b1_N, b1atm_N, b2_N, b2atm_N, $
       b1_S, b1atm_S, b2_S, b2atm_S, $
       psiatm_N, psiatm_S, $
       im, $
       bsup3_n,  bsup3_s, $
       hratm_N, hratm_S, $
       bxg_N, byg_N, bxg_S, byg_S
  close,unit1
  free_lun,unit1
  print,'data0001.dat loaded!!!'
  stop

end
