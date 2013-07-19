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
