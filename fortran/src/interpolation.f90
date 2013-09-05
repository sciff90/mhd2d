module interpolation
implicit none
contains
!Subroutine to do a quadratic fit

!V are the function values at known points X
!U is the array containing the x values that the interpolation values are wanted at
!Y is the output array
subroutine quad_inter(V,X,U,Y)
    implicit none
    real,dimension(:) :: U,X,V,Y   
    real:: A,B,C,a0,a1,a2    !Y = A*x^2 + B*x +C

    !Counters
    integer ii,jj,kk;
    
    do ii = 1,size(U)

            do jj = 1, Size(X)
                    if (X(jj)>=U(ii)) then
                            exit
                    end if
             end do
            
        jj = jj-1        

        a0 = V(jj)/((X(jj)-X(jj+1))*(X(jj)-X(jj+2)))
        a1 = V(jj+1) / ((X(jj+1)-X(jj))*(X(jj+1)-X(jj+2)))
        a2 = V(jj+2) / ((X(jj+2)-X(jj))*(X(jj+2)-X(jj+1)))

        A = a0+a1+a2
        B = -(a0*(X(jj+1)+X(jj+2)) + a1*(X(jj)+X(jj+2))+a2*(X(jj)+x(jj+1)))
        C = a0*X(jj+1)*X(jj+2) + a1*X(jj)*X(jj+2) + a2*X(jj)*X(jj+1)

        Y(ii) = A*U(ii)**2+ B*U(ii) + C

   end do
   return
end subroutine
    
end module
