program main
       use global
       implicit none 
       real etime
       real elapsed(2)
       real total

       call gen_grid()
       !call gen_va()
       !call gen_basisFn()
       !call gen_numfactors()
       !call iterate()

       write (*, *) "Program finished in",total       
end program main

