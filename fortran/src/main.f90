program main
       use global
       implicit none 
       real etime
       real elapsed(2)
       real total
       call system('ulimit -s unlimited')
       call gen_grid()
       call get_va()
       call get_facts()
       call get_basisfn()
       !call gen_numfactors()
       !call iterate()

       write (*, *) "Program finished in",total       
end program main

