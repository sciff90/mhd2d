program main
       use global
       implicit none
       double precision etime
       double precision elapsed(2)
       double precision total
       call system('ulimit -s unlimited')
       call gen_grid()
       call get_va()
       call get_facts()
       call get_basisfn()
       call write_params()
       call iterate()

       write (*, *) "Program finished in",total
end program main

