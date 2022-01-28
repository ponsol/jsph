
module sphcont

  contains

   function checkstop() result(oval)
    use partmod
    integer  :: oval
      oval = 1

       if (spvars%giter >= spvars%gitermax) oval = 0

       if ( abs(spvars%gtime - spvars%gtmax)/spvars%gtmax <= 1.e-7 ) oval = 0

       if ( spvars%gtime >= spvars%gtmax ) oval = 0

   end function
end module 


program sph

 use partmod
 use geom
 use sphsub
 use problem
 use sphcont
 implicit none

 integer :: iter
 integer :: stopiter

       call initgeom()
       call initprob()
       call initfinal()
       call getnblist()
       call fillvars()


       call outallpart()
       call outonecut(1,1,1)

    stopiter = 1
    do while ( stopiter  == 1 .and. spvars%gitermax > 0 )

       spvars%giter = spvars%giter + 1
       call get_dden()    
       call get_dvel()
       call get_de()
       call getdt(spvars%gdt) 

        write(*,*) '#iter ', spvars%giter, '  dt ', spvars%gdt, ' time ', spvars%gtime
       if ( spvars%gdt + spvars%gtime >= spvars%gtmax ) then
            spvars%gdt = spvars%gtmax - spvars%gtime 
       endif

       call allupdate()

       !call boundary()
       spvars%gtime = spvars%gtime + spvars%gdt ;

       if ( mod(spvars%giter, 1) == 0 ) then
         write(*,*) '#iter ', spvars%giter, '  dt ', spvars%gdt, ' time ', spvars%gtime
       endif  
       stopiter =  checkstop()
    enddo

       call outallpart()
       call outonecut(1,1,1)

end program sph


