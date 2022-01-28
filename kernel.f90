module kernel
   use partmod
   implicit none

   contains
 
     function kern(posi,posj, hval) result(oval)
      real(kind=8), intent(in) :: posi(3), posj(3), hval
      real (kind=8) :: oval

      integer :: kopt, ndim

      real (kind=8) :: r, h,norm, q


        h = hval

        ndim =spvars%ndim
        kopt =spvars%kopt

        r = sqrt ( dot_product( posi- posj , posi- posj ) )

        q= r/h

        oval = 0.0
        norm = 0.0

       if ( kopt == 1 ) then
          norm = 1.0/h/sqrt(3.14159)
          if ( q <= 3.0 ) then
           oval = exp( - q*q )
          endif
          oval = oval*norm

       elseif ( kopt == 2 ) then

          if ( ( q >= 0.0 ) .and. ( q <= 1.0) ) then
           oval = 1.0 - 3.0/2.0*q*q * (1.0 - q/2.0)
          elseif ( ( q > 1 ) .and. ( q <= 2.0) ) then
               oval = (2.0 - q)**3.0  / 4.0
          else
               oval = 0.0
          endif
          norm  = 2.0/3.0/h
          oval = oval *norm

       !Monaghan & Lattanzio 1985
       elseif ( kopt == 3 ) then
          if ( ( q >= 0.0 ) .and. ( q <= 1.0) ) then
             oval = 4.0 - 6.0*q*q + 3.0*q**3.0
          elseif ( ( q > 1 ) .and. ( q <= 2.0) ) then
             oval = (2.0 - q)**3.0  
          elseif ( q > 2 ) then
             oval = 0.0
          endif

          norm  = 1.0/4.0/3.14159
          norm  = norm / h**3.0

          oval = oval *norm


       ! ML1985
       elseif ( kopt == 4 ) then

          if ( ( q >= 0.0 ) .and. ( q < 0.5) ) then
             oval = 6.0*q*q*q - 6.0*q*q + 1.0 
          elseif ( ( q >= 0.5 ) .and. ( q <= 1.0) ) then
             oval = 2.0*(1.0 - q)**3.0  
          elseif ( q > 1 ) then
             oval = 0.0
          endif

          if ( ndim ==1 ) norm  = 4.0/3.0 / h
          if ( ndim ==2 ) norm  = 40.0/7.0/3.14159 / (h*h)
          if ( ndim ==3 ) norm  = 8.0/3.14159  / (h*h*h)
         
          oval = oval *norm

       ! MG1983  Super Gaussian
       elseif ( kopt == 5 ) then

          oval = exp( - q*q )
          oval = oval * ( 3.0 / 2.0 - q*q  )

          norm  = 1.0 / h / sqrt( 3.14159 )
         
          oval = oval *norm
       endif


     end function
 
     function gradkern(posi, posj, hval) result(oval)
      real(kind=8), intent(in) :: posi(3), posj(3), hval
      real (kind=8) :: oval(3), tval

      integer :: kopt, ndim

      real (kind=8) :: x, y, z
      real (kind=8) :: r, h,norm, q
      real (kind=8) :: dr(3)
      real (kind=8) :: f1, f2




        ndim = spvars%ndim
        kopt = spvars%kopt

        r = sqrt( dot_product( posi - posj ,  posi - posj  ) )
        dr = posi - posj
        !dr = abs (dr)
        h = hval
        q= r/h




        oval = 0.0
        norm = 0.0
        tval = 0.0

      if ( kopt == 1 ) then

          norm = 1.0/h/sqrt(3.14159)

      !spline
      elseif ( kopt == 3 ) then

         if ( ( q >= 0.0) .and. (q <= 1.0 ) ) then
             tval = 3.0*q*(4.0 - 3.0*q)
         elseif ( ( q > 1) .and. (q <= 2.0 ) ) then
             tval = 3.0*(2.0 - q)**2.0
         elseif ( (q > 2.0 ) ) then
             tval = 0.0
         endif
         norm =  -1.0/4/3.1415
         norm =  norm/h**4.0
         tval =  norm*tval  

       ! ML1985
       elseif ( kopt == 4 ) then

          if ( ( q >= 0.0 ) .and. ( q < 0.5) ) then
             tval = 3.0*q*q - 2.0*q 
          elseif ( ( q >= 0.5 ) .and. ( q <= 1.0) ) then
             tval = -(1.0 - q)**2.0  
          elseif ( q > 1 ) then
             tval = 0.0
          endif

          if ( ndim ==1 ) norm  = 4.0/3.0 / h
          if ( ndim ==2 ) norm  = 40.0/7.0/3.14159 / (h*h)
          if ( ndim ==3 ) norm  = 8.0/3.14159  / (h*h*h)

          norm = norm * 6.0 / h

          tval = norm * tval

       ! MG1983  Super Gaussian
       elseif ( kopt == 5 ) then

          f1 = exp( - q*q )
          f2 = ( 3.0 / 2.0 - q*q )

          tval = - 2*q/h * (f1*f2  + f1) 

          norm  = 1.0 / h / sqrt( 3.14159 )
          tval = tval *norm
       endif


      
       if ( r > 0 ) then
         oval = tval * (dr/r)
       else
         if ( ndim  == 1 ) dr = (/ 1.0, 0.0, 0.0 /)
         if ( ndim  == 2 ) dr = (/ 1.0, 1.0, 0.0 /)
         if ( ndim  == 3 ) dr = (/ 1.0, 1.0, 1.0 /)
         oval = tval * dr
       endif



    end function
       
end module
