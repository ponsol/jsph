module problem
   use partmod
   use eos
   use kernel
   implicit none

   contains

   subroutine problemgeom()

      spvars%gnptot    = 200
      spvars%gnbdymax   = 300
      spvars%ndim    = 1
      spvars%gnbdy   = spvars%gnbdymax 
   end subroutine 

   subroutine problemgetdt(dt) 
       real (kind=8 ),intent(inout) :: dt

       dt = 0.01d0

   end subroutine 



   subroutine initprob()
     integer :: i, j
     real (kind=8) :: dxl, dxr
     real (kind=8) :: xmin, xmax, ymin, ymax, xc;
     integer :: nl, nr, isen
     
     real (kind=8) :: rhol, rhor, vell, velr, pl, pr, enul, enur
     real (kind=8) :: agamma
     real (kind=8) :: dt, tmax
     real (kind=8) :: hl, hr
     real (kind=8) :: xav, asum

     integer :: nbl, nbr

     xmin = 0.0; xmax = 1.0;
     ymin = 0.0; ymax = 0.0;
     xc = 0.5;
     agamma = 1.4
     isen = 1

     spvars%kopt = 4   !kernel
     spvars%isartvis = 1
     spvars%isentropy = isen
     spvars%dtmethod = 2
     spvars%isproblemdt = 0

     spvars%alpha = 0.5
     spvars%beta  = 0.5


     if ( 1== 0 ) then
        rhol = 1.0; rhor = 1.0
        vell = 0.0; velr = 0.0
        pl = 1.0d0; pr = 1.0d0
     endif

     if ( 1== 1 ) then
       rhol = 1.0; rhor = 0.125
       vell = 0.0; velr = 0.0
       pl = 1.0d0; pr = 0.1d0
     endif

     tmax = 0.15 ; dt = 0.0005

     spvars%gtmax = tmax
     spvars%gitermax = 100
    

     nl = rhol /( rhor + rhol) *spvars%gnptot;
     nr = spvars%gnptot - nl;
     dxl = (xc-xmin)/(nl)
     dxr = (xmax -xc)/(nr)
     dxr = 1.0397*dxr

     hl = 0.15 
     hr = hl
     nbl = ceiling (3.0*hl/dxl)
     nbr = ceiling (3.0*hr/dxr)

     print*, '#nl, nr', nl, nr
     print*, '#dxl dxr ', dxl, dxr
     print*, '#hl, hr ', hl, hr
     print*, '#nbl, nbr ', nbl, nbr
     print*, '#dt, tmax, itermax', dt, tmax, spvars%gitermax

     if ( nbl + nbr > spvars%gnbdymax ) then
       write(*,*) 'Error: Insufficient boundary particles'
       write(*,*) 'Required boundary particles=', nbl+nbr
       write(*,*) 'Availabls boundary particles=', spvars%gnbdymax
       stop
     else
      spvars%gnbdy = nbl + nbr
     endif

       enul = eos_e_prho(isen, pl, rhol, agamma)
       enur = eos_e_prho(isen, pr, rhor, agamma)

     do i = 1, nl

         part(i)%pos(1) = xmin + dxl *(i-1)

         part(i)%pos(2) = 0.0
         part(i)%pos(3) = 0.0
         part(i)%vel(1) = vell
         part(i)%vel(2) = 0.0; part(i)%vel(3) = 0.0
         part(i)%en  =  enul
         part(i)%p  =  pl
         part(i)%mass =  1.0
         part(i)%den =  rhol
         part(i)%h = hl  
     enddo
 

     do j = 1, nr
         i = (spvars%gnptot - nr) + j

         part(i)%pos(1) = xc + dxr * j

         part(i)%pos(2) = 0.0
         part(i)%pos(3) = 0.0
         part(i)%vel(1) = velr
         part(i)%vel(2) = 0.0; part(i)%vel(3) = 0.0
         part(i)%en =  enur
         part(i)%p =  pr
         part(i)%den =  rhor
         part(i)%h = hr
         part(i)%mass =  1.0 
     enddo


    !right boundary
     do j = 1, nbr
        i  = j + spvars%gnptot
         part(i)%pos(1) = xc + dxr*nr + dxr*j
         part(i)%pos(2) = 0.0
         part(i)%pos(3) = 0.0
         part(i)%vel(1) = velr
         part(i)%vel(2) = 0.0; part(i)%vel(3) = 0.0
         part(i)%en =  enur
         part(i)%p =  pr
         part(i)%den =  rhor
         part(i)%h = hr
         part(i)%mass =  1.0
     enddo
 

    !left boundary
     do j = 1, nbl
        i  = j + spvars%gnptot + nbr

         part(i)%pos(1) = xmin - dxl*j

         part(i)%pos(2) = 0.0
         part(i)%pos(3) = 0.0
         part(i)%vel(1) = vell
         part(i)%vel(2) = 0.0; part(i)%vel(3) = 0.0
         part(i)%en  =  enul
         part(i)%p  =  pl
         part(i)%den =  rhol
         part(i)%h = hl  
         part(i)%mass =  1.0
     enddo

      !scale mass
      xav = 0.0
      do j= 1,  spvars%gnptot + spvars%gnbdy
         i = 1
           xav = xav + 0.5 * ( kern( part(i)%pos, part(j)%pos, part(i)%h ) &
                          + kern( part(i)%pos, part(j)%pos, part(j)%h ) )
      enddo

      do j= 1,  spvars%gnptot + spvars%gnbdy
          part(j)%mass = part(j)%mass / xav
      enddo


   if ( 1 == 0 ) then
      !correct mass
      !left mass
       asum = 0.0
       do j= 1,  spvars%gnptot + spvars%gnbdy
              i = 1
              xav = 0.5 * ( kern( part(i)%pos, part(j)%pos, part(i)%h ) &
                          + kern( part(i)%pos, part(j)%pos, part(j)%h ) )
              asum = asum + part(j)%mass*xav
       enddo
       print*, 'left mass correction', rhol / asum
       do i = 1, nl
         part(i)%mass = part(i)%mass * rhol /asum
       enddo
       do j = 1, nbl
          i  = j + spvars%gnptot + nbr
         part(i)%mass = part(i)%mass * rhol /asum
       enddo


      !right mass
       asum = 0.0
       do j= 1,  spvars%gnptot + spvars%gnbdy
              i = spvars%gnptot
              xav = 0.5 * ( kern( part(i)%pos, part(j)%pos, part(i)%h ) &
                          + kern( part(i)%pos, part(j)%pos, part(j)%h ) )
              asum = asum + part(j)%mass*xav
       enddo
       print*, 'right mass correction', rhor / asum
       do j = 1, nr
         i = (spvars%gnptot - nr) + j
         part(i)%mass = part(i)%mass * rhor /asum
       enddo
       do j = 1, nbr
         i  = j + spvars%gnptot
         part(i)%mass = part(i)%mass * rhor /asum
       enddo
   endif

 
    !init left boundary  density
        asum = 0.0
       do j= 1,  spvars%gnptot + spvars%gnbdy
            i = 1
              xav = 0.5 * ( kern( part(i)%pos, part(j)%pos, part(i)%h ) &
                          + kern( part(i)%pos, part(j)%pos, part(j)%h ) )
              asum = asum + part(j)%mass*xav
       enddo
       do j = 1, nbl
          i  = j + spvars%gnptot + nbr
         part(i)%den = asum
       enddo
  
    !init right boundary  density
        asum = 0.0
       do j= 1,  spvars%gnptot + spvars%gnbdy
             i  = spvars%gnptot
              xav = 0.5 * ( kern( part(i)%pos, part(j)%pos, part(i)%h ) &
                          + kern( part(i)%pos, part(j)%pos, part(j)%h ) )
              asum = asum + part(j)%mass*xav
       enddo
       do j = 1, nbr
         i  = j + spvars%gnptot
         part(i)%den = asum
       enddo
 

    
     spvars%agamma = agamma
     spvars%gxmin = xmin
     spvars%gxmax = xmax
     spvars%gymin = ymin
     spvars%gymax = ymax

   end subroutine


end module problem
