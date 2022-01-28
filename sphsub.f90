module sphsub
   use partmod
   use artvismod
   use kernel
   use problem
   implicit none

   contains
 

     subroutine getdt(dt)
      real (kind=8), intent(inout) :: dt

       if ( spvars%isproblemdt == 1 ) then
         call problemgetdt(dt) 
       elseif ( spvars%dtmethod == 1 ) then
         call getdt_classic(dt)                
       elseif ( spvars%dtmethod == 2 ) then
         call getdt_m2(dt)                
       endif

     end subroutine

     function getdelv(pos, ivel) result(oval)
       real (kind=8 ) :: pos(3), ivel(3)
       real (kind=8 ) :: oval
       real (kind=8 ) :: tval
       integer :: i

       oval = 0.0
       do i = 1, spvars%gnptot + spvars%gnbdy
        tval = dot_product(  ivel - part(i)%vel, gradkern(pos, part(i)%pos, part(i)%h ) )
        oval = oval +  part(i)%mass / part(i)%den * tval
       enddo

     end function 

 
     subroutine getdt_classic(dt)
      real (kind=8), intent(inout) :: dt
      real (kind=8) :: h, cs
      real (kind=8) :: dt1, dt1s
      real (kind=8) :: delv, avcs
      integer :: i

         do i= 1, spvars%gnptot 
            h = part(i)%h
            cs = eos_cs(part(i)%p, part(i)%den, spvars%agamma)

            if ( artvis%opt == 1 ) then
             delv = getdelv(part(i)%pos, part(i)%vel )
             avcs = h*delv + 1.2*artvis%alpha * cs

             if ( delv < 0.0 ) avcs = avcs +  1.2*artvis%beta*h *delv 

             cs = cs + avcs
            endif

            dt1 = h / sqrt(cs)

            if ( i == 1 ) dt1s = dt1
            if ( dt1 < dt1s ) dt1s = dt1

         enddo
 
            dt = 0.3 * dt1

     end subroutine
 
 
     subroutine getdt_m2(dt)
      real (kind=8), intent(inout) :: dt
      real (kind=8) :: h, cs, vel(3)
      real (kind=8) :: m1, m2, m1s, m2s
      real (kind=8) :: dt1, dt2, dt1s, dt2s
      integer :: i

         do i= 1, spvars%gnptot 
            h = part(i)%h
            vel = part(i)%vel
            cs = eos_cs(part(i)%p, part(i)%den, spvars%agamma)
            m1  = dot_product(vel,vel) + cs*cs

            if ( i == 1 ) m1s = m1
            if ( m1 > m1s ) m1s = m1

            dt1 = h / sqrt(m1)
            if ( i == 1 ) dt1s = dt1
            if ( dt1 < dt1s ) dt1s = dt1


            m2 = dot_product( part(i)%accel, part(i)%accel )
            if ( i == 1 ) m2s = m2
            if ( m2 > m2s ) m2s = m2

            dt2 = sqrt( h / sqrt(m2) )
            if ( i == 1 ) dt2s = dt2
            if ( dt2 < dt2s ) dt2s = dt2

         enddo
 
            dt = 0.3 * min ( h/sqrt(m1s),  sqrt( h / sqrt(m2s) ) )
            !dt = 0.3 * min ( dt1, dt2)

     end subroutine
 

     subroutine getdensity()
        real (kind=8) :: asum, x1, xav
        integer :: i, j, jj
        real (kind=8) :: asum1

         !calculate density
         do i= 1, spvars%gnptot 
           asum  = 0.0
           asum1  = 0.0

!           do jj= 1, part(i)%nb 
!                   j=  part(i)%nblist(jj) 
           do j= 1,  spvars%gnptot + spvars%gnbdy

              xav = 0.5 * ( kern( part(i)%pos, part(j)%pos, part(i)%h ) &
                          + kern( part(i)%pos, part(j)%pos, part(j)%h ) )
              asum1 = asum1 + xav
              asum = asum + part(j)%mass*xav
           enddo 
            part(i)%den = asum
         enddo
 
     end subroutine

     !calculate en variable from pressure
     subroutine get_en_from_p()
        real (kind=8) :: asum
        real (kind=8) :: asum1
        integer :: i

         do i= 1,  spvars%gnptot + spvars%gnbdy

            part(i)%en = eos_e_prho(spvars%isentropy, part(i)%p , part(i)%den , spvars%agamma )
         enddo

     end subroutine

  
     subroutine initfinal()
        write(*,*) '#total number of particles', spvars%gnptot+spvars%gnbdy
        artvis%alpha = spvars%alpha
        artvis%beta = spvars%beta
        artvis%opt = spvars%isartvis
     end subroutine


     subroutine fillvars()

        call getdensity()
        call get_en_from_p()

     end subroutine

 

     subroutine get_dden()
         real (kind=8) :: asum

         real (kind=8) :: kvec1(3), kvec2(3), kvec(3), vec(3), tval 

         integer :: i, j, jj

         do i= 1, spvars%gnptot 

              asum = 0.0
             !do jj= 1, part(i)%nb 
             !      j=  part(i)%nblist(jj) 
             do j= 1, spvars%gnptot +spvars%gnbdy

                 kvec1 = gradkern( part(i)%pos, part(j)%pos, part(i)%h )
                 kvec2 = gradkern( part(i)%pos, part(j)%pos, part(j)%h )

                 kvec = (kvec1+kvec2)/2.0
                 vec =   part(i)%vel - part(j)%vel

                !print*, 'in den 0',  kvec
                !print*, 'in den 00',  vec
                tval = dot_product(kvec, vec)
                !print*, 'in den 1',  tval
                tval = tval *  part(j)%mass 

              asum = asum + tval
            enddo

              part(i)%delden = asum

         enddo 

     end subroutine

     subroutine get_dvel()
         integer :: i, j, jj
         real (kind=8) :: asum(3)
         real (kind=8) ::  kvec1(3), kvec2(3), s1, s2, pt,  ap

            
         do i= 1, spvars%gnptot 
               asum = 0.0
               ap = 0.0
             !do jj= 1, part(i)%nb 
             !      j=  part(i)%nblist(jj) 
             do j= 1, spvars%gnptot +spvars%gnbdy

                 kvec1 = gradkern( part(i)%pos, part(j)%pos, part(i)%h )
                 kvec2 = gradkern( part(i)%pos, part(j)%pos, part(j)%h )

               s1 =   part(i)%p / part(i)%den**2.0
               s2 =   part(j)%p / part(j)%den**2.0
               pt = s1 + s2
               if ( spvars%isartvis > 0 ) then
                 ap =  artvis_p( part(i), part(j) )
                 pt = pt + ap
               endif
                
               asum = asum - part(j)%mass * pt * (kvec1 + kvec2) / 2.0d0 

             enddo
                

                part(i)%accel = asum

         enddo 


     end subroutine



     subroutine get_de()
         integer :: i, j, jj
         real (kind=8) :: asum, tval
         real (kind=8) :: kvec1(3), kvec2(3), vec(3)
         real (kind=8) :: s1, s2, ap, agamma
 

         agamma = spvars%agamma

         do i= 1, spvars%gnptot 

                     asum = 0.0
                     ap = 0.0
             !      do jj= 1, part(i)%nb 
             !          j=  part(i)%nblist(jj) 
             do j= 1, spvars%gnptot +spvars%gnbdy

                        ap = 0.0
                      if ( spvars%isartvis > 0 ) then
                         ap =  artvis_p( part(i), part(j) )
                      endif

                      if ( spvars%isentropy == 1 ) then

                         vec  =  part(i)%vel - part(j)%vel
                         kvec1 = gradkern( part(i)%pos, part(j)%pos, part(i)%h )
                         kvec2 = gradkern( part(i)%pos, part(j)%pos, part(j)%h )

                         tval = 0.5 * part(j)%mass * ap * dot_product(vec, (kvec1+kvec2)/2.0 )
                         tval = tval* (agamma -1.0)/part(i)%den**(agamma -1.0)

                         asum = asum + tval

                      elseif ( spvars%pmethod == 1 ) then

                           vec  =  part(i)%vel - part(j)%vel 
                           kvec1 = gradkern( part(i)%pos, part(j)%pos, part(i)%h )
                           kvec2= gradkern( part(i)%pos, part(j)%pos, part(i)%h )

                           s1 = dot_product(vec, kvec1) * ( part(i)%p / part(i)%den**2.0 + ap/2.0)

                           s2 = dot_product(vec, kvec2) * ( part(j)%p / part(j)%den**2.0 + ap/2.0 )

                           asum = asum + (s1 + s2)*0.5*part(j)%mass

                      endif !entropy

                   enddo !j particles

                 part(i)%dele = asum

         enddo 


     end subroutine


     subroutine allupdate()
         integer :: i
         real (kind=8) :: dt

          dt = spvars%gdt  
         do i= 1, spvars%gnptot 


           part(i)%vel = part(i)%vel  + dt* part(i)%accel

           part(i)%pos = part(i)%pos  + dt* part(i)%vel  

           part(i)%en = part(i)%en  + dt* part(i)%dele

           call getdensity()
           part(i)%p =   eos_p_erho(spvars%isentropy, part(i)%en , part(i)%den, spvars%agamma)

           part(i)%delden = 0.0
           part(i)%accel = 0.0
           part(i)%delpos = 0.0
           part(i)%dele = 0.0

         enddo

     end subroutine


     subroutine boundary()
        integer :: i, j

         return

         do j= 1, spvars%gnbdy/2
             i = spvars%gnptot+ spvars%gnbdy/2 + j
             part(i)%p = part(spvars%gnptot)%p
             part(i)%vel = part(spvars%gnptot)%vel
             part(i)%en = part(spvars%gnptot)%en
             part(i)%den = part(spvars%gnptot)%den
         enddo
 
     end subroutine

   function isboundary(pos,h) result(oval)
      integer :: oval(3)
      real (kind=8) :: pos(3), h
      real (kind=8) :: xp, yp, zp
      real (kind=8) :: xbmin, xbmax, ybmin, ybmax, zbmin, zbmax

          oval = 0

          xp = pos(1)
          yp = pos(2)
          zp = pos(3)

          xbmin = spvars%gxmin
          ybmin = spvars%gymin
          zbmin = spvars%gzmin

          xbmax = spvars%gxmax
          ybmax = spvars%gymax
          zbmax = spvars%gzmax

         if ( xp - h < xbmin ) oval(1) = -1
         if ( xp + h > xbmax ) oval(1) =  1
         if ( yp - h < ybmin ) oval(2) = -1
         if ( yp + h > ybmax ) oval(2) =  1
         if ( zp - h < zbmin ) oval(3) = -1
         if ( zp + h > zbmax ) oval(3) =  1

   end function

   subroutine getnblist()
     integer :: i, j
     integer :: m
     real (kind=8) :: x, y, z, r, ha
     real (kind=8) :: xt, yt, zt


     do i = 1, spvars%gnptot 

       m=0
       do j = 1, spvars%gnptot + spvars%gnbdy
          x= part(i)%pos(1)- part(j)%pos(1)
          y= part(i)%pos(2)- part(j)%pos(2)
          z= part(i)%pos(3)- part(j)%pos(3)
          
          r = sqrt( x*x + y*y + z*z )

          ha = part(i)%h 

           if ( r <= ha ) then
              m = m+1
              part(i)%nblist(m) = j
              !print*, 'type ', part(j)%ptype, part(i)%pos(1), part(j)%pos(1)
           endif
            if ( m == spvars%maxnb ) exit

       enddo !j
       part(i)%nb =  m

     enddo !i
   end subroutine 

   subroutine outonecut(dir, n, m)
      integer :: dir, n, m
      integer :: ncut, i, j

      real ( kind=8 ) :: amin , amax, da, a
      real ( kind=8 ) :: rho, p, vel(3)
      real ( kind=8 ) :: tp

      real ( kind=8 ) :: pos(3)
      real ( kind=8 ) :: kav, sumk
      integer :: ios
      integer,save :: f11open = 0

      ncut = 100


      if ( f11open == 0 ) then
        open(11, file='pcut', status='replace', iostat=ios)
        f11open = 1
        if ( ios /= 0 ) return
      else
        open(11, file='pcut', status='old', position='append', iostat=ios)
        if ( ios /= 0 ) return
      endif

      write(11,*) '# iter ', spvars%giter,  'time ', spvars%gtime
      write(11,*) '# x ',  ' den   ', ' p ', '  v1 '

      amin = spvars%gxmin; amax = spvars%gxmax; 
      da = (amax - amin)/(ncut -1 )
      do  i=1, ncut
         a = amin + (i-1)*da
         rho=0.0; p=0.0; vel = 0.0; sumk = 0.0;

         pos(1) = a; pos(2) = 0.0; pos(3) = 0.0;

         do j = 1, spvars%gnptot+spvars%gnbdy 
            kav =  kern(pos, part(j)%pos, part(j)%h )
            sumk = sumk + kav
            rho = rho + part(j)%mass*kav
            p = p +  part(j)%p*part(j)%mass*kav /  part(j)%den
            vel = vel +  part(j)%vel*part(j)%mass*kav/ part(j)%den
         enddo

         !write(11,*) a, rho/sumk, p/sumk, vel(1)/sumk
         write(11,*) a, rho, p, vel(1)
      enddo
      write(11,*) ''
      write(11,*) ''
      close(11)

   end subroutine 

   subroutine outallpart()
      integer :: i

      real ( kind=8 ) :: x, y, z
      real ( kind=8 ) :: rho, p, v1, v2, v3
      real ( kind=8 ) :: tp
      integer :: ios
      integer, save :: f10open = 0

      if ( f10open == 0 ) then
          open(10, file='pall', status='replace', iostat=ios)
          f10open = 1
          if ( ios /= 0 ) return
      else
          open(10, file='pall', status='old', position='append', iostat=ios)
          if ( ios /= 0 ) return
      endif

      write(10,*) '# iter ', spvars%giter,  'time ', spvars%gtime
      write(10,*) '#x, rho, p, v1'
    
      do  i=1, spvars%gnptot + spvars%gnbdy
            rho = part(i)%den
            p =  part(i)%p
            v1 = part(i)%vel(1)
            v2 = part(i)%vel(2)
            v3 = part(i)%vel(3)

            x = part(i)%pos(1)
            y = part(i)%pos(2)
            z = part(i)%pos(3)
            !write(10,*) x, y, rho, p, v1
            write(10,*) x, rho, p, v1
      enddo
      write(10,*) ''
      write(10,*) ''
      close(10)

   end subroutine 

end module sphsub

