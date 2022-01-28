module artvismod
  use partmod
  use eos
  implicit none
 
  type artvistype


   real (kind=8) :: alpha
   real (kind=8) :: beta
   real (kind=8) :: eps

   integer :: opt = 1

   !opt = 1 ; M89 
   !opt = 2 ; MG83 

  endtype

  type (artvistype) :: artvis


   contains

     function artvis_p(parti, partj) result(oval)

          real (kind=8) :: oval
          type (parttype) :: parti, partj

          oval = 0.0

          if ( artvis%opt == 1 ) then
            oval =  artvis_p_m89(parti, partj)
          elseif ( artvis%opt == 2 ) then
            oval =  artvis_p_mg83(parti, partj)
          endif

     end function artvis_p

     !M89 viscosicy  
     function artvis_p_m89(parti, partj) result(oval)

          real (kind=8) :: oval
          type (parttype) :: parti, partj

          real (kind=8) :: csav, denav, hav
          real (kind=8) :: mu, vp, pvs, dpos(3), rsq



          oval = 0.0
          pvs =  dot_product( (parti%vel - partj%vel), (parti%pos - partj%pos ) )

          if ( pvs < 0.0 ) then

             csav = 0.5 * ( eos_cs(parti%en, parti%den, spvars%agamma) &
                            +  eos_cs(partj%en, partj%den, spvars%agamma)  &
                          )
             denav = (parti%den + partj%den)/2.0

             hav = 0.5*(parti%h + partj%h ) 

             dpos = parti%pos - partj%pos
             rsq = dot_product(dpos, dpos)

             mu = hav * pvs / (rsq + artvis%eps * hav**2.0 )

             vp = (-artvis%alpha *csav *mu  + artvis%beta * mu*mu )/denav
             oval = vp ;
          endif
     end function


     !MG83 viscosicy  
     function artvis_p_mg83(parti, partj) result(oval)
          real (kind=8) :: oval
          type (parttype) :: parti, partj

          real (kind=8) :: beta

          beta = artvis%beta ;
          artvis%beta = 0.0
          oval =  artvis_p_m89(parti, partj)
          artvis%beta = beta;

     end function




end module
