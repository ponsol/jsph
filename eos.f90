
module eos
   implicit none

  contains

    function eos_cs(p, rho, agamma) result(oval)
         real (kind=8) :: p, rho, agamma
         real (kind=8) :: oval
         oval = eos_cs_prho(p, rho, agamma)

    end function


    function eos_cs_prho(p, rho, agamma) result(oval)
         real (kind=8) :: p, rho, agamma
         real (kind=8) :: oval

         oval = sqrt ( agamma* p/rho )
    end function

    function eos_cs_erho(isen, en, rho, agamma) result(oval)
         real (kind=8) :: en, rho, agamma
         real (kind=8) :: oval
         real (kind=8) :: p
         integer :: isen

         p = eos_p_erho(isen, en, rho, agamma)
         oval =  eos_cs_prho(p, rho, agamma)

    end function

    function eos_p_erho(isen, en, rho, agamma) result(oval)
         real (kind=8) :: en, rho, agamma
         real (kind=8) :: oval
         integer :: isen

        if ( isen == 1 ) then
             oval  =   eos_p_erho_isen(en , rho, agamma)
        else
            oval = (agamma -1.0) * en * rho
        endif

    end function


    function eos_p_erho_isen(en, rho, agamma) result(oval)
         real (kind=8) :: en, rho, agamma
         real (kind=8) :: oval
         oval =  en * rho ** agamma
    end function



    function eos_e_prho_isen(p, rho, agamma) result(oval)
         real (kind=8) :: p, rho, agamma
         real (kind=8) :: oval
          oval =  p / rho**agamma
    end function


    function eos_e_prho(isen, p, rho, agamma) result(oval)
         real (kind=8) :: p, rho, agamma
         real (kind=8) :: oval
         integer :: isen

         if ( isen == 1 ) then
          oval = eos_e_prho_isen(p, rho, agamma)
         else
          oval =  p / (agamma -1.0) /  rho
         endif

    end function





end module eos
