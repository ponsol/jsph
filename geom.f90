
module geom
   use partmod
   use problem
   use artvismod
   implicit none

    contains

    subroutine initgeom

       spvars%ndim = 3
       spvars%gnptot = 100
       spvars%gnbdy = 0
       spvars%gitermax = 1
       spvars%isartvis = 0
       spvars%isproblemdt = 1
       spvars%isentropy = 0
       spvars%pmethod = 1

       spvars%xlng = 10
       spvars%xrng = spvars%xlng;

       spvars%ylng = spvars%xlng;
       spvars%yrng = spvars%ylng;

       spvars%zlng = spvars%ylng;
       spvars%zrng = spvars%zlng;


       artvis%alpha = 0.5
       artvis%beta =  0.5 
       artvis%eps = 0.01

     call problemgeom()
     call initvars()

       part(1:spvars%gnptot)%ptype = 1
       part(spvars%gnptot+1:spvars%gnptot+spvars%gnbdy)%ptype = 2

    end subroutine initgeom

end module geom

