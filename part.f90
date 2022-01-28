module partmod
   implicit none

    integer, parameter :: t_maxnb = 100

  type parttype
     real (kind=8) :: pos(3)
     real (kind=8) :: vel(3)
     real (kind=8) :: en
     real (kind=8) :: p
     real (kind=8) :: den
     real (kind=8) :: h
     real (kind=8) :: mass
     integer       :: nb
     integer       :: ptype
     integer       :: nblist(t_maxnb)

     real (kind=8) :: delpos(3)
     real (kind=8) :: accel(3)
     real (kind=8) :: dele
     real (kind=8) :: delden
     real (kind=8) :: delh

  endtype

   type (parttype) ,allocatable :: part(:)

  type spvartype
   integer  :: giter
   integer  :: gitermax 
   integer  :: gnbdymax !maximum bdy particles
   integer  :: gnptot   !live particles 

   integer  :: gnbdy    !used bdy particles
   integer  :: xlng, xrng
   integer  :: ylng, yrng
   integer  :: zlng, zrng

   integer  :: kopt

   integer  :: ndim

   integer  :: pmethod
   integer  :: isentropy
   integer  :: isartvis
   integer  :: isproblemdt
   integer  :: dtmethod = 1

   integer       :: maxnb = t_maxnb

   real (kind=8) :: alpha, beta
   real (kind=8) :: gtmax
   real (kind=8) :: gtime = 0.0
   real (kind=8) :: gdt
   real (kind=8) :: agamma
   real (kind=8) :: gxmin, gxmax
   real (kind=8) :: gymin, gymax
   real (kind=8) :: gzmin, gzmax
  endtype

   type (spvartype) :: spvars

  contains

    subroutine initvars ()
       integer :: astat

       allocate(part(spvars%gnptot+spvars%gnbdymax), stat= astat )
       if ( astat /= 0 ) stop "Memory allocation error for: part"

    end subroutine 
 
end module partmod
