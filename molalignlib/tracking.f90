! purpose: determines the most convenient atom to start runs over the structure
!          of each fragment in a molecule. (returns nfrag, fragidcs)
module tracking

use stdio
use molecule
use bounds
use sorting

implicit none
private
public find_molfrags

type(atomlist_type), allocatable :: adjlists(:)

contains

subroutine find_molfrags (mol)
   type(mol_type), intent(inout) :: mol
   ! Local variables
   integer :: iatom, nfrag
   integer, allocatable :: fragidcs(:), fragsize(:)
   logical, allocatable :: tracked(:)

   allocate (fragidcs(size(mol%atoms)))
   allocate (fragsize(size(mol%atoms)))
   allocate (tracked(size(mol%atoms)))

   adjlists = mol%get_adjlists()

   ! initialization

   nfrag = 0
   tracked(:) = .false.
   fragidcs(:) = 0
   fragsize(:) = 0

   ! detect fragments and populate frag arrays
   iatom = 1
   do while (iatom <= size(mol%atoms))
      if (tracked(iatom)) then
         iatom = iatom + 1
      else
         nfrag = nfrag + 1
         call recrun(tracked, iatom, nfrag, fragidcs, fragsize)
         iatom = 1
      end if
   end do

   ! register nfrag and fragidcs in the mol strucutre
   call mol%set_molfrags(nfrag, fragidcs)

end subroutine

recursive subroutine recrun (tracked, iatom, nfrag, fragidcs, fragsize)
! runs recursivelly over the structure and populates arrays
   logical, intent(inout) :: tracked(:)
   integer, intent(in) :: iatom, nfrag
   integer, intent(inout) :: fragidcs(:), fragsize(:)
   ! Local variables
   integer :: i

   if (tracked(iatom)) return

   tracked(iatom) = .true.
   fragidcs(iatom) = nfrag
   fragsize(nfrag) = fragsize(nfrag) + 1
   
   do i = 1, size(adjlists(iatom)%atomidcs)
      call recrun(tracked, adjlists(iatom)%atomidcs(i), nfrag, fragidcs, fragsize)
   end do

end subroutine

subroutine set_molfrags(self, nfrag, fragidcs)
   class(mol_type), intent(inout) :: self
   integer, intent(in) :: nfrag
   integer, intent(in) :: fragidcs(:)

   self%molfrags = atompartition(nfrag, fragidcs)

end subroutine

function get_molfrags(self) result(molfrags)
   class(mol_type), intent(in) :: self
   ! Local variables
   integer :: i
   type(atomlist_type), allocatable :: molfrags(:)

   allocate (molfrags(size(self%molfrags)))

   do i = 1, size(self%molfrags)
      molfrags(i)%atomidcs = self%molfrags(i)%atomidcs
   end do

end function

function get_molfragroots(self) result(fragroots)
   class(mol_type), intent(in) :: self
   ! Result variable
   integer, allocatable :: fragroots(:)
   ! Local variables
   integer :: h, i
   integer :: iatom, root_atom
   integer :: eltypepop_i, eltypepop_min
   integer, allocatable :: atomeltypes(:)

   allocate (fragroots(size(self%molfrags)))

   atomeltypes = self%eltypes%indices
   do h = 1, size(self%molfrags)
      eltypepop_min = huge(eltypepop_min)
      do i = 1, size(self%molfrags(h)%atomidcs)
         iatom = self%molfrags(h)%atomidcs(i)
         eltypepop_i = size(self%eltypes%parts(atomeltypes(iatom))%items)
         if (eltypepop_i < eltypepop_min) then
            root_atom = iatom
            eltypepop_min = eltypepop_i
         end if
         fragroots(h) = root_atom
      end do
   end do

end function

end module
