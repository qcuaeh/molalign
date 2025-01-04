! purpose: determines the most convenient atom to start runs over the structure
!          of each fragment in a molecule. (returns nfrag, fragidcs)
module tracking
use parameters
use sorting
use molecule
use partition
use common_types

implicit none
private
public find_molfrags

type(atom_type), allocatable :: atoms(:)

contains

subroutine find_molfrags( mol, eltypes, molfrags)
   type(mol_type), intent(in) :: mol
   type(partition_type), intent(in) :: eltypes
   type(intlist_type), allocatable, intent(out) :: molfrags(:)
   ! Local variables
   integer :: i, nfrag
   logical, allocatable :: tracked(:)
   integer, allocatable :: fragszs(:), fragidcs(:,:)
   integer, allocatable :: order(:)

   allocate (fragszs(size(mol%atoms)))
   allocate (fragidcs(size(mol%atoms), size(mol%atoms)))
   allocate (tracked(size(mol%atoms)))
   atoms = mol%atoms

   ! initialization

   nfrag = 0
   fragszs(:) = 0
   tracked(:) = .false.

   ! detect fragments and populate frag arrays
   i = 1
   do while (i <= size(mol%atoms))
      if (tracked(i)) then
         i = i + 1
      else
         nfrag = nfrag + 1
         call recrun( tracked, i, nfrag, fragszs, fragidcs)
         i = 1
      end if
   end do

   ! Order molecular fragments
   do i = 1, nfrag
      order = sorted_order(eltypes%parts(eltypes%idcs(fragidcs(:fragszs(i), i)))%part_size)
      fragidcs(:fragszs(i), i) = fragidcs(order, i)
!      write (stderr, *) fragidcs1(:fragszs1(i), i)
!      write (stderr, *)
   end do

   allocate (molfrags(nfrag))

   do i = 1, nfrag
      molfrags(i)%n = fragidcs(:fragszs(i), i)
   end do

end subroutine

recursive subroutine recrun( tracked, iatom, nfrag, fragszs, fragidcs)
! runs recursivelly over the structure and populates arrays
   logical, intent(inout) :: tracked(:)
   integer, intent(in) :: iatom, nfrag
   integer, intent(inout) :: fragidcs(:,:), fragszs(:)
   ! Local variables
   integer :: i

   if (tracked(iatom)) return

   tracked(iatom) = .true.
   fragszs(nfrag) = fragszs(nfrag) + 1
   fragidcs(fragszs(nfrag), nfrag) = iatom
   
   do i = 1, size(atoms(iatom)%adjlist)
      call recrun( tracked, atoms(iatom)%adjlist(i), nfrag, fragszs, fragidcs)
   end do

end subroutine

!subroutine set_molfrags(self, nfrag, fragidcs)
!   class(mol_type), intent(inout) :: self
!   integer, intent(in) :: nfrag
!   integer, intent(in) :: fragidcs(:)
!
!   self%molfrags = atompartition(nfrag, fragidcs)
!
!end subroutine
!
!function get_molfrags(self) result(molfrags)
!   class(mol_type), intent(in) :: self
!   ! Local variables
!   integer :: i
!   type(atomlist_type), allocatable :: molfrags(:)
!
!   allocate (molfrags(size(self%molfrags)))
!
!   do i = 1, size(self%molfrags)
!      molfrags(i)%atomidcs = self%molfrags(i)%atomidcs
!   end do
!
!end function
!
!function get_molfragroots(self) result(fragroots)
!   class(mol_type), intent(in) :: self
!   ! Result variable
!   integer, allocatable :: fragroots(:)
!   ! Local variables
!   integer :: h, i
!   integer :: iatom, root_atom
!   integer :: eltypepop_i, eltypepop_min
!   integer, allocatable :: atomeltypes(:)
!
!   allocate (fragroots(size(self%molfrags)))
!
!   atomeltypes = self%eltypes%idcs
!   do h = 1, size(self%molfrags)
!      eltypepop_min = huge(eltypepop_min)
!      do i = 1, size(self%molfrags(h)%atomidcs)
!         iatom = self%molfrags(h)%atomidcs(i)
!         eltypepop_i = size(self%eltypes%parts(atomeltypes(iatom))%items)
!         if (eltypepop_i < eltypepop_min) then
!            root_atom = iatom
!            eltypepop_min = eltypepop_i
!         end if
!         fragroots(h) = root_atom
!      end do
!   end do
!
!end function

end module
