module types
use kinds
use discrete
use rotation
use translation
use adjacency
use alignment
implicit none
private

type :: Block
   integer :: blklen
   integer, allocatable :: blkidx(:)
end type

type, public :: Equiv
   integer :: eqvlen
   integer, allocatable :: eqvidx(:)
end type

type :: Atom
   character(:), allocatable :: label
   integer :: znum
   integer :: ztype
   real(wp) :: weight
   real(wp) :: coords(3)
   integer :: nadj
   integer, allocatable :: adjidx(:)
contains
   procedure :: print => print_atom
end type

type, public :: Molecule
   integer :: natom
   integer :: nblock
   integer :: nequiv
   character(:), allocatable :: title
   type(Atom), allocatable :: atoms(:)
   type(Block), allocatable :: blocks(:)
   type(Equiv), allocatable :: equivs(:)
contains
   procedure :: get_znums
   procedure :: get_ztypes
   procedure :: get_weights
   procedure :: get_coords
   procedure :: set_coords
   procedure :: get_labels
   procedure :: get_adjmat
   procedure :: mirror_coords
   procedure, private :: matrix_rotate_coords
   procedure, private :: quater_rotate_coords
   generic, public :: rotate_coords => matrix_rotate_coords, quater_rotate_coords
   procedure :: translate_coords
   procedure :: permutate_atoms
   procedure :: print => print_molecule
end type

contains

subroutine mirror_coords(self)

   class(Molecule), intent(inout) :: self
   real(wp) :: coords(3, self%natom)

   coords = self%get_coords()
   coords(1, :) = -coords(1, :)
   call self%set_coords(coords)

end subroutine mirror_coords

subroutine matrix_rotate_coords(self, rotmat)

   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: rotmat(3, 3)

   call self%set_coords(matrix_rotated(self%natom, self%get_coords(), rotmat))

end subroutine matrix_rotate_coords

subroutine quater_rotate_coords(self, rotquat)

   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: rotquat(4)

   call self%set_coords(quater_rotated(self%natom, self%get_coords(), rotquat))

end subroutine quater_rotate_coords

subroutine translate_coords(self, travec)

   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: travec(3)

   call self%set_coords(translated(self%natom, self%get_coords(), travec))

end subroutine translate_coords

subroutine permutate_atoms(self, order)

   class(Molecule), intent(inout) :: self
   integer, intent(in) :: order(:)
   integer i, k, invorder(self%natom)

   invorder = inverse_perm(order)
   self%atoms = self%atoms(order(:))
   do i = 1, self%natom
      do k = 1, self%atoms(i)%nadj
         self%atoms(i)%adjidx(k) = invorder(self%atoms(i)%adjidx(k))
      end do
   end do

end subroutine permutate_atoms

function get_znums(self) result(znums)

   class(Molecule), intent(in) :: self
   integer :: znums(self%natom)
   integer i

   do i = 1, self%natom
      znums(i) = self%atoms(i)%znum
   end do

end function get_znums

function get_ztypes(self) result(ztypes)

   class(Molecule), intent(in) :: self
   integer :: ztypes(self%natom)
   integer i

   do i = 1, self%natom
      ztypes(i) = self%atoms(i)%ztype
   end do

end function get_ztypes

function get_weights(self) result(weights)

   class(Molecule), intent(in) :: self
   real(wp) :: weights(self%natom)
   integer i

   do i = 1, self%natom
      weights(i) = self%atoms(i)%weight
   end do

end function get_weights

function get_coords(self) result(coords)

   class(Molecule), intent(in) :: self
   real(wp) :: coords(3, self%natom)
   integer i

   do i = 1, self%natom
      coords(:, i) = self%atoms(i)%coords
   end do

end function get_coords

subroutine set_coords(self, coords)

   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: coords(3, self%natom)
   integer i

   do i = 1, self%natom
      self%atoms(i)%coords = coords(:, i)
   end do

end subroutine set_coords

function get_labels(self) result(labels)

   class(Molecule), intent(in) :: self
   character(wl) :: labels(self%natom)
   integer i

   do i = 1, self%natom
      labels(i) = self%atoms(i)%label
   end do

end function get_labels

function get_adjmat(self) result(adjmat)

   class(Molecule), intent(in) :: self
   type(Atom) :: iatom
   logical :: adjmat(self%natom, self%natom)
   integer i, k

   adjmat(:, :) = .false.

   do i = 1, self%natom
      iatom = self%atoms(i)
      do k = 1, iatom%nadj
         adjmat(i, iatom%adjidx(k)) = .true.
      end do
   end do

end function get_adjmat

subroutine print_atom(self, ind, outLvl)
   class(Atom), intent(in) :: self
   integer, intent(in) :: ind, outLvl
   character(255) :: frmt
   character(2) :: num

! *** code to manage unitlbl pending ***
   
   select case (outLvl)

   case (1)
      if (self%nadj == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (*, frmt) ind, ": ", self%label(:2), self%znum, self%ztype, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') self%nadj
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (*, frmt) ind, ": ", self%label(:2), self%znum, self%ztype, &
               self%weight, " {", self%coords(:), " } [", &
               self%adjidx(:self%nadj), " ]"
      end if

!  case (2)
!
   case default
      if (self%nadj == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (*, frmt) ind, ": ", self%label(:2), self%znum, self%ztype, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') self%nadj
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (*, frmt) ind, ": ", self%label(:2), self%znum, self%ztype, &
               self%weight, " {", self%coords(:), " } [", &
               self%adjidx(:self%nadj), " ]"
      end if
   end select

end subroutine print_atom

subroutine print_molecule(self, molname)
   class(Molecule), intent(in) :: self
   character(*), intent(in) :: molname
   integer :: i

! *** code to manage unitlbl pending ***
   
   write (*, '(a,1x,a,3x,a,i0,a)') "Contents of molecule structure:", molname, &
                                   "(", self%natom, " atoms)"
   write (*, '(a,a4,a5,a6,a7,2a17)') "ind:", "lbl", "znum", "ztype", "weight", &
                     "{ coords }", "[ adjlist ]"

   do i = 1, self%natom
      call self%atoms(i)%print(i, 1)
   end do

end subroutine print_molecule

end module
