module types
use kinds
use bounds
use discrete
use rotation
use translation
use adjacency
use alignment
use strutils
implicit none
private

type, public :: Part
   integer, allocatable :: atomidcs(:)
end type

type, public :: Partition
   type(Part), allocatable :: parts(:)
contains
   procedure :: init => partition_init
   procedure :: get_lenlist => partition_get_lenlist
end type

type :: Atom
   integer, private :: elnum
   integer, private :: label
   integer, private :: eltype
   integer, private :: mnatype
   integer, private :: fragidx
   integer, allocatable, private :: mnaid(:)
   real(wp), private :: weight
   real(wp), private :: coords(3)
   integer, allocatable, private :: adjlist(:)
   integer, allocatable :: atomneimnatypepartlens(:)
contains
   procedure :: print => print_atom
end type

type, public :: Molecule
   integer :: natom
   integer :: nfrag
   character(:), allocatable :: title
   type(Atom), allocatable :: atoms(:)
   type(Partition) :: eltypepartition
   type(Partition) :: mnatypepartition
   type(Partition) :: molfragpartition
contains
   procedure :: print => print_molecule
   procedure :: set_elnums
   procedure :: set_labels
   procedure :: set_coords
   procedure :: set_weights
   procedure :: set_adjlists
   procedure :: set_eltypes
   procedure :: set_mnatypes
   procedure :: set_molfrags
   procedure :: set_atomneimnatypepartlens
   procedure :: get_natom
   procedure :: get_neltype
   procedure :: get_nmnatype
   procedure :: get_eltypepartlens
   procedure :: get_mnatypepartlens
   procedure :: mirror_coords
   procedure :: translate_coords
   procedure :: permutate_atoms
   procedure :: bonded
   procedure :: remove_bond
   procedure :: add_bond
   procedure :: get_fragparts
   procedure :: get_fragroots
   procedure :: get_eltypeparts
   procedure, private :: get_default_nadjs
   procedure, private :: get_default_adjlists
   procedure, private :: get_default_adjmat
   procedure, private :: get_default_coords
   procedure, private :: get_default_elnums
   procedure, private :: get_default_atommnatypes
   procedure, private :: get_default_labels
   procedure, private :: get_default_atomeltypes
   procedure, private :: get_default_atomweights
   procedure, private :: get_default_natomneimnatypes
   procedure, private :: get_default_atomneimnatypepartlens
   procedure, private :: get_sorted_nadjs
   procedure, private :: get_sorted_adjlists
   procedure, private :: get_sorted_adjmat
   procedure, private :: get_sorted_coords
   procedure, private :: get_sorted_elnums
   procedure, private :: get_sorted_atommnatypes
   procedure, private :: get_sorted_labels
   procedure, private :: get_sorted_atomeltypes
   procedure, private :: get_sorted_atomweights
   procedure, private :: get_sorted_natomneimnatypes
   procedure, private :: get_sorted_atomneimnatypepartlens
   procedure, private :: matrix_rotate_coords
   procedure, private :: quater_rotate_coords
   generic :: get_nadjs => get_default_nadjs, get_sorted_nadjs
   generic :: get_natomneimnatypes => get_default_natomneimnatypes, get_sorted_natomneimnatypes
   generic :: get_atomeltypes => get_default_atomeltypes, get_sorted_atomeltypes
   generic :: get_atommnatypes => get_default_atommnatypes, get_sorted_atommnatypes
   generic :: get_atomneimnatypepartlens => get_default_atomneimnatypepartlens, get_sorted_atomneimnatypepartlens
   generic :: get_coords => get_default_coords, get_sorted_coords
   generic :: get_adjlists => get_default_adjlists, get_sorted_adjlists
   generic :: get_adjmat => get_default_adjmat, get_sorted_adjmat
   generic :: get_elnums => get_default_elnums, get_sorted_elnums
   generic :: get_labels => get_default_labels, get_sorted_labels
   generic :: get_atomweights => get_default_atomweights, get_sorted_atomweights
   generic :: rotate_coords => matrix_rotate_coords, quater_rotate_coords
end type

contains

integer function get_natom(self) result(natom)
   class(Molecule), intent(in) :: self

   natom = size(self%atoms)

end function

integer function get_neltype(self) result(neltype)
   class(Molecule), intent(in) :: self

   neltype = size(self%eltypepartition%parts)

end function

integer function get_nmnatype(self) result(nmnatype)
   class(Molecule), intent(in) :: self

   nmnatype = size(self%mnatypepartition%parts)

end function

subroutine mirror_coords(self)
   class(Molecule), intent(inout) :: self
   ! Local variables
   real(wp), allocatable :: coords(:, :)

   allocate(coords(3, self%natom))

   coords = self%get_coords()
   coords(1, :) = -coords(1, :)
   call self%set_coords(coords)

end subroutine

subroutine matrix_rotate_coords(self, rotmat)
   class(Molecule), intent(inout) :: self
   ! Local variables
   real(wp), intent(in) :: rotmat(3, 3)

   call self%set_coords(matrix_rotated(self%natom, self%get_coords(), rotmat))

end subroutine

subroutine quater_rotate_coords(self, rotquat)
   class(Molecule), intent(inout) :: self
   ! Local variables
   real(wp), intent(in) :: rotquat(4)

   call self%set_coords(quater_rotated(self%natom, self%get_coords(), rotquat))

end subroutine

subroutine translate_coords(self, travec)
   class(Molecule), intent(inout) :: self
   ! Local variables
   real(wp), intent(in) :: travec(3)

   call self%set_coords(translated(self%natom, self%get_coords(), travec))

end subroutine

subroutine permutate_atoms(self, atomorder)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i, k
   integer, allocatable :: atomunorder(:)

   allocate(atomunorder(self%natom))

   atomunorder = inverse_mapping(atomorder)
   self%atoms = self%atoms(atomorder(:))
   do i = 1, self%natom
      do k = 1, size(self%atoms(i)%adjlist)
         self%atoms(i)%adjlist(k) = atomunorder(self%atoms(i)%adjlist(k))
      end do
   end do

end subroutine

subroutine set_elnums(self, atomelnums)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomelnums(:)

   self%atoms(:)%elnum = atomelnums(:)

end subroutine

function get_default_elnums(self) result(atomelnums)
   class(Molecule), intent(in) :: self
   integer, allocatable :: atomelnums(:)

   atomelnums = self%atoms(:)%elnum

end function

function get_sorted_elnums(self, atomorder) result(atomelnums)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer, allocatable :: atomelnums(:)

   atomelnums = self%atoms(atomorder)%elnum

end function

subroutine set_labels(self, atomlabels)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomlabels(:)

   self%atoms(:)%label = atomlabels(:)

end subroutine

function get_default_labels(self) result(atomlabels)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer, allocatable :: atomlabels(:)

   atomlabels = self%atoms(:)%label

end function

function get_sorted_labels(self, atomorder) result(atomlabels)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer, allocatable :: atomlabels(:)

   atomlabels = self%atoms(atomorder)%label

end function

subroutine set_eltypes(self, neltype, atomeltypes)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: neltype
   integer, intent(in) :: atomeltypes(:)

   self%atoms(:)%eltype = atomeltypes(:)
   call self%eltypepartition%init(neltype, atomeltypes)

end subroutine

subroutine set_mnatypes(self, nmnatype, atommnatypes)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: nmnatype
   integer, intent(in) :: atommnatypes(:)

   self%atoms(:)%mnatype = atommnatypes(:)
   call self%mnatypepartition%init(nmnatype, atommnatypes)

end subroutine

subroutine set_molfrags(self, nfrag, fragidcs)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: nfrag
   integer, intent(in) :: fragidcs(:)

   self%atoms(:)%fragidx = fragidcs(:)
   call self%molfragpartition%init(nfrag, fragidcs)

end subroutine

subroutine partition_init(self, npart, partidcs)
   class(Partition), intent(inout) :: self
   integer, intent(in) :: npart
   integer, intent(in) :: partidcs(:)
   ! Local variables
   integer :: i
   integer, allocatable :: n(:)

   allocate(n(npart))
   allocate(self%parts(npart))

   n(:) = 0
   do i = 1, size(partidcs)
      n(partidcs(i)) = n(partidcs(i)) + 1
   end do

   do i = 1, npart
      allocate(self%parts(i)%atomidcs(n(i)))
   end do

   n(:) = 0
   do i = 1, size(partidcs)
      n(partidcs(i)) = n(partidcs(i)) + 1
      self%parts(partidcs(i))%atomidcs(n(partidcs(i))) = i
   end do

end subroutine

function partition_get_lenlist(self) result(lenlist)
   class(Partition), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: lenlist(:)

   allocate(lenlist(size(self%parts)))

   do i = 1, size(self%parts)
      lenlist(i) = size(self%parts(i)%atomidcs)
   end do

end function

function get_mnatypepartlens(self) result(mnatypepartlens)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer, allocatable :: mnatypepartlens(:)

   mnatypepartlens = self%mnatypepartition%get_lenlist()

end function

function get_eltypeparts(self, atomorder) result(parts)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   integer, allocatable :: atomunorder(:)
   type(Part), allocatable :: parts(:)

   atomunorder = inverse_mapping(atomorder)
   allocate(parts(size(self%eltypepartition%parts)))

   do i = 1, size(self%eltypepartition%parts)
      parts(i)%atomidcs = atomunorder(self%eltypepartition%parts(i)%atomidcs(:))
   end do

end function

function get_eltypepartlens(self) result(eltypepartlens)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer, allocatable :: eltypepartlens(:)

   eltypepartlens = self%eltypepartition%get_lenlist()

end function

function get_fragparts(self, atomorder) result(parts)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   integer, allocatable :: atomunorder(:)
   type(Part), allocatable :: parts(:)

   atomunorder = inverse_mapping(atomorder)
   allocate(parts(size(self%molfragpartition%parts)))

   do i = 1, size(self%molfragpartition%parts)
      parts(i)%atomidcs = atomunorder(self%molfragpartition%parts(i)%atomidcs(:))
   end do

end function

function get_fragroots(self, atomorder) result(fragroots)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   integer, allocatable :: order(:)
   integer, allocatable :: fragroots(:)
   integer, allocatable :: atommnatypes(:)
   integer, allocatable :: mnatypepartlens(:)
   type(Part), allocatable :: fragparts(:)

   fragparts = self%get_fragparts(atomorder)
   atommnatypes = self%get_atommnatypes(atomorder)
   mnatypepartlens = self%get_mnatypepartlens()
   allocate(fragroots(size(fragparts)))

   do i = 1, size(fragparts)
      order = sorted_order(mnatypepartlens(atommnatypes(fragparts(i)%atomidcs)))
      fragroots(i) = fragparts(i)%atomidcs(order(1))
   end do

end function

function get_default_atomeltypes(self) result(atomeltypes)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer, allocatable :: atomeltypes(:)

   atomeltypes = self%atoms(:)%eltype

end function

function get_sorted_atomeltypes(self, atomorder) result(atomeltypes)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer, allocatable :: atomeltypes(:)

   atomeltypes = self%atoms(atomorder)%eltype

end function

function get_default_atommnatypes(self) result(atommnatypes)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer, allocatable :: atommnatypes(:)

   atommnatypes = self%atoms(:)%mnatype

end function

function get_sorted_atommnatypes(self, atomorder) result(atommnatypes)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer, allocatable :: atommnatypes(:)

   atommnatypes = self%atoms(atomorder)%mnatype

end function

subroutine set_weights(self, weights)
   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: weights(self%natom)

   self%atoms(:)%weight = weights(:)

end subroutine

function get_default_atomweights(self) result(weights)
   class(Molecule), intent(in) :: self
   ! Local variables
   real(wp), allocatable :: weights(:)

   weights = self%atoms(:)%weight

end function

function get_sorted_atomweights(self, atomorder) result(weights)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   real(wp), allocatable :: weights(:)

   weights = self%atoms(atomorder)%weight

end function

subroutine set_coords(self, coords)
   class(Molecule), intent(inout) :: self
   real(wp), intent(in) :: coords(3, self%natom)
   ! Local variables
   integer :: i

   do i = 1, self%natom
      self%atoms(i)%coords = coords(:, i)
   end do

end subroutine

function get_default_coords(self) result(coords)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i
   real(wp), allocatable :: coords(:, :)

   allocate(coords(3, self%natom))

   do i = 1, self%natom
      coords(:, i) = self%atoms(i)%coords(:)
   end do

end function

function get_sorted_coords(self, atomorder) result(coords)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   real(wp), allocatable :: coords(:, :)

   allocate(coords(3, self%natom))

   do i = 1, self%natom
      coords(:, i) = self%atoms(atomorder(i))%coords(:)
   end do

end function

function get_default_adjmat(self) result(adjmat)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i, k
   type(Atom) :: iatom
   logical, allocatable :: adjmat(:, :)

   allocate(adjmat(self%natom, self%natom))

   adjmat(:, :) = .false.

   do i = 1, self%natom
      iatom = self%atoms(i)
      do k = 1, size(iatom%adjlist)
         adjmat(i, iatom%adjlist(k)) = .true.
      end do
   end do

end function

function get_sorted_adjmat(self, atomorder) result(adjmat)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i, k
   type(Atom) :: atom_i
   integer, allocatable :: atomunorder(:)
   logical, allocatable :: adjmat(:, :)

   atomunorder = inverse_mapping(atomorder)
   allocate(adjmat(self%natom, self%natom))
   adjmat(:, :) = .false.

   do i = 1, self%natom
      atom_i = self%atoms(atomorder(i))
      do k = 1, size(atom_i%adjlist)
         adjmat(i, atomunorder(atom_i%adjlist(k))) = .true.
      end do
   end do

end function

function get_default_nadjs(self) result(nadjs)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: nadjs(:)

   allocate(nadjs(size(self%atoms)))
   do i = 1, size(self%atoms)
      nadjs(i) = size(self%atoms(i)%adjlist)
   end do

end function

function get_sorted_nadjs(self, atomorder) result(nadjs)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   integer, allocatable :: nadjs(:)

   allocate(nadjs(size(self%atoms)))

   do i = 1, size(self%atoms)
      nadjs(i) = size(self%atoms(atomorder(i))%adjlist)
   end do

end function

subroutine set_adjlists(self, nadjs, adjlists)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: nadjs(:)
   integer, intent(in) :: adjlists(:, :)
   ! Local variables
   integer :: i

   do i = 1, self%natom
      self%atoms(i)%adjlist = adjlists(:nadjs(i), i)
   end do

end subroutine

function get_default_adjlists(self) result(adjlists)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: adjlists(:, :)

   allocate(adjlists(maxcoord, self%natom))

   do i = 1, self%natom
      adjlists(:size(self%atoms(i)%adjlist), i) = self%atoms(i)%adjlist
   end do

end function

function get_sorted_adjlists(self, atomorder) result(adjlists)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   type(Atom) :: atom_i
   integer, allocatable :: atomunorder(:)
   integer, allocatable :: adjlists(:, :)

   atomunorder = inverse_mapping(atomorder)
   allocate(adjlists(maxcoord, self%natom))

   do i = 1, self%natom
      atom_i = self%atoms(atomorder(i))
      adjlists(:size(atom_i%adjlist), i) = atomunorder(atom_i%adjlist)
   end do

end function

subroutine set_atomneimnatypepartlens(self, natomneimnatypes, atomneimnatypepartlens)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: natomneimnatypes(:)
   integer, intent(in) :: atomneimnatypepartlens(:, :)
   ! Local variables
   integer :: i

   do i = 1, size(self%atoms)
      if (.not. allocated(self%atoms(i)%atomneimnatypepartlens)) then
         allocate(self%atoms(i)%atomneimnatypepartlens(natomneimnatypes(i)))
      end if
      self%atoms(i)%atomneimnatypepartlens = atomneimnatypepartlens(:natomneimnatypes(i), i)
   end do

end subroutine

function get_default_atomneimnatypepartlens(self) result(atomneimnatypepartlens)
   class(Molecule), intent(inout) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: atomneimnatypepartlens(:, :)

   allocate(atomneimnatypepartlens(maxcoord, size(self%atoms)))

   do i = 1, size(self%atoms)
      atomneimnatypepartlens(:size(self%atoms(i)%atomneimnatypepartlens), i) = self%atoms(i)%atomneimnatypepartlens
   end do

end function

function get_sorted_atomneimnatypepartlens(self, atomorder) result(atomneimnatypepartlens)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   type(Atom) :: atom_i
   integer, allocatable :: atomneimnatypepartlens(:, :)

   allocate(atomneimnatypepartlens(maxcoord, size(self%atoms)))

   do i = 1, size(self%atoms)
      atom_i = self%atoms(atomorder(i))
      atomneimnatypepartlens(:size(atom_i%atomneimnatypepartlens), i) = atom_i%atomneimnatypepartlens
   end do

end function

function get_default_natomneimnatypes(self) result(natomneimnatypes)
   class(Molecule), intent(inout) :: self
   ! Local variables
   integer :: i
   integer, allocatable :: natomneimnatypes(:)

   allocate(natomneimnatypes(size(self%atoms)))

   do i = 1, size(self%atoms)
      natomneimnatypes(i) = size(self%atoms(i)%atomneimnatypepartlens)
   end do

end function

function get_sorted_natomneimnatypes(self, atomorder) result(natomneimnatypes)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: atomorder(:)
   ! Local variables
   integer :: i
   type(Atom) :: atom_i
   integer, allocatable :: natomneimnatypes(:)

   allocate(natomneimnatypes(size(self%atoms)))

   do i = 1, size(self%atoms)
      atom_i = self%atoms(atomorder(i))
      natomneimnatypes(i) = size(atom_i%atomneimnatypepartlens)
   end do

end function

subroutine print_atom(self, ind, outLvl)
   class(Atom), intent(in) :: self
   integer, intent(in) :: ind
   integer, intent(in), optional :: outLvl
   ! Local variables
   character(255) :: frmt
   character(2) :: num
   integer :: outLevel

! *** code to manage unitlbl pending ***

   outLevel = 1
   if (present(outLvl)) outLevel = outLvl
   
   select case (outLevel)

   case (1)
      if (size(self%adjlist) == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (stderr, frmt) ind, ": ", self%elnum, self%label, self%eltype, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') size(self%adjlist)
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%elnum, self%label, self%eltype, &
               self%weight, " {", self%coords(:), " } [", &
               self%adjlist(:size(self%adjlist)), " ]"
      end if

!  case (2)
!
   case default
      if (size(self%adjlist) == 0) then
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a)'
         write (stderr, frmt) ind, ": ", self%elnum, self%label, self%eltype, &
                     self%weight, " {", self%coords(:), " }"
      else
         write (num, '(i0)') size(self%adjlist)
         frmt = '(i3,2a,2i3,f7.2,a,3f8.3,a,'//num//'i3,a)'
         write (stderr, frmt) ind, ": ", self%elnum, self%label, self%eltype, &
               self%weight, " {", self%coords(:), " } [", &
               self%adjlist(:size(self%adjlist)), " ]"
      end if
   end select

end subroutine

subroutine print_molecule(self)
   class(Molecule), intent(in) :: self
   ! Local variables
   integer :: i

! *** code to manage unitlbl pending ***
   
   write (stderr, '(a,i0,a)') "Contents of molecule structure:   (", &
                                         self%natom, " atoms)"
   write (stderr, '(2a)') 'Title: ', self%title
   write (stderr, '(a,a4,a5,a12,a7,2a14)') "ind:", "lbl", "elnum", "eltype", &
                                          "weight", "{ coords }", "[ adjlist ]"

   do i = 1, self%natom
      call self%atoms(i)%print(i)
   end do

end subroutine

function bonded(self, idx1, idx2) result(isbond)
   class(Molecule), intent(in) :: self
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: i
   logical :: isbond, found1, found2
   integer, allocatable :: adjlist1(:), adjlist2(:)

   allocate(adjlist1(maxcoord), adjlist2(maxcoord))

! copy arrays of adjlist
   adjlist1 = self%atoms(idx1)%adjlist
   adjlist2 = self%atoms(idx2)%adjlist

! initialization
   found1 = .false.
   found2 = .false.

! check cross reference of idx1 and idx2 in both adjlists
   do i = 1, size(self%atoms(idx1)%adjlist)
      if (idx2 == adjlist1(i)) then
         found1 = .true.
         exit
      end if
   end do
   do i = 1, size(self%atoms(idx2)%adjlist)
      if (idx1 == adjlist2(i)) then
         found2 = .true.
         exit
      end if
   end do

! report or stop program
   if (found1 .and. found2) then
      isbond = .true.
   else if (found1 .or. found2) then
      write (stderr, '(a,i0,2x,i0)') 'Inconsistent bond for atoms: ', idx1, idx2
      stop
   else
      isbond = .false.
   end if

end function

subroutine remove_bond(self, idx1, idx2)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: i, pos1, pos2, nadj1, nadj2
   integer, allocatable :: adjlist1(:), adjlist2(:)

   allocate(adjlist1(maxcoord), adjlist2(maxcoord))

! copy adjlist arrays
   nadj1 = size(self%atoms(idx1)%adjlist)
   nadj2 = size(self%atoms(idx2)%adjlist)
   adjlist1 = self%atoms(idx1)%adjlist
   adjlist2 = self%atoms(idx2)%adjlist

! initialization
   pos1 = 0   ! position of idx2 in adjlist of atom 1
   pos2 = 0   ! position of idx1 in adjlist of atom 2

! find position of idx2 and idx1 in adjlist of atoms idx1 and idx2, resp.
   do i = 1, nadj1
      if (idx2 == adjlist1(i)) then
         pos1 = i
         exit
      end if
   end do
   do i = 1, nadj2
      if (idx1 == adjlist2(i)) then
         pos2 = i
         exit
      end if
   end do

! delete idx2 and idx1 from the ajdlists where they appear
   if ((pos1 /= 0) .and. (pos2 /= 0)) then
      nadj1 = nadj1 - 1
      do i = pos1, nadj1
         adjlist1(i) = adjlist1(i+1)
      end do
      nadj2 = nadj2 - 1
      do i = pos2, nadj2
         adjlist2(i) = adjlist2(i+1)
      end do
! update neighbor arrays for atoms idx1 and idx2
      self%atoms(idx1)%adjlist = adjlist1(:nadj1)
      self%atoms(idx2)%adjlist = adjlist2(:nadj2)
   else
      write (stderr, '(a,i0,2x,i0)') 'Error: atoms not bonded: ', idx1, idx2
   end if

end subroutine

subroutine add_bond(self, idx1, idx2)
   class(Molecule), intent(inout) :: self
   integer, intent(in) :: idx1, idx2
   ! Local variables
   integer :: pos1, pos2, nadj1, nadj2
   integer, allocatable :: adjlist1(:), adjlist2(:)

   allocate(adjlist1(maxcoord), adjlist2(maxcoord))

! copy array of adjlist
   nadj1 = size(self%atoms(idx1)%adjlist)
   nadj2 = size(self%atoms(idx2)%adjlist)
   adjlist1 = self%atoms(idx1)%adjlist
   adjlist2 = self%atoms(idx2)%adjlist

! initialization
   pos1 = nadj1
   pos2 = nadj2

   if (.not. self%bonded(idx1, idx2)) then
! indices in adjlist are supposed to be sorted; inserting new indices
!      size(self%atoms(idx1)%adjlist) = size(self%atoms(idx1)%adjlist) + 1
      nadj1 = nadj1 + 1
! find position to insert idx2 and shift indices greater than idx2
      do while ((pos1 >= 1) .and. (idx2 < adjlist1(pos1)))
         adjlist1(pos1+1) = adjlist1(pos1)
         pos1 = pos1 - 1
      end do
      adjlist1(pos1+1) = idx2
      
      nadj2 = nadj2 + 1
! find position to insert idx1 and shift indices greater than idx1
      do while ((pos2 >= 1) .and. (idx1 < adjlist2(pos2)))
         adjlist2(pos2+1) = adjlist2(pos2)
         pos2 = pos2 - 1
      end do
      adjlist2(pos2+1) = idx1
! update neighbor arrays for atoms in idx1 and idx2
      self%atoms(idx1)%adjlist = adjlist1(:nadj1)
      self%atoms(idx2)%adjlist = adjlist2(:nadj2)
   else
      write (stderr, '(a,i0,2x,i0)') "Error: atoms already bonded: ", idx1, idx2
   end if

end subroutine

end module
