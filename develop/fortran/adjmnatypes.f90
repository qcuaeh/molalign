function get_adjmnatypes(self, mnatypes) result(adjmnatypes)
   class(molecule_type), intent(in) :: self
   type(partition_type), intent(in) :: mnatypes
   ! Result variable
   type(partition_type), allocatable :: adjmnatypes(:)
   ! Local variables
   integer :: h, i
   type(atom_type) :: atom
   integer, allocatable :: adjlistsubset(:)

   do i = 1, size(self%atoms)
      atom = self%atoms(i)
      allocate (adjmnatypes(i)%parts(size(mnatypes%parts)))
      do h = 1, size(mnatypes%parts)
         adjlistsubset = intersection(atom%adjlist, mnatypes%parts(h)%indices, size(self%atoms))
         allocate (adjmnatypes(i)%parts(h)%indices(size(adjlistsubset)))
         adjmnatypes(i)%parts(h)%indices = adjlistsubset
      end do
   end do

end function
