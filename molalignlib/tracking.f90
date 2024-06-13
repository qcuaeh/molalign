! purpose: determines the most convenient atom to start runs over the structure
!          of each fragment in a molecule. (returns nfrag, fragroots)
module tracking

use stdio
use types
use bounds
use sorting
implicit none
private
public findmolfrags

contains

!CZGC: new definition
subroutine findmolfrags (mol)
!subroutine findmolfrags (mol, nfrag, fragroots)
!CZGC: old definition
!subroutine findmolfrags (natom, nadjs, adjlists, ntype, typelenlist, &
!                  nequiv, equivlenlist, nfrag, fragroots)

   type(Molecule), intent(inout) :: mol
   integer :: natom, ntype, nequiv
   integer, dimension(mol%natom) :: nadjs, typelenlist, equivlenlist
   integer, dimension(maxcoord, mol%natom) :: adjlists
!   integer, intent(in) :: natom, ntype, nequiv
!   integer, dimension(:), intent(in) :: nadjs, typelenlist, equivlenlist
!   integer, dimension(:, :), intent(in) :: adjlists
!   integer, intent(out) :: nfrag, fragroots(mol%natom)
   integer :: nfrag, fragroots(mol%natom)

   integer :: i, h, n, f, firstn, prevn, offset
   integer, dimension(mol%natom) :: order
   integer, dimension(mol%natom) :: blkidx
   integer, dimension(mol%natom) :: eqvidx
   integer, dimension(mol%natom) :: fragid
   integer, dimension(mol%natom) :: fragsize
   integer, dimension(mol%natom,mol%natom) :: fragind, fragblock, fragequiv, fragcoord
   logical :: track(mol%natom), printInfo = .false.
   character(80) :: fmtstr
   character(3) :: numat

!CZGC: temporal variable
   natom = mol%natom
   ntype = mol%get_ntype()
   typelenlist = mol%get_typelenlist()
   nequiv = mol%get_nequiv()
   equivlenlist = mol%get_equivlenlist()
   nadjs = mol%get_nadjs()
   adjlists = mol%get_adjlists()

   ! set atoms block indices

   offset = 0
   do h = 1, ntype
      blkidx(offset+1:offset+typelenlist(h)) = h
      offset = offset + typelenlist(h)
   end do

   ! set atoms equivalence indices

   offset = 0
   do h = 1, nequiv
      eqvidx(offset+1:offset+equivlenlist(h)) = h
      offset = offset + equivlenlist(h)
   end do

   ! initialization

   nfrag = 0
   fragroots(:) = 1

   order = [ (n, n = 1, natom) ]
   track(:) = .false.
   fragid(:) = 0
   fragsize(:) = 0
   fragind(:,:) = 0
   fragblock(:,:) = 0
   fragequiv(:,:) = 0
   fragcoord(:,:) = 0

   ! detect fragments and populate frag arrays
   ! 'order' could be shuffled at this point to randomize the fragments order
   n = 1
   do while ( n <= natom )
      if ( track(order(n)) ) then
         n = n + 1
      else
         nfrag = nfrag + 1
         call recrun (order(n), track)
         n = 1
      end if
   end do

   ! print initial information
   if ( printInfo ) then
      write (numat, '(i0)') natom
      fmtstr = '(i2,x,a,'//numat//'i3)'
      print *, "fragments detected: ", nfrag
      print fmtstr, 0, "sequence ini:   ", (n, n = 1, natom)
      print fmtstr, 0, "fragids ini:    ", fragid(:natom)
      do f = 1, nfrag
         print fmtstr, f, "frag inds ini:  ", fragind(:fragsize(f),f)
         print fmtstr, f, "frag block size:", fragblock(:fragsize(f),f)
         print fmtstr, f, "frag equiv size:", fragequiv(:fragsize(f),f)
         print fmtstr, f, "frag nadjs:    ", fragcoord(:fragsize(f),f)
      end do
   end if

   ! sorts sequentially by typelenlist, equivlenlist, and nadjs for each fragment
   order = [ (n, n = 1, natom) ]
   do f = 1, nfrag
      firstn = fragsize(f)
      call groupminsort (fragblock(:firstn,f), order(:firstn), firstn)
      fragind(:firstn,f) = fragind(order(:firstn),f)
      fragblock(:firstn,f) = fragblock(order(:firstn),f)
      fragequiv(:firstn,f) = fragequiv(order(:firstn),f)
      fragcoord(:firstn,f) = fragcoord(order(:firstn),f)
      call groupminsort (fragequiv(:firstn,f), order(:firstn), firstn)
      fragind(:firstn,f) = fragind(order(:firstn),f)
      fragblock(:firstn,f) = fragblock(order(:firstn),f)
      fragequiv(:firstn,f) = fragequiv(order(:firstn),f)
      fragcoord(:firstn,f) = fragcoord(order(:firstn),f)
      prevn = firstn
      call groupminsort (fragcoord(:firstn,f), order(:firstn), firstn)
      if ( firstn /= prevn ) order(:prevn) = invsort (order(:prevn), prevn)
      fragind(:prevn,f) = fragind(order(:prevn),f)
      fragblock(:firstn,f) = fragblock(order(:firstn),f)
      fragequiv(:firstn,f) = fragequiv(order(:firstn),f)
      fragcoord(:prevn,f) = fragcoord(order(:prevn),f)
      fragroots(f) = fragind(1,f)
   end do

   ! print final vectors
   if ( printInfo ) then
      write (numat, '(i0)') natom
      fmtstr = '(i2,x,a,'//numat//'i3)'
      do f = 1, nfrag
         print fmtstr, f, "frag inds fin:  ", fragind(:fragsize(f),f)
         print fmtstr, f, "frag block size:", fragblock(:fragsize(f),f)
         print fmtstr, f, "frag equiv size:", fragequiv(:fragsize(f),f)
         print fmtstr, f, "frag nadjs:    ", fragcoord(:fragsize(f),f)
         print '(i2,1x,a,1x,i10)', f, "fragment start ind:", fragroots(f)
      end do
   end if

   ! register nfrag and fragroots in the mol strucutre
   if ( mol%nfrag == 0 ) then
      allocate(mol%fragroots(mol%natom))
   end if
   mol%nfrag = nfrag
   mol%fragroots = fragroots

contains

   ! runs recursivelly over the structure and populates arrays
   recursive subroutine recrun (n, track)
      integer, intent(in) :: n
      logical, dimension(:), intent(inout) :: track
      integer i

      if ( track(n) ) return
      track(n) = .true.

      fragsize(nfrag) = fragsize(nfrag) + 1
      fragid(n) = nfrag
      fragind(fragsize(nfrag),nfrag) = n
      fragblock(fragsize(nfrag),nfrag) = typelenlist(blkidx(n))
      fragequiv(fragsize(nfrag),nfrag) = equivlenlist(eqvidx(n))
      fragcoord(fragsize(nfrag),nfrag) = nadjs(n)
      
      do i = 1, nadjs(n)
         call recrun (adjlists(i, n), track)
      end do
   end subroutine recrun

   ! sorts vec and returns the order and number of elements in the first group
   subroutine groupminsort (vec, order, firstn)
!        integer, intent(in) :: n
      integer, dimension(:), intent(in) :: vec
      integer, dimension(:), intent(out) :: order
      integer, intent(inout) :: firstn
      integer i, n, firstitem

      order = sorted_order (vec, firstn)
      n = firstn
      firstitem = vec(order(1))
      do i = 2, n
         if ( vec(order(i)) /= firstitem ) then
            firstn = i - 1
            return
         end if
      end do
   end subroutine groupminsort

   ! inverts the order of a vector
   function invsort (vec, n) result(invec)
      integer, intent(in) :: vec(:), n
      integer :: i, p, invec(n)
      p = n
      do i = 1, n
         invec(p) = vec(i)
         p = p - 1
      end do
   end function invsort

end subroutine findmolfrags

end module
