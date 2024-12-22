module eigen
use parameters

implicit none

external ssyev
external dsyev

private
public leasteigval
public leasteigvec

interface leasteigval
   module procedure leasteigval_sp
   module procedure leasteigval_dp
end interface

interface leasteigvec
   module procedure leasteigvec_sp
   module procedure leasteigvec_dp
end interface

contains

function leasteigval_sp(a) result(leasteigval)
   real(real32), intent(inout) :: a(4, 4)
   integer :: info
   real(real32) :: leasteigval
   real(real32) :: w(4), work(25)
   call ssyev('N', 'U', 4, a, 4, w, work, 25, info)
   leasteigval = w(1)
end function

function leasteigval_dp(a) result(leasteigval)
   real(real64), intent(inout) :: a(4, 4)
   integer :: info
   real(real64) :: leasteigval
   real(real64) :: w(4), work(25)
   call dsyev('N', 'U', 4, a, 4, w, work, 25, info)
   leasteigval = w(1)
end function

function leasteigvec_sp(a) result(leasteigvec)
   real(real32), intent(inout) :: a(4, 4)
   integer :: info
   real(real32) :: w(4), work(25)
   real(real32) :: leasteigvec(4)
   call ssyev('V', 'U', 4, a, 4, w, work, 25, info)
   leasteigvec = a(:, 1)
end function

function leasteigvec_dp(a) result(leasteigvec)
   real(real64), intent(inout) :: a(4, 4)
   integer :: info
   real(real64) :: leasteigvec(4)
   real(real64) :: w(4), work(25)
   call dsyev('V', 'U', 4, a, 4, w, work, 25, info)
   leasteigvec = a(:, 1)
end function 

end module
