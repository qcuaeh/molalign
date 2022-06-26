module options

integer, parameter :: wp = 8
integer, parameter :: lbllen = 16
integer, parameter :: optlen = 128
integer, parameter :: ttllen = 256
integer, parameter :: maxcoord = 16

integer maxtrial, maxmatch
logical live, remap, iterative, biased, trialing, matching, testing
character(optlen) weighter, outformat
real(wp) biasscale, tolerance

end module
