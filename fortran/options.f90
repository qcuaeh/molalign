module options

integer, parameter :: wp = 8
integer maxtrial, maxmatch
logical live, iterative, biased, trialing, matching, testing
character(32) weighting, outformat
real(wp) scaling, tolerance

end module
