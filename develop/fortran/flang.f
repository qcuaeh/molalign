integer a, b
(integer, allocatable, rank=1) array1, array2
(integer, kind=int64, allocatable, rank=2) int64_array2d
(integer, kind=int64, rank=1, dimensions=[2]) int64_array_3x4
(type=Permutation, pointer, rank=1) permutation_array_pointer
((integer, pointer), rank=1) integer_pointer_array
((type=Permutation, pointer), rank=1) permutation_pointer_array
((integer, allocatable), allocatable) nested_integer_list ! List of lists of integers
(((integer, allocatable), allocatable), allocatable) nested_nested_integer_list ! List of lists of lists of integers

allocate( int64_array2d, [3, 4])
