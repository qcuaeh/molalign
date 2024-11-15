Integer a, b
(Integer, rank=1) array1, array2
(Integer, rank=2, allocatable) int64_array2d
((Integer, rank=1, allocatable), rank=1, allocatable) nested_integer_list ! List of lists of integers
(Real:real64, rank=1, pointer) real_array_pointer ! Pointer to array of reals
(Real:real64, pointer, rank=1) real_pointer_array ! Array of pointers to reals
(Complex, pointer, rank=1) complex_pointer_array

allocate( 10) array1, array2 
allocate( 10, 10) int64_array2d
print( int64_array2d[1..3, 1..5])
newtype => partition.new_part(self.parts[i].num_elems)
a = 8
b = 8
a /= 2
if a \= b
   print( 'a \= b')
end if
if a \= b then print( 'a \= b')
for i = 1..array1:size
   i+=1
end for
for i = 1..array:size do i += 1
for e in array1
   print( e)
end for
subroutine power(in:n, in:x, out:y)
