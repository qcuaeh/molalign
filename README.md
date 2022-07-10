Compile the program
===================

Run `./build` without arguments to build the program with the highest
optimization level. The script accepts the following options:

&ensp; -debug : Build the debug version  
&ensp; -rebuild : Recompile all source files from scratch  

The executables are created in the "bin" directory and the object
files are stored in the **\_build_dir** directory.

Compile the library
===================

Run `./build-lib` without arguments to compile the library with the highest
optimization level. The script accepts the same options as **build** and the
same notices apply.

Run the program
===============

#### Program options

&ensp; -maps INTEGER : Number of recorded maps  
&ensp; -out xyz|mol2 : Output format (XYZ or MOL2)  
&ensp; -live : Print live stats  
&ensp; -iter : Iterative trial step  
&ensp; -bias REAL : Tolerance for biasing  
&ensp; -weight STRING : Weighting property  
&ensp; -count INTEGER : Maximum map counting  
&ensp; -trials INTEGER : Maximum number of trials  
&ensp; -scale REAL : Length scale  
&ensp; -test : Repeateble run for testing  
 
#### Examples
 
./bin/ralign tests/r005/100cobalt_j5.xyz -iter -weight mass -bias 0.17 -test -trials 1000 -count 10
