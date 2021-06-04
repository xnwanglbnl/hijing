The compile script with gortran is following:

gfortran -O -fno-automatic  main.f  hijing2.f  hipyset1.4.f \
 -lnsl -lg2c -L ./  -lkernlib -lpdflib804 -lmathlib -lpacklib \
 -o main.out

There are four lib files needed in this compiling: libkernlib.a, libpdflib804.a, libmathlib.a, libpacklib.a. One can find them on CERNLIB website: https://cernlib.web.cern.ch/cernlib/

The default directory of these lib files in this script are current working directory. 
If the lib files are put on  any other directory, one can just amend " -L  ./" as " -L   any-other/" in this script to tell compiler the correct directory.

As an example, one main program dAu.f is presented here. One can simply run the script file runp_gfortran to compile it using the following command:

./runp_gfortran dAu

After compiling, one can get a excutable file named dAu.out. Then the following command can give the final result file:
./dAu.out

 
