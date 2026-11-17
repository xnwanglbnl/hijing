HIJING/BBbar v2.0 (2008 code, code version 2.1) -- build and run on Linux
=========================================================================
Model: V. Topor Pop et al., Phys. Rev. C 70, 064906 (2004),
       https://arxiv.org/abs/nucl-th/0407095

1. Install the build tools (Debian/Ubuntu):
     sudo apt-get install build-essential gfortran perl git xutils-dev
   The program itself needs only gfortran and make.  perl, git and
   xutils-dev (imake, makedepend) are needed for the CERNLIB build in
   step 2 only.

2. Build CERNLIB (packlib + kernlib):
   The required cernlib libraries can be installed separately, that can be obtained from
   https://github.com/apc-llc/cernlib  if the following steps do not work. 

     cd cernlib_build
     sh build_cernlib_packlib.sh
     cd ..
   (If CERNLIB is already installed (this is the best option), skip this and use
    make CERNLIB_DIR=/path/to/cernlib/lib  in step 3.)


3. Build:
     make

4. Run:
     ./hbb2test < au200test.dat > au200test.log
   Output: auau200test.paw (HBOOK histograms), auau200test.num, and the log.
