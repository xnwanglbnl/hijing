# HIJING/BBbar v2.0 (2008 code) -- Linux / gfortran build
# CERNLIB_DIR: directory with libpacklib.a and libkernlib.a
#   make                                (default: built by cernlib_build/build_cernlib_packlib.sh)
#   make CERNLIB_DIR=/cern/2006/lib     (existing CERNLIB)
# Original 2008 Makefile (f77 + CERNLIB 2005): Makefile_original_2008

CERNLIB_DIR ?= $(CURDIR)/cernlib_build/cernlib/2006/lib

PROG = hbb2test
OBJ  = xhijing.bb2.1.o hijing1.382.bb2.1.o hipyset.bb2.1.o junction2.1.o gfortran_getarg_shim.o

FC      = gfortran
FFLAGS  = -O -std=legacy -fno-automatic -finit-local-zero -ffixed-line-length-132
CC      = gcc
CFLAGS  = -O
LIBS    = -L$(CERNLIB_DIR) -lpacklib -lkernlib -Wl,--allow-multiple-definition

all: $(PROG)

$(PROG): $(OBJ)
	$(FC) -o $@ $(OBJ) $(LIBS)

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(PROG)
