OBJ = hijing1.34.bb.o hipyset.bb.o xhijing.bb.o junction.o
PROG = hijbb

#FC = f77 -O
FC = gfortran -O

LIBS = -L/star/data02/pwg/ashishpandav/ROOTAREA/HIJING_B/CERN_LIB/cernlib/2006/lib -lpacklib -lpawlib -lkernlib 
#########################
all: $(PROG)

$(PROG): $(OBJ) 
	$(FC) -o $@ $(OBJ) $(LIBS) 
