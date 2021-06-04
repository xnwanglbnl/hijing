hijing.o:hijing.f
	f77 -c hijing.f -o hijing.o

hipyset.o:hipyset.f
	f77 -c hipyset.f -o hipyset.o

test.o:test.f
	f77 -c test.f -o test.o

test: test.o hijing.o hipyset.o
	f77 -o test test.o hijing.o hipyset.o -L/usr/cern/lib -lpdflib -lkernlib

