CC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
          CFLAGS = -I ./eigen-eigen/
else
          CFLAGS = -I ./eigen-eigen/
endif

RM = /bin/rm -f 
all: main 
main: bezier.o
    $(CC) -O3 $(CFLAGS) -o bezier bezier.o $(LDFLAGS)
raytrace.o: bezier.cpp
    $(CC) -O3 $(CFLAGS) -c bezier.cpp -o bezier.o
clean: 
    $(RM) *.o as3



