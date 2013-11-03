CC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -I./include/ -I/usr/X11/include -DOSX \
			-I ./eigen-eigen/
	LDFLAGS = -framework GLUT -framework OpenGL \
			-L "/System/Library/Frameworks/OpenGL.framework/Libraries" \
			-lGL -lGLU -lm -lstdc++
else
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -Iglut-3.7.6-bin \
			-I ./eigen-eigen/
	LDFLAGS = -lglut -lGLU
endif

RM = /bin/rm -f 
all: main 
main: bezier.o
	$(CC) -O3 $(CFLAGS) -o bezier bezier.o $(LDFLAGS) 
bezier.o: bezier.cpp 
	$(CC) -O3 $(CFLAGS) -c bezier.cpp -o bezier.o
clean: 
	$(RM) *.o as3
