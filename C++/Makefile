INCFLAGS = -IQuaternions -I/opt/local/include -I/usr/local/include
LIBFLAGS = -L/opt/local/lib -L/usr/local/lib
ifdef GSL_HOME
	INCFLAGS = -I${GSL_HOME}/include ${INCFLAGS}
	LIBFLAGS = -L${GSL_HOME}/lib ${LIBFLAGS}
endif
# Set compiler name and optimization flags here, if desired
C++ = g++
OPT = -O3 -fopenmp -Wall -Wno-deprecated -lgsl
## DON'T USE -ffast-math in OPT


#############################################################################
## The following are pretty standard and probably won't need to be changed ##
#############################################################################

# Tell 'make' not to look for files with the following names
.PHONY : all

all : PNEvolution.o

# This is how to build those object files
%.o : %.cpp %.hpp
	$(C++) $(OPT) -c $(INCFLAGS) $< -o $@
