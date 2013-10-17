INCFLAGS = -IQuaternions -I/opt/local/include -I/usr/local/include
LIBFLAGS = -L/opt/local/lib -L/usr/local/lib -lgsl
ifdef GSL_HOME
	INCFLAGS = -I${GSL_HOME}/include ${INCFLAGS}
	LIBFLAGS = -L${GSL_HOME}/lib ${LIBFLAGS}
endif
# Set compiler name and optimization flags here, if desired
C++ = g++
OPT = -O3 -fopenmp -Wall -Wno-deprecated
## DON'T USE -ffast-math in OPT


#############################################################################
## The following are pretty standard and probably won't need to be changed ##
#############################################################################

# Tell 'make' not to look for files with the following names
.PHONY : all

all : PNEvolution.o EvolvePN PNEvolution_Q.o EvolvePN_Q

# This is how to build those object files
PNEvolution.o : PNEvolution.cpp PNEvolution.hpp PNApproximants.ipp
	$(C++) $(OPT) -c $(INCFLAGS) $< -o $@

# This is how to build those object files
PNEvolution_Q.o : PNEvolution_Q.cpp PNEvolution.hpp PNApproximants_Q.ipp
	$(C++) $(OPT) -c $(INCFLAGS) $< -o $@

# This is how to build those object files
%.o : %.cpp %.hpp
	$(C++) $(OPT) -c $(INCFLAGS) $< -o $@

EvolvePN : EvolvePN.cpp PNEvolution.cpp PNEvolution.hpp PNApproximants.ipp
	$(C++) $(OPT) $(INCFLAGS) $(LIBFLAGS) PNEvolution.o Quaternions/Quaternions.o Quaternions/IntegrateAngularVelocity.o EvolvePN.cpp -o EvolvePN

EvolvePN_Q : EvolvePN_Q.cpp PNEvolution_Q.cpp PNEvolution.hpp PNApproximants_Q.ipp
	$(C++) $(OPT) $(INCFLAGS) $(LIBFLAGS) PNEvolution_Q.o Quaternions/Quaternions.o Quaternions/IntegrateAngularVelocity.o EvolvePN_Q.cpp -o EvolvePN_Q
