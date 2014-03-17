INCFLAGS = -IQuaternions -I/usr/local/include
LIBFLAGS = -L/usr/local/lib -lgsl
test -d /opt/local/include && INCFLAGS = -IQuaternions -I/opt/local/include -I/usr/local/include
test -d /opt/local/lib && LIBFLAGS = -L/opt/local/lib -L/usr/local/lib -lgsl
ifdef GSL_HOME
	INCFLAGS = -I${GSL_HOME}/include ${INCFLAGS}
	LIBFLAGS = -L${GSL_HOME}/lib ${LIBFLAGS}
endif
# Set compiler name and optimization flags here, if desired
C++ = g++
OPT = -O3 -Wall -Wno-deprecated
## DON'T USE -ffast-math in OPT


#############################################################################
## The following are pretty standard and probably won't need to be changed ##
#############################################################################

# Tell 'make' not to look for files with the following names
.PHONY : python all clean cpp

python :
	python setup.py install --user

all : cpp python

cpp : PNEvolution.o EvolvePN PNEvolution_Q.o EvolvePN_Q PNWaveformModes.o

# This is how to build those object files
PNEvolution.o : PNEvolution.cpp PNEvolution.hpp PNApproximants.ipp
	$(C++) $(OPT) -c $(INCFLAGS) $< -o $@

# This is how to build those object files
PNEvolution_Q.o : PNEvolution_Q.cpp PNEvolution.hpp PNApproximants_Q.ipp
	$(C++) $(OPT) -c $(INCFLAGS) $< -o $@

# This is how to build those object files
%.o : %.cpp %.hpp
	$(C++) $(OPT) -c $(INCFLAGS) $< -o $@

EvolvePN : EvolvePN.cpp PNEvolution.o PNEvolution_Q.o PNEvolution.hpp PNApproximants.ipp
	$(C++) $(OPT) $(INCFLAGS) $(LIBFLAGS) PNEvolution.o PNEvolution_Q.o Quaternions/Quaternions.o Quaternions/IntegrateAngularVelocity.o EvolvePN.cpp -o EvolvePN

EvolvePN_Q : EvolvePN_Q.cpp PNEvolution.o PNEvolution_Q.o PNEvolution.hpp PNApproximants_Q.ipp
	$(C++) $(OPT) $(INCFLAGS) $(LIBFLAGS) PNEvolution.o PNEvolution_Q.o Quaternions/Quaternions.o Quaternions/IntegrateAngularVelocity.o EvolvePN_Q.cpp -o EvolvePN_Q

clean:
	/bin/rm -f *.o EvolvePN EvolvePN_Q
	/bin/rm -rf build _PNEvolution.so PNEvolution.py PNEvolution/__init__.py

# Mike likes to use python's virtual environments to build with
# various versions and arrangements of python.  As such, he finds that
# the following is, unfortunately, necessary for him.  This makes him
# happy.  It is probably not necessary for most users, however.
MikeHappy :
	$(VIRTUAL_ENV)/bin/python setup.py install --prefix=$(VIRTUAL_ENV)
