# Include your header file location
CXXFLAGS += $(shell larlite-config --includes)
CXXFLAGS += -I. $(shell root-config --cflags) -g

#include your shared object lib location
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += -L$(LARLITE_LIBDIR) -lAdrien_Chimera
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -g

# platform-specific options
OSNAME = $(shell uname -s)
include $(LARLITE_BASEDIR)/Makefile/Makefile.${OSNAME}

#add your program below with a space after the previous one.
#This makefile compiles all binaries specified below.
PROGRAMS = GenerateChimeraEvents

all= 	$(PROGRAMS)

$(PROGRAMS):
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cc -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM

clean:
	rm -f $(PROGRAMS)