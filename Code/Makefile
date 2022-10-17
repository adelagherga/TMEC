# Makefile for programs for Bennett-Gherga-Rechnitzer algorithm

# eclib is a requirement.  Define ECLIB_BASE to be the installation
# directory, for example (if installed in the usual place) /usr/local,
# or if installed anywhere else set the install directory here:

ECLIB_BASE=/usr/local
#ECLIB_BASE=$(HOME)/eclib
INCDIR = $(ECLIB_BASE)/include
LIBDIR = $(ECLIB_BASE)/lib

GCC=g++ -std=c++11 -fmax-errors=1
CC = $(GCC)

# NB If used with a multithreaded build of eclib then you MUST define
# USE_BOOST=1 below so that the correct compiler and linker stuff is
# appended below.  Otherwise set USE_BOOST=0 (or do not set it).

USE_BOOST=1
ifeq ($(USE_BOOST), 1)
 BOOST_ASIO_LIB = -lboost_system-mt
 BOOST_CPPFLAGS =   -DECLIB_MULTITHREAD -DHAVE_STDCXX_0X=/\*\*/ -DHAVE_TR1_UNORDERED_MAP=/\*\*/ -DHAVE_STDCXX_0X=/\*\*/ -DHAVE_UNORDERED_MAP=/\*\*/# -pthread -I/usr/include
 BOOST_SYSTEM_LIB = -lboost_system
 BOOST_THREAD_LIB = -lboost_thread
 BOOST_LDFLAGS = -L/usr/lib -pthread $(BOOST_SYSTEM_LIB) $(BOOST_THREAD_LIB)
endif

# to disable checking of assert() use the following:
#OPTFLAG = -DNDEBUG -O3 -fPIC
# to enable checking of assert() use the following:
OPTFLAG = -O3 -Wall -fPIC

CFLAGS = -c -g $(OPTFLAG) $(BOOST_CPPFLAGS) -I$(INCDIR)
LFLAGS = -lec -lntl -lstdc++  -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)

all: tests

sources: ccs headers
	chmod a+r *.h *.cc

ccs: cubic_utils.cc TME.cc N2TME.cc TME_compare.cc egros.cc CurvesNj0.cc CurvesNj1728.cc # S2TME.cc

headers: cubic_utils.h TME.h egros.h

%.o:   %.cc
	$(CC) $(CFLAGS) $<

TESTS = N2TME TME_compare CurvesNj0 CurvesNj1728 # S2TME
tests: $(TESTS)

clean:
	rm -f $(TESTS)
	rm -f *.o *~

N2TME: N2TME.o cubic_utils.o TME.o
	$(CC) -o N2TME N2TME.o cubic_utils.o TME.o $(LFLAGS)

S2TME: S2TME.o cubic_utils.o
	$(CC) -o S2TME S2TME.o cubic_utils.o $(LFLAGS)

TME_compare: TME_compare.o cubic_utils.o TME.o
	$(CC) -o TME_compare TME_compare.o cubic_utils.o TME.o $(LFLAGS)

CurvesNj0: CurvesNj0.o egros.o
	$(CC) -o CurvesNj0 CurvesNj0.o egros.o $(LFLAGS)

CurvesNj1728: CurvesNj1728.o egros.o
	$(CC) -o CurvesNj1728 CurvesNj1728.o egros.o $(LFLAGS)


CurvesNj1728.o: CurvesNj1728.cc egros.h
CurvesNj0.o: CurvesNj0.cc egros.h
N2TME.o: N2TME.cc cubic_utils.h TME.h
S2TME.o: S2TME.cc cubic_utils.h
TME_compare.o: TME_compare.cc cubic_utils.h TME.h
cubic_utils.o: cubic_utils.h cubic_utils.cc
TME.o: cubic_utils.h TME.h TME.cc

