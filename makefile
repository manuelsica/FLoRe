# Makefile per compilare il tool FLoRe con supporto CFL
# Utilizza g++ con standard C++17, abilitando -Wall e -O3 per ottimizzazione.

CXX       = g++
CXXFLAGS  = -std=c++17 -Wall -O3
OBJS = main.o logging.o read.o util.o cfl.o icfl.o cfl_icfl.o filter.o index.o overlap.o jsonoutput.o profiling.o

all: FLORE_bin set-exec

FLORE_bin: $(OBJS)
	$(CXX) $(CXXFLAGS) -o FLORE_bin $(OBJS)

main.o: main.cpp logging.hpp read.hpp util.hpp filter.hpp index.hpp overlap.hpp jsonoutput.hpp profiling.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp

logging.o: logging.cpp logging.hpp
	$(CXX) $(CXXFLAGS) -c logging.cpp

read.o: read.cpp read.hpp
	$(CXX) $(CXXFLAGS) -c read.cpp

util.o: util.cpp util.hpp cfl.hpp icfl.hpp cfl_icfl.hpp
	$(CXX) $(CXXFLAGS) -c util.cpp

cfl.o: cfl.cpp cfl.hpp
	$(CXX) $(CXXFLAGS) -c cfl.cpp

icfl.o: icfl.cpp icfl.hpp
	$(CXX) $(CXXFLAGS) -c icfl.cpp

# nuovo modulo CFL→ICFL
cfl_icfl.o: cfl_icfl.cpp cfl_icfl.hpp cfl.hpp icfl.hpp
	$(CXX) $(CXXFLAGS) -c cfl_icfl.cpp

# Regola per il modulo filter
filter.o: filter.cpp filter.hpp
	$(CXX) $(CXXFLAGS) -c filter.cpp

index.o: index.cpp index.hpp util.hpp cfl.hpp icfl.hpp filter.hpp overlap.hpp
	$(CXX) $(CXXFLAGS) -c index.cpp

overlap.o: overlap.cpp overlap.hpp index.hpp
	$(CXX) $(CXXFLAGS) -c overlap.cpp

jsonoutput.o: jsonoutput.cpp jsonoutput.hpp
	$(CXX) $(CXXFLAGS) -c jsonoutput.cpp

profiling.o: profiling.cpp profiling.hpp
	$(CXX) $(CXXFLAGS) -c profiling.cpp

set-exec:
	chmod +x FLORE_bin 2>/dev/null || true

clean:
	rm -f $(OBJS) FLORE_bin