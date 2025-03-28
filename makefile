CXX = g++
CXXFLAGS = -std=c++17 -Wall -O3

OBJS = main.o util.o overlap.o read.o

all: FLORE_bin set-exec

FLORE_bin: $(OBJS)
	$(CXX) $(CXXFLAGS) -o FLORE_bin $(OBJS)

main.o: main.cpp util.hpp overlap.hpp read.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp

util.o: util.cpp util.hpp
	$(CXX) $(CXXFLAGS) -c util.cpp

overlap.o: overlap.cpp overlap.hpp
	$(CXX) $(CXXFLAGS) -c overlap.cpp

read.o: read.cpp read.hpp
	$(CXX) $(CXXFLAGS) -c read.cpp

set-exec:
	chmod +x FLORE 2>/dev/null || true

clean:
	rm -f $(OBJS) FLORE_bin
