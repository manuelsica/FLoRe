CXX = g++
CXXFLAGS = -std=c++17 -Wall -O3

OBJS = main.o logging.o read.o util.o index.o overlap.o jsonoutput.o

all: FLORE_bin set-exec

FLORE_bin: $(OBJS)
	$(CXX) $(CXXFLAGS) -o FLORE_bin $(OBJS)

main.o: main.cpp logging.hpp read.hpp util.hpp index.hpp overlap.hpp jsonoutput.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp

logging.o: logging.cpp logging.hpp
	$(CXX) $(CXXFLAGS) -c logging.cpp

read.o: read.cpp read.hpp
	$(CXX) $(CXXFLAGS) -c read.cpp

util.o: util.cpp util.hpp
	$(CXX) $(CXXFLAGS) -c util.cpp

index.o: index.cpp index.hpp util.hpp
	$(CXX) $(CXXFLAGS) -c index.cpp

overlap.o: overlap.cpp overlap.hpp index.hpp
	$(CXX) $(CXXFLAGS) -c overlap.cpp

jsonoutput.o: jsonoutput.cpp jsonoutput.hpp
	$(CXX) $(CXXFLAGS) -c jsonoutput.cpp

set-exec:
	chmod +x FLORE_bin 2>/dev/null || true

clean:
	rm -f $(OBJS) FLORE_bin