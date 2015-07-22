CXX = g++
CXXFLAGS = -std=c++0x -Wall -O2
PREFIX ?= /usr/local

OBJ = seq.o search.o stats.o fileroutines.o readroutines.o

all: rm_reads extractor remove

fileroutines.o: fileroutines.h fileroutines.cpp
	$(CXX) -c $(CXXFLAGS) fileroutines.cpp

readroutines.o: readroutines.h readroutines.cpp
	$(CXX) -c $(CXXFLAGS) readroutines.cpp

seq.o: seq.h seq.cpp
	$(CXX) -c $(CXXFLAGS) seq.cpp

search.o: search.h search.cpp
	$(CXX) -c $(CXXFLAGS) search.cpp

stats.o: stats.h stats.cpp
	$(CXX) -c $(CXXFLAGS) stats.cpp

rm_reads: rm_reads.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o rm_reads rm_reads.cpp $(OBJ)

extractor: extractor.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o extractor extractor.cpp $(OBJ)

remove: remove.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o remove remove.cpp $(OBJ)

separate: separate.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o separate separate.cpp $(OBJ)

.PHONY: install uninstall clean

clean:
	rm rm_reads extractor remove
	rm *.o

install:
	cp rm_reads $(PREFIX)/bin/rm_reads
	cp extractor $(PREFIX)/bin/extractor
	cp remove $(PREFIX)/bin/remove

uninstall:
	rm $(PREFIX)/bin/rm_reads $(PREFIX)/bin/extractor $(PREFIX)/bin/remove

