CXX = g++
CXXFLAGS = -std=c++0x -Wall -O2
PREFIX ?= /usr/local/bin

OBJ = seq.o search.o stats.o fileroutines.o readroutines.o

all: rm_reads extract remove separate

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

extract: extract.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o extract extract.cpp $(OBJ)

remove: remove.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o remove remove.cpp $(OBJ)

separate: separate.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o separate separate.cpp $(OBJ)

.PHONY: install uninstall clean

clean:
	rm rm_reads extract remove separate
	rm *.o

install:
	cp rm_reads $(PREFIX)/rm_reads
	cp separate $(PREFIX)/separate
	cp extract $(PREFIX)/extract
	cp remove $(PREFIX)/remove
	cp cookiecutter $(PREFIX)/cookiecutter
	chmod +x $(PREFIX)/cookiecutter 

uninstall:
	rm $(PREFIX)/rm_reads $(PREFIX)/extract \
		$(PREFIX)/remove $(PREFIX)/separate \
		$(PREFIX)/cookiecutter

