CXX= g++
CXXFLAGS = -std=c++0x -Wall
OPT = -O2
DEBUG = -g -O0 -D DEBUG
all: *.cpp
	$(CXX) $(CXXFLAGS) $(OPT) search.cpp seq.cpp stats.cpp rm_reads.cpp -o rm_reads
	$(CXX) $(CXXFLAGS) $(OPT) search.cpp seq.cpp stats.cpp extractor.cpp -o extractor
	$(CXX) $(CXXFLAGS) $(OPT) search.cpp seq.cpp stats.cpp remove.cpp -o remove
.PHONY: clean
clean:
	rm -rf rm_reads

