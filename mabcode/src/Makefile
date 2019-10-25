SHELL = /bin/sh

CPP = g++
LINK := g++ 

INCLUDE = -I.
CFLAGS += -std=c++11
CFLAGS += -O3


all: mab
build: mab

debug: CFLAGS +=-DDEBUG -g -G
debug: ssppr

OBJ = $(patsubst %.cpp,%.o,$(wildcard *.cpp))

%.o: %.cpp 
	$(CPP) $< $(CFLAGS) $(INCLUDE) -c -o $@

%.o: %.cu
	$(CPP) $< $(CFLAGS) $(INCLUDE) -c -o $@
	
mab: $(OBJ)
	$(CPP) $(CFLAGS) $(INCLUDE) -o $@ $^ 

clean:
	rm -rf mab *.o
	