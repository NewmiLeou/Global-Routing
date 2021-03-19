# HW5 Makefile

# declare variable
CC=g++
LDFLAGS=-std=c++11 -O3 -lm
SOURCES=main.cpp
OBJECTS=$(SOURCES:.c=.o)
BINDIR=../bin/
EXE=hw5


all: $(SOURCES) $(EXE)
$(EXE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(BINDIR)$@

clean:
	rm -rf *.o $(BINDIR)$(EXE)