CC=g++
CFLAGS=-c -Wall
LDFLAGS=
ROOTFLAGS=`root-config --cflags --glibs --libs --evelibs`
SOURCES=test.cpp EGammaSF.cc 
BINARIES=test.o EGammaSF.o
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(ROOTFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(ROOTFLAGS) $< -o $@

clean:
	rm -f $(BINARIES) 