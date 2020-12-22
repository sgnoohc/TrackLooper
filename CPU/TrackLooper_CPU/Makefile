# Simple makefile

EXES=bin/doAnalysis bin/mtv bin/tracklet bin/sdl

SUBDIRS := efficiency
SOURCES=$(wildcard src/*.cc) $(wildcard SDL/*.cc)
OBJECTS=$(SOURCES:.cc=.o)
HEADERS=$(SOURCES:.cc=.h)

CC          = g++
CXX         = g++
CXXFLAGS    = -g -O2 -Wall -fPIC -Wshadow -Woverloaded-virtual
LD          = g++
LDFLAGS     = -g -O2
SOFLAGS     = -g -shared
CXXFLAGS    = -g -O2 -Wall -fPIC -Wshadow -Woverloaded-virtual
LDFLAGS     = -g -O2
ROOTLIBS    = $(shell root-config --libs)
ROOTCFLAGS  = $(shell root-config --cflags)
CXXFLAGS   += $(ROOTCFLAGS) -ISDL -I$(shell pwd) -Isrc
CFLAGS      = $(ROOTCFLAGS) -Wall -Wno-unused-function -g -O2 -fPIC -fno-var-tracking -ISDL -I$(shell pwd) -Isrc
EXTRACFLAGS = $(shell rooutil-config)
EXTRAFLAGS  = -fPIC -ITMultiDrawTreePlayer -Wunused-variable -lTMVA -lEG -lGenVector -lXMLIO -lMLP -lTreePlayer

all: $(EXES) $(SUBDIRS)

bin/doAnalysis: bin/doAnalysis.o $(OBJECTS)
	$(LD) $(CXXFLAGS) $(LDFLAGS) $^ $(ROOTLIBS) $(EXTRAFLAGS) -o $@

bin/mtv: bin/mtv.o $(OBJECTS)
	$(LD) $(CXXFLAGS) $(LDFLAGS) $^ $(ROOTLIBS) $(EXTRAFLAGS) -o $@

bin/tracklet: bin/tracklet.o $(OBJECTS)
	$(LD) $(CXXFLAGS) $(LDFLAGS) $^ $(ROOTLIBS) $(EXTRAFLAGS) -o $@

bin/sdl: bin/sdl.o $(OBJECTS)
	$(LD) $(CXXFLAGS) $(LDFLAGS) $^ $(ROOTLIBS) $(EXTRAFLAGS) -o $@

%.o: %.cc
	$(CC) $(CFLAGS) $(EXTRACFLAGS) $< -c -o $@

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	rm -f $(OBJECTS) bin/*.o $(EXES)
	rm -f efficiency/doAnalysis efficiency/*.o

.PHONY: all $(SUBDIRS)
