CC = g++
CFLAGS = -O3 -c -Wall -D CONSIDER_BOTH_ORIENTATION #-D OUTPUT_STRING_PATH 
 #-D OUTPUT_STRING_PATH
DFLAGS = -g -Wall
LDFLAGS = -lz -lm -lpthread

SOURCES_SORT = sort.cpp jreader.cpp
#SOURCES_DEB = debruijin.cpp 
SOURCES = preprocess.cpp bwt.cpp main.cpp ui.cpp

OBJECTS_SORT = $(SOURCES_SORT:.cpp=.o)
#OBJECTS_DEB = $(SOURCES_DEB:.cpp=.o)
OBJECTS = $(SOURCES:.cpp=.o)


EXECUATION_SORT = deSPI-sort
#EXECUATION_DEB = deb
EXECUATION = deSPI 

DEBUG_SORT = deSPI-sort.out
#DEBUG_DEB = deb.out
DEBUG = deSPI.out

all: $(SOURCES_SORT) $(EXECUATION_SORT) $(SOURCES) $(EXECUATION)
	mkdir -p ../bin && cp $(EXECUATION_SORT) $(EXECUATION) ../deSPI-index ../deSPI-download ../bin
$(EXECUATION_SORT): $(OBJECTS_SORT)
	$(CC)  $(OBJECTS_SORT) -o $@ $(LDFLAGS)

#$(EXECUATION_DEB): $(OBJECTS_DEB)

#	$(CC)  $(OBJECTS_DEB) -o $@ $(LDFLAGS)
$(EXECUATION): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
.c.o:
	$(CC) $(CFLAGS) $< -o $@
debug:
	$(CC) $(DFLAGS) $(SOURCES_SORT) -o $(DEBUG_SORT) $(LDFLAGS)

#	$(CC) $(DFLAGS) $(SOURCES_DEB) -D OUTPUT_STRING_PATH -o $(DEBUG_DEB) $(LDFLAGS)
	$(CC) $(DFLAGS) $(SOURCES) -o $(DEBUG) $(LDFLAGS)
clean:
	rm -rf *.o $(EXECUATION_SORT)  $(DEBUG_SORT) $(EXECUATION) $(DEBUG)


 
