CC = g++
#CFLAGS = #-O3i
CFLAGS = -c -Wall -O3#-D CONSIDER_BOTH_ORIENTATION #-D OUTPUT_STRING_PATH 
 #-D OUTPUT_STRING_PATH
DFLAGS = -g -Wall
LDFLAGS = -lz -lm -lpthread

#SOURCES_SORT = sort.cpp jreader.cpp
#SOURCES_DEB = debruijin.cpp 
SOURCES = sort.cpp bwt.cpp main.cpp ui.cpp jreader.cpp special_kmers.cpp

#OBJECTS_SORT = $(SOURCES_SORT:.cpp=.o)
#OBJECTS_DEB = $(SOURCES_DEB:.cpp=.o)
OBJECTS = $(SOURCES:.cpp=.o)


#EXECUATION_SORT = deSPI-sort
#EXECUATION_DEB = deb
EXECUATION = deSPI 

#DEBUG_SORT = deSPI-sort.out
#DEBUG_DEB = deb.out
DEBUG = deSPI.out

all:  $(SOURCES) $(EXECUATION)
	mkdir -p ../bin && cp $(EXECUATION) ../scripts/* ../third_party/jdbg ../bin  
$(EXECUATION): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
.c.o:
	$(CC) $(CFLAGS) $< -o $@
debug:
	$(CC) $(DFLAGS) $(SOURCES) -o $(DEBUG) $(LDFLAGS)
clean:
	rm -rf *.o    $(EXECUATION) $(DEBUG) ../bin/*


 
