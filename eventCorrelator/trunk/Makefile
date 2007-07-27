
#############################################################################
####


#Global Thingies
CC 	= g++
LD	= g++
SOFLAGS	= -shared
OBJSUF	= o
SRCSUF	= cxx
OPT	= -O2 -g --debug  -Wall -fPIC

ROOTINCLUDES      = -I$(ROOTSYS)/include 
INCLUDES	= -I$(EVENT_READER_DIR) 
CXXFLAGS	= $(EXCEPTION) $(OPT) $(CXXOPT) -fPIC $(INCLUDES) $(ROOTINCLUDES)

ROOTLIBS      = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs) -lMinuit -lTreePlayer -lMathMore
ROOTGLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)
LIBS		= -lz -lm $(ROOTLIBS) -lfftw3



#ROOT stuff

ROOT_LIBRARY = libAnitaPlotter.so 
LIB_OBJS = PrettyAnitaEvent.o FFTtools.o FFTWComplex.o CorrelationSummary.o plotDict.o
CLASS_HEADERS =  PrettyAnitaEvent.h

all : $(ROOT_LIBRARY) 

plotDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(INCLUDES) $(CLASS_HEADERS) LinkDef.h



$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIB_OBJS) -o $@

%.$(OBJSUF) : %.$(SRCSUF)
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(TEST)
