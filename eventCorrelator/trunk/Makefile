
#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#Site Specific  Flags
SYSINCLUDES	= -I/sw/include
SYSLIBS         = -L/sw/lib/ -lgsl 



ifdef ANITA_UTIL_INSTALL_DIR
ANITA_UTIL_LIB_DIR=${ANITA_UTIL_INSTALL_DIR}/lib
ANITA_UTIL_INC_DIR=${ANITA_UTIL_INSTALL_DIR}/include
LD_ANITA_UTIL=-L$(ANITA_UTIL_LIB_DIR)
INC_ANITA_UTIL=-I$(ANITA_UTIL_INC_DIR)
else
ANITA_UTIL_LIB_DIR=
ANITA_UTIL_INC_DIR=
ifdef EVENT_READER_DIR
LD_ANITA_UTIL=-L$(EVENT_READER_DIR)
INC_ANITA_UTIL=-I$(EVENT_READER_DIR)
endif
endif

#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES) $(INC_ANITA_UTIL)
LDFLAGS      += -g $(ROOTLDFLAGS) 
LIBS          = $(ROOTLIBS) -lMathMore -lMinuit $(SYSLIBS) $(LD_ANITA_UTIL) -lAnitaEvent -lRootFftwWrapper -lfftw3
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)


#Now the bits we're actually compiling


#ROOT stuff

ROOT_LIBRARY = libAnitaPlotter.${DLLSUF}
LIB_OBJS = PrettyAnitaEvent.o CorrelationSummary.o UsefulAdu5Pat.o  plotDict.o
CLASS_HEADERS =  PrettyAnitaEvent.h CorrelationSummary.h UsefulAdu5Pat.h

all : $(ROOT_LIBRARY) 

plotDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(CXXFLAGS) $(CLASS_HEADERS) LinkDef.h

progs: makeCorrelationRunTree makeGoodCorrelationRunTree makeInitialGoodCorrelationRunTree


makeCorrelationRunTree : $(ROOT_LIBRARY) makeCorrelationRunTree.$(SRCSUF)
	@echo "<**Compiling**> "  
	$(LD)  $(CXXFLAGS) $(LDFLAGS) makeCorrelationRunTree.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@

makeGoodCorrelationRunTree : $(ROOT_LIBRARY) makeGoodCorrelationRunTree.$(SRCSUF)
	@echo "<**Compiling**> "  
	$(LD)  $(CXXFLAGS) $(LDFLAGS) makeGoodCorrelationRunTree.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@

makeInitialGoodCorrelationRunTree : $(ROOT_LIBRARY) makeInitialGoodCorrelationRunTree.$(SRCSUF)
	@echo "<**Compiling**> "  
	$(LD)  $(CXXFLAGS) $(LDFLAGS) makeInitialGoodCorrelationRunTree.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@


#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	$(LD) $(SOFLAGS) $^ $(LIBS) $(OutPutOpt) $@
ifeq ($(MACOSX_MINOR),4)
	ln -sf $@ $(subst .$(DLLSUF),.so,$@)
else
	$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $(LIBS) $^ \
	 $(OutPutOpt) $(subst .$(DLLSUF),.so,$@)
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIB_OBJS) $(LIBS) -o $@
endif

%.$(OBJSUF) : %.$(SRCSUF)
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


install: $(ROOT_LIBRARY)
ifeq ($(PLATFORM),macosx)
	cp $(ROOT_LIBRARY) $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY)) $(ANITA_UTIL_LIB_DIR)
else
	cp $(ROOT_LIBRARY) $(ANITA_UTIL_LIB_DIR)
endif
	cp  $(CLASS_HEADERS) $(ANITA_UTIL_INC_DIR)

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
