#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#Site Specific  Flags
SYSINCLUDES	= 
SYSLIBS         = #-lprofiler -ltcmalloc
DLLSUF = ${DllSuf}
OBJSUF = ${ObjSuf}
SRCSUF = ${SrcSuf}


ifdef ANITA_UTIL_INSTALL_DIR
ANITA_UTIL_LIB_DIR=${ANITA_UTIL_INSTALL_DIR}/lib
ANITA_UTIL_INC_DIR=${ANITA_UTIL_INSTALL_DIR}/include
LD_ANITA_UTIL=-L$(ANITA_UTIL_LIB_DIR) -lAnitaEvent
INC_ANITA_UTIL=-I$(ANITA_UTIL_INC_DIR)
ANITA_UTIL_CALIB_DIR=$(ANITA_UTIL_INSTALL_DIR)/share/anitaCalib
else
ANITA_UTIL_LIB_DIR=/usr/local/lib
ANITA_UTIL_INC_DIR=/usr/local/include
ANITA_UTIL_CALIB_DIR=/usr/local/share/anitaCalib
ifdef EVENT_READER_DIR
LD_ANITA_UTIL=-L$(EVENT_READER_DIR)  -lAnitaEvent
INC_ANITA_UTIL=-I$(EVENT_READER_DIR)
endif
endif


#Toggles the FFT functions on and off
USE_FFT_TOOLS=1

ifdef USE_FFT_TOOLS
FFTLIBS = -lRootFftwWrapper -lfftw3
FFTFLAG = -DUSE_FFT_TOOLS
else
FFTLIBS =
FFTFLAG =
endif

#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS) $(FFTFLAG) $(SYSINCLUDES) $(INC_ANITA_UTIL)
LDFLAGS      += -g $(ROOTLDFLAGS) 

LIBS          = $(ROOTLIBS) -lMathMore -lMinuit $(SYSLIBS) $(LD_ANITA_UTIL) $(FFTLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#Now the bits we're actually compiling


#ROOT stuff

ROOT_LIBRARY = libAnitaCorrelator.${DLLSUF}
LIB_OBJS = PrettyAnitaEvent.o CorrelationSummary.o UsefulAdu5Pat.o  correlatorDict.o BedmapReader.o RampdemReader.o
CLASS_HEADERS =  PrettyAnitaEvent.h CorrelationSummary.h UsefulAdu5Pat.h BedmapReader.h RampdemReader.h

all : $(ROOT_LIBRARY) 

correlatorDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(CXXFLAGS) $(CLASS_HEADERS) LinkDef.h

progs: makeCorrelationRunTree makeSimpleDST #makeGoodCorrelationRunTree makeInitialGoodCorrelationRunTree  makeHPolCorrelationRunTree


makeCorrelationRunTree : $(ROOT_LIBRARY) makeCorrelationRunTree.$(SRCSUF)
	@echo "<**Compiling**> "  
	$(LD)  $(CXXFLAGS) $(LDFLAGS) makeCorrelationRunTree.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@

makeSimpleDST : $(ROOT_LIBRARY) makeSimpleDST.$(SRCSUF)
	@echo "<**Compiling**> "  
	$(LD)  $(CXXFLAGS) $(LDFLAGS) makeSimpleDST.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@

makeGoodCorrelationRunTree : $(ROOT_LIBRARY) makeGoodCorrelationRunTree.$(SRCSUF)
	@echo "<**Compiling**> "  
	$(LD)  $(CXXFLAGS) $(LDFLAGS) makeGoodCorrelationRunTree.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@

makeHPolCorrelationRunTree : $(ROOT_LIBRARY) makeHPolCorrelationRunTree.$(SRCSUF)
	@echo "<**Compiling**> "  
	$(LD)  $(CXXFLAGS) $(LDFLAGS) makeHPolCorrelationRunTree.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@


#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIB_OBJS) -o $@
endif

%.$(OBJSUF) : %.$(SRCSUF)
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


install: $(ROOT_LIBRARY)
ifeq ($(PLATFORM),macosx)
	install -c -m 755 $(ROOT_LIBRARY) $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY)) $(ANITA_UTIL_LIB_DIR)
else
	install -c -m 755 $(ROOT_LIBRARY) $(ANITA_UTIL_LIB_DIR)
endif
	install -c -m 644  $(CLASS_HEADERS) $(ANITA_UTIL_INC_DIR)
	install -d $(ANITA_UTIL_CALIB_DIR)
	for file in data/*.asc data/*.bin data/*.hdr; do install -c -m 644 "$${file}" $(ANITA_UTIL_CALIB_DIR); done

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
	@rm -f makeGoodCorrelationRunTree makeHPolCorrelationRunTree makeCorrelationRunTree makeSimpleDST
