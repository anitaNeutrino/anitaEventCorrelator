#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch


### Package subdirectories
LIBDIR=lib
BUILDDIR=build
INCLUDEDIR=include
BINDIR=bin


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

ROOT_LIBRARY = $(LIBDIR)/libAnitaCorrelator.${DLLSUF}
DICT=correlatorDict
LIB_OBJS =  $(addprefix $(BUILDDIR)/, PrettyAnitaEvent.o CorrelationSummary.o CorrelationSummaryAnita3.o UsefulAdu5Pat.o BedmapReader.o RampdemReader.o $(DICT).o )
CLASS_HEADERS =  $(addprefix $(INCLUDEDIR)/, PrettyAnitaEvent.h CorrelationSummary.h CorrelationSummaryAnita3.h UsefulAdu5Pat.h BedmapReader.h RampdemReader.h )

all : $(ROOT_LIBRARY) 

$(DICT).C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c -p -I$(shell $(RC) --incdir) $(INC_ANITA_UTIL) $(CLASS_HEADERS) LinkDef.h

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

$(LIB_OBJS): | $(BUILDDIR) 

$(BINDIR): 
	mkdir -p $(BINDIR)

$(BUILDDIR): 
	mkdir -p $(BUILDDIR)

$(LIBDIR): 
	mkdir -p $(LIBDIR)


#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) | $(LIBDIR)
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $(LIBS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $(LIB_OBJS) -o $@
endif
	@if [ $(shell root-config --version | cut -c1) -ge 6 ]; then \
	cp $(BUILDDIR)/*.pcm $(LIBDIR); \
	fi; # Additional install command for ROOTv6



$(BUILDDIR)/%.$(OBJSUF) : src/%.$(SRCSUF) $(CLASS_HEADERS) Makefile | $(BUILDDIR) 
	@echo "<**Compiling**> "$<
	$(CXX) -I$(INCLUDEDIR) $(CXXFLAGS)  -c $< -o  $@

$(BUILDDIR)/%.$(OBJSUF) : $(BUILDDIR)/%.C
	@echo "<**Compiling**> "$<
	$(CXX) -I$(INCLUDEDIR) -I./ $(CXXFLAGS) -c $< -o  $@


#eventDict.C: $(CLASS_HEADERS)
$(BUILDDIR)/$(DICT).C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c -p -I$(shell $(RC) --incdir) $(INC_ANITA_UTIL) $(CLASS_HEADERS) LinkDef.h


install: $(ROOT_LIBRARY)
ifeq ($(PLATFORM),macosx)
	install -c -m 755 $(ROOT_LIBRARY) $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY)) $(ANITA_UTIL_LIB_DIR)
else
	install -c -m 755 $(ROOT_LIBRARY) $(ANITA_UTIL_LIB_DIR)
endif
	install -c -m 644  $(CLASS_HEADERS) $(ANITA_UTIL_INC_DIR)
	install -d $(ANITA_UTIL_CALIB_DIR)

	@if [ $(shell root-config --version | cut -c1) -ge 6 ]; then \
	install -c -m 755 $(BUILDDIR)/$(DICT)_rdict.pcm $(ANITA_UTIL_LIB_DIR) ;\
	fi # Additional install command for ROOTv6

	for file in data/*.asc data/*.bin data/*.hdr; do install -c -m 644 "$${file}" $(ANITA_UTIL_CALIB_DIR); done

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
	@rm -f makeGoodCorrelationRunTree makeHPolCorrelationRunTree makeCorrelationRunTree makeSimpleDST
