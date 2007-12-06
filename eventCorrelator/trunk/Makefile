
#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#Site Specific  Flags
SYSINCLUDES	= 
SYSLIBS         = 
OUTLIBDIR          = /home/rjn/install/lib
OUTINCDIR          = /home/rjn/install/includes

#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES) -I$(EVENT_READER_DIR)
LDFLAGS      += $(ROOTLDFLAGS) -L$(EVENT_READER_DIR)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
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



#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifeq ($(MACOSX_MINOR),4)
	ln -sf $@ $(subst .$(DLLSUF),.so,$@)
else
	$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
	 $(OutPutOpt) $(subst .$(DLLSUF),.so,$@)
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
	@cp $(ROOT_LIBRARY) $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY)) $(OUTLIBDIR)
else
	@cp $(ROOT_LIBRARY) $(OUTLIBDIR)
endif
	@cp  $(CLASS_HEADERS) $(OUTINCDIR)

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
