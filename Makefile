PROGRAM=eventplane

CC=g++
CFLAGS=-g -Wall -fPIC -I$(O2_ROOT)/include `root-config --cflags`
LDFLAGS=-L$(O2_ROOT)/lib -lpthread `root-config --glibs`
INCLUDE=-I$(O2_ROOT)/include/GPU -I$(O2_ROOT)/include/SimulationDataFormat -I$(MS_GSL_ROOT)/include

SOFLAGS=-shared

HEADERS=src/Histos.h \
		src/DataManager.h \
		src/Corrections.h \
		src/Eventplane.h

SOURCES=$(HEADERS:.h=.cxx)
OBJECTS=$(HEADERS:.h=.o)

all: $(PROGRAM)

$(PROGRAM) : $(OBJECTS) $(PROGRAM).cc
	@echo "Linking $@ ..."
	$(CC) -lEG -lPhysics -L$(PWD) $(LDFLAGS) -lO2GPUCommon $(PROGRAM).cc $(CFLAGS) $(OBJECTS) $(INCLUDE) $(FFTWINC) -o $(PROGRAM)
	@echo "Done!"

%.cxx:

%: %.cxx
	$(LINK.cc) $^ $(CFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
	$(COMPILE.cc) -lO2DataFormatsFV0 -lO2DataFormatsFT0 -lO2SimulationDataFormat -lGPUCommon -lGPUCommonRtypes -lgsl -lgslblas $(OUTPUT_OPTION) $(FFTWINC) $(LDFLAGS) $(INCLUDE) $(CFLAGS) $<

clean:
	rm -f $(PROGRAM) src/*.o

epDict.cc: $(HEADERS)
	@echo "Generating dictionary ..."
	@rootcint epDict.cc -c $(HEADERS)