IDIR       = .
ODIR       = obj
SDIR       = .
RSDIR      = $(CMSSW_BASE)/src/SusyAnaTools/Tools

TTIDIR     = $(CMSSW_BASE)/src/TopTagger/TopTagger/include
TPIDIR     = $(CMSSW_BASE)/src/TopTagger/CfgParser/include

OPENCV_DIRECTORY=$(CMSSW_BASE)/src/opencv

CXX        = g++

CXXFLAGS  += -I. -I$(CMSSW_BASE)/src -I$(OPENCV_DIRECTORY)/include/ -I$(OPENCV_DIRECTORY)/modules/core/include -I$(OPENCV_DIRECTORY)/modules/video/include -I$(OPENCV_DIRECTORY)/modules/objdetect/include -I$(OPENCV_DIRECTORY)/modules/ml/include/ -I$(OPENCV_DIRECTORY)/modules/photo/include/ -I$(OPENCV_DIRECTORY)/modules/imgproc/include/ -std=c++0x

## Optimization flag
CXXFLAGS += -g #-O3
## Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

## Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD         = g++
LDFLAGS    =

ROOTLIBS   = $(shell root-config --glibs)
MT2LIB     = -L$(CMSSW_BASE)/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2
TOPLIB     = -L$(CMSSW_BASE)/lib/${SCRAM_ARCH}/ -lTopTaggerTopTagger -lTopTaggerCfgParser
OPENVCLIBS = -L$(OPENCV_DIRECTORY)/lib/ -lopencv_ml -lopencv_core

#OBJS       = $(patsubst %, $(ODIR)/%, $(OBJ))

PROGRAMS = DataMC mkNjetsWeightHist

all: mkobj sampPyWrap $(PROGRAMS)

mkobj:
	@mkdir -p obj

#code to compile shared library to link samples to python                                                                                                                               
sampPyWrap: $(ODIR)/samplesModule.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(RSDIR)/samples.cc $(RSDIR)/samples.h
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(RSDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<


$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(SDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(SDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(SDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RSDIR) -I$(TTIDIR) -I$(TPIDIR) -o $@ -c $<

DataMC: $(ODIR)/DataMC.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/samples.o $(ODIR)/baselineDef.o $(ODIR)/customize.o $(ODIR)/searchBins.o
	$(LD) $^ $(TOPLIB) $(OPENVCLIBS) $(MT2LIB) $(ROOTLIBS) -o $@

mkNjetsWeightHist: $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/samples.o $(ODIR)/baselineDef.o $(ODIR)/customize.o $(ODIR)/mkNjetsWeightHist.o
	$(LD) $^ $(TOPLIB) $(OPENVCLIBS) $(MT2LIB) $(ROOTLIBS) -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.so $(ODIR)/*.d $(PROGRAMS) core 

-include $(ODIR)/*.d
