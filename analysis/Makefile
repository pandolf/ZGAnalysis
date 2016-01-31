CC = g++
CFLAGS = -Wall -c -g -fpermissive 

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I.. -I$(ROOTSYS)/include -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src -I$(CMSSW_BASE)/src/CommonTools 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -L$(CMSSW_BASE)/lib/slc6_amd64_gcc491/ -L/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_7/lib/slc6_amd64_gcc491/ -lHtml -lMathCore -lGenVector -lMinuit -lEG  -lPhysics -L/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/zlib/1.2.8-cms/lib -lz -lTreePlayer -lRooFitCore -lRooFit


all: runZGAnalysis

runZGAnalysis: runZGAnalysis.cpp ZGSample.o ZGConfig.o ZGCommonTools.o
	$(CC) -Wall $(INCLUDES) -o runZGAnalysis runZGAnalysis.cpp ZGSample.o ZGConfig.o ZGCommonTools.o  $(ROOTFLAG) $(EXTRALIBS)

convertSNTtrees: convertSNTtrees.cpp
	$(CC) -Wall $(INCLUDES) -o convertSNTtrees convertSNTtrees.cpp $(ROOTFLAG) $(EXTRALIBS)

drawZG: drawZG.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawZG drawZG.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawLimitPlot: drawLimitPlot.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawLimitPlot drawLimitPlot.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawGenEfficiency: drawGenEfficiency.cpp ZGDrawTools.o 
	$(CC) -Wall $(INCLUDES) -o drawGenEfficiency drawGenEfficiency.cpp ZGDrawTools.o $(ROOTFLAG) $(EXTRALIBS)

computeGenEfficiency: computeGenEfficiency.cpp ZGCommonTools.o 
	$(CC) -Wall $(INCLUDES) -o computeGenEfficiency computeGenEfficiency.cpp ZGCommonTools.o  $(ROOTFLAG) $(EXTRALIBS)

computeGenAcceptance: computeGenAcceptance.cpp ZGCommonTools.o ZGDrawTools.o
	$(CC) -Wall $(INCLUDES) -o computeGenAcceptance computeGenAcceptance.cpp ZGCommonTools.o ZGDrawTools.o  $(ROOTFLAG) $(EXTRALIBS)






ZGConfig.o: ../src/ZGConfig.cc 
	$(CC) -Wall $(INCLUDES) -c ../src/ZGConfig.cc $(ROOTFLAG) $(EXTRALIBS)

ZGSample.o: ../src/ZGSample.cc 
	$(CC) -Wall $(INCLUDES) -c ../src/ZGSample.cc $(ROOTFLAG) $(EXTRALIBS)

ZGDrawTools.o: ../src/ZGDrawTools.cc 
	$(CC) -Wall $(INCLUDES) -c ../src/ZGDrawTools.cc $(ROOTFLAG) $(EXTRALIBS)

ZGCommonTools.o: ../src/ZGCommonTools.cc 
	$(CC) -Wall $(INCLUDES) -c ../src/ZGCommonTools.cc $(ROOTFLAG) $(EXTRALIBS)







clean:
	rm *.o 

