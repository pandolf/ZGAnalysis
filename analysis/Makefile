CC = g++
CFLAGS = -Wall -c -g -fpermissive 

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I.. -I$(ROOTSYS)/include -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src -I$(CMSSW_BASE)/src/CommonTools 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -L$(CMSSW_BASE)/lib/slc6_amd64_gcc491/ -L/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_7/lib/slc6_amd64_gcc491/ -lHtml -lMathCore -lGenVector -lMinuit -lEG  -lPhysics -L/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/zlib/1.2.8-cms/lib -lz -lTreePlayer -lRooFitCore -lRooFit


all: runZGAnalysis

runZGAnalysis: runZGAnalysis.cpp ZGSample.o ZGConfig.o ZGCommonTools.o RoccoR.o rochcor2016.o EnergyScaleCorrection_class.o
	$(CC) -Wall $(INCLUDES) -o runZGAnalysis runZGAnalysis.cpp ZGSample.o ZGConfig.o ZGCommonTools.o RoccoR.o rochcor2016.o EnergyScaleCorrection_class.o $(ROOTFLAG) $(EXTRALIBS)

convertSNTtrees: convertSNTtrees.cpp
	$(CC) -Wall $(INCLUDES) -o convertSNTtrees convertSNTtrees.cpp $(ROOTFLAG) $(EXTRALIBS)

optimizeCut: optimizeCut.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o optimizeCut optimizeCut.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

stitchZGMC: stitchZGMC.cpp 
	$(CC) -Wall $(INCLUDES) -o stitchZGMC stitchZGMC.cpp $(ROOTFLAG) $(EXTRALIBS)

compareStitching: compareStitching.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o compareStitching compareStitching.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawOptimizeCut: drawOptimizeCut.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawOptimizeCut drawOptimizeCut.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

compareMCsignals: compareMCsignals.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o compareMCsignals compareMCsignals.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

compareSmearing: compareSmearing.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o compareSmearing compareSmearing.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)


drawZG: drawZG.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawZG drawZG.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawQ: drawQ.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawQ drawQ.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawSignals: drawSignals.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawSignals drawSignals.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawGenSignals: drawGenSignals.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawGenSignals drawGenSignals.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawLimitPlot: drawLimitPlot.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawLimitPlot drawLimitPlot.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawLimitPlot_8TeV: drawLimitPlot_8TeV.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawLimitPlot_8TeV drawLimitPlot_8TeV.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawLimitPlotQgamma: drawLimitPlotQgamma.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawLimitPlotQgamma drawLimitPlotQgamma.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawCombinedLimitPlot: drawCombinedLimitPlot.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawCombinedLimitPlot drawCombinedLimitPlot.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawCombinedLimitPlot_all: drawCombinedLimitPlot_all.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawCombinedLimitPlot_all drawCombinedLimitPlot_all.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawCombinedPvaluePlot: drawCombinedPvaluePlot.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawCombinedPvaluePlot drawCombinedPvaluePlot.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawCombinedPvaluePlot_all: drawCombinedPvaluePlot_all.cpp ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawCombinedPvaluePlot_all drawCombinedPvaluePlot_all.cpp ZGDrawTools.o ZGConfig.o $(ROOTFLAG) $(EXTRALIBS)

drawGenEfficiency: drawGenEfficiency.cpp ZGDrawTools.o 
	$(CC) -Wall $(INCLUDES) -o drawGenEfficiency drawGenEfficiency.cpp ZGDrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawAcceptanceTimesEfficiency: drawAcceptanceTimesEfficiency.cpp ZGDrawTools.o 
	$(CC) -Wall $(INCLUDES) -o drawAcceptanceTimesEfficiency drawAcceptanceTimesEfficiency.cpp ZGDrawTools.o $(ROOTFLAG) $(EXTRALIBS)

computeGenEfficiency: computeGenEfficiency.cpp ZGCommonTools.o muresolution_run2.o rochcor2015.o
	$(CC) -Wall $(INCLUDES) -o computeGenEfficiency computeGenEfficiency.cpp ZGCommonTools.o muresolution_run2.o rochcor2015.o $(ROOTFLAG) $(EXTRALIBS)

computeGenAcceptance: computeGenAcceptance.cpp ZGCommonTools.o ZGDrawTools.o
	$(CC) -Wall $(INCLUDES) -o computeGenAcceptance computeGenAcceptance.cpp ZGCommonTools.o ZGDrawTools.o  $(ROOTFLAG) $(EXTRALIBS)

compareMassShape: compareMassShape.cpp ZGCommonTools.o ZGDrawTools.o
	$(CC) -Wall $(INCLUDES) -o compareMassShape compareMassShape.cpp ZGCommonTools.o ZGDrawTools.o  $(ROOTFLAG) $(EXTRALIBS)

generateSignalToys: generateSignalToys.cpp ZGCommonTools.o ZGDrawTools.o
	$(CC) -Wall $(INCLUDES) -o generateSignalToys generateSignalToys.cpp ZGCommonTools.o ZGDrawTools.o  $(ROOTFLAG) $(EXTRALIBS)

fitGenSignalShapes: fitGenSignalShapes.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o RooDoubleCBShape.o
	$(CC) -Wall $(INCLUDES) -o fitGenSignalShapes fitGenSignalShapes.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o RooDoubleCBShape.o $(ROOTFLAG) $(EXTRALIBS)

fitSignalShapes: fitSignalShapes.cpp ZGCommonTools.o ZGDrawTools.o RooDoubleCBShape.o ZGConfig.o 
	$(CC) -Wall $(INCLUDES) -o fitSignalShapes fitSignalShapes.cpp ZGCommonTools.o ZGDrawTools.o RooDoubleCBShape.o ZGConfig.o  $(ROOTFLAG) $(EXTRALIBS)

computeSyst: computeSyst.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o computeSyst computeSyst.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o  $(ROOTFLAG) $(EXTRALIBS)

drawCutEfficiency: drawCutEfficiency.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o
	$(CC) -Wall $(INCLUDES) -o drawCutEfficiency drawCutEfficiency.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o  $(ROOTFLAG) $(EXTRALIBS)

drawSignalEfficiency: drawSignalEfficiency.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o
	$(CC) -Wall $(INCLUDES) -o drawSignalEfficiency drawSignalEfficiency.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o  $(ROOTFLAG) $(EXTRALIBS)

drawSignalEfficiency_convVeto: drawSignalEfficiency_convVeto.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o
	$(CC) -Wall $(INCLUDES) -o drawSignalEfficiency_convVeto drawSignalEfficiency_convVeto.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o  $(ROOTFLAG) $(EXTRALIBS)

draw74vs76: draw74vs76.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o
	$(CC) -Wall $(INCLUDES) -o draw74vs76 draw74vs76.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o  $(ROOTFLAG) $(EXTRALIBS)

drawGenSignalShapes: drawGenSignalShapes.cpp ZGDrawTools.o 
	$(CC) -Wall $(INCLUDES) -o drawGenSignalShapes drawGenSignalShapes.cpp ZGDrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawTriggerTurnOn: drawTriggerTurnOn.cpp ZGDrawTools.o 
	$(CC) -Wall $(INCLUDES) -o drawTriggerTurnOn drawTriggerTurnOn.cpp ZGDrawTools.o $(ROOTFLAG) $(EXTRALIBS)

runFTest: runFTest.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o
	$(CC) -Wall $(INCLUDES) -o runFTest runFTest.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o  $(ROOTFLAG) $(EXTRALIBS)

drawBiasFunction: drawBiasFunction.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o
	$(CC) -Wall $(INCLUDES) -o drawBiasFunction drawBiasFunction.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o  $(ROOTFLAG) $(EXTRALIBS)

drawPrivateVsOfficial: drawPrivateVsOfficial.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o
	$(CC) -Wall $(INCLUDES) -o drawPrivateVsOfficial drawPrivateVsOfficial.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o  $(ROOTFLAG) $(EXTRALIBS)

drawMiniIsoEfficiency: drawMiniIsoEfficiency.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o
	$(CC) -Wall $(INCLUDES) -o drawMiniIsoEfficiency drawMiniIsoEfficiency.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o  $(ROOTFLAG) $(EXTRALIBS)

computeMuonScaleSyst: computeMuonScaleSyst.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o rochcor2015.o muresolution_run2.o
	$(CC) -Wall $(INCLUDES) -o computeMuonScaleSyst computeMuonScaleSyst.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o rochcor2015.o muresolution_run2.o $(ROOTFLAG) $(EXTRALIBS)

drawMuonScaleSyst: drawMuonScaleSyst.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o rochcor2015.o muresolution_run2.o
	$(CC) -Wall $(INCLUDES) -o drawMuonScaleSyst drawMuonScaleSyst.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o rochcor2015.o muresolution_run2.o $(ROOTFLAG) $(EXTRALIBS)

drawSignalHLTEfficiency: drawSignalHLTEfficiency.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o rochcor2015.o muresolution_run2.o
	$(CC) -Wall $(INCLUDES) -o drawSignalHLTEfficiency drawSignalHLTEfficiency.cpp ZGCommonTools.o ZGDrawTools.o ZGConfig.o ZGSample.o rochcor2015.o muresolution_run2.o $(ROOTFLAG) $(EXTRALIBS)



provaroofit: provaroofit.C 
	$(CC) -Wall $(INCLUDES) -o provaroofit provaroofit.C  $(ROOTFLAG) $(EXTRALIBS)





ZGConfig.o: ../src/ZGConfig.cc 
	$(CC) -Wall $(INCLUDES) -c ../src/ZGConfig.cc $(ROOTFLAG) $(EXTRALIBS)

ZGSample.o: ../src/ZGSample.cc 
	$(CC) -Wall $(INCLUDES) -c ../src/ZGSample.cc $(ROOTFLAG) $(EXTRALIBS)

ZGDrawTools.o: ../src/ZGDrawTools.cc 
	$(CC) -Wall $(INCLUDES) -c ../src/ZGDrawTools.cc $(ROOTFLAG) $(EXTRALIBS)

ZGCommonTools.o: ../src/ZGCommonTools.cc 
	$(CC) -Wall $(INCLUDES) -c ../src/ZGCommonTools.cc $(ROOTFLAG) $(EXTRALIBS)

RooDoubleCBShape.o: ../src/RooDoubleCBShape.cc 
	$(CC) -Wall $(INCLUDES) -c ../src/RooDoubleCBShape.cc $(ROOTFLAG) $(EXTRALIBS)




muresolution_run2.o: ../src/muresolution_run2.cc
	$(CC) -Wall $(INCLUDES) -c ../src/muresolution_run2.cc $(ROOTFLAG) $(EXTRALIBS)

rochcor2015.o: ../src/rochcor2015.cc
	$(CC) -Wall $(INCLUDES) -c ../src/rochcor2015.cc $(ROOTFLAG) $(EXTRALIBS)

rochcor2016.o: ../src/rochcor2016.cc
	$(CC) -Wall $(INCLUDES) -c ../src/rochcor2016.cc $(ROOTFLAG) $(EXTRALIBS)

RoccoR.o: ../src/RoccoR.cc
	$(CC) -Wall $(INCLUDES) -c ../src/RoccoR.cc $(ROOTFLAG) $(EXTRALIBS)

EnergyScaleCorrection_class.o: ../src/EnergyScaleCorrection_class.cc
	$(CC) -Wall $(INCLUDES) -c ../src/EnergyScaleCorrection_class.cc $(ROOTFLAG) $(EXTRALIBS)


clean:
	rm *.o 

