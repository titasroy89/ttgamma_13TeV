all: makeTemplates makeSkim signalAcceptance

makeTemplates: Histogrammer.o HistCollect.o Selector.o EventPick.o EventTree.o JECvariation.o BTagCalibrationStandalone.o makeTemplates.cpp OverlapRemove.cpp PUReweight.cpp PUReweight.h  JECvariation.cpp BTagCalibrationStandalone.cpp 
	g++ -o makeTemplates `root-config --libs` -I`root-config --incdir` -std=c++11 EventTree.o Selector.o EventPick.o Histogrammer.o HistCollect.o JECvariation.o  BTagCalibrationStandalone.o makeTemplates.cpp OverlapRemove.cpp  BTagCalibrationStandalone.cpp PUReweight.cpp JetMETObjects/FactorizedJetCorrector.o JetMETObjects/JetCorrectorParameters.o JetMETObjects/SimpleJetCorrector.o JetMETObjects/JetCorrectionUncertainty.o JetMETObjects/SimpleJetCorrectionUncertainty.o

makeSkim: Selector.o EventPick.o EventTree.o JECvariation.o makeSkim.cpp OverlapRemove.cpp
	g++ -o makeSkim `root-config --libs` -I`root-config --incdir` -std=c++11 EventTree.o EventPick.o Selector.o JECvariation.o OverlapRemove.cpp makeSkim.cpp JetMETObjects/FactorizedJetCorrector.o JetMETObjects/JetCorrectorParameters.o JetMETObjects/SimpleJetCorrector.o JetMETObjects/JetCorrectionUncertainty.o JetMETObjects/SimpleJetCorrectionUncertainty.o 
signalAcceptance: Selector.o EventPick.o EventTree.o signalAcceptance.cpp OverlapRemove.cpp
	g++ -o signalAcceptance `root-config --libs` -I`root-config --incdir` -std=c++11 EventTree.o EventPick.o Selector.o OverlapRemove.cpp signalAcceptance.cpp

EventTree.o: EventTree.cpp EventTree.h
	g++ -c -I`root-config --incdir` -std=c++11 EventTree.cpp

Selector.o: EventTree.o Selector.cpp Selector.h
	g++ -c -I`root-config --incdir` -std=c++11 Selector.cpp

JECvariation.o: EventTree.o Selector.o JECvariation.cpp JECvariation.h
	g++ -c -I`root-config --incdir` -std=c++11 JECvariation.cpp

BTagCalibrationStandalone.o: BTagCalibrationStandalone.cpp BTagCalibrationStandalone.h
	g++ -c -I`root-config --incdir` -std=c++11 BTagCalibrationStandalone.cpp

EventPick.o: EventPick.cpp EventPick.h
	g++ -c -I `root-config --incdir` -std=c++11 EventPick.cpp 

Histogrammer.o: Histogrammer.cpp Histogrammer.h
	g++ -c -I `root-config --incdir` -std=c++11 Histogrammer.cpp

HistCollect.o: HistCollect.cpp HistCollect.h
	g++ -c -I`root-config --incdir` -std=c++11 HistCollect.cpp

clean:
	rm EventTree.o Selector.o JECvariation.o Histogrammer.o HistCollect.o EventPick.o
