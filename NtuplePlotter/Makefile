all: makeTemplates Selector.o EventTree.o EventPick.o makeSkim

makeTemplates: Histogrammer.o HistCollect.o Selector.o EventPick.o EventTree.o BTagCalibrationStandalone.o makeTemplates.cpp OverlapRemove.cpp PUReweight.cpp PUReweight.h   
	g++ -o makeTemplates `root-config --libs` -I`root-config --incdir` -std=c++11 EventTree.o Selector.o EventPick.o Histogrammer.o HistCollect.o   makeTemplates.cpp OverlapRemove.cpp PUReweight.cpp  BTagCalibrationStandalone.o  

makeSkim: Selector.o EventPick.o EventTree.o  makeSkim.cpp OverlapRemove.cpp
	g++ -o makeSkim `root-config --libs` -I`root-config --incdir` -std=c++11 EventTree.o EventPick.o Selector.o  OverlapRemove.cpp makeSkim.cpp 


EventTree.o: EventTree.cpp EventTree.h
	g++ -c -I`root-config --incdir` -std=c++11 EventTree.cpp

BTagCalibrationStandalone.o: BTagCalibrationStandalone.cpp BTagCalibrationStandalone.h
	g++ -c -I`root-config --incdir` -std=c++11 BTagCalibrationStandalone.cpp

Selector.o: EventTree.o Selector.cpp Selector.h
	g++ -c -I`root-config --incdir` -std=c++11 Selector.cpp


EventPick.o: EventPick.cpp EventPick.h
	g++ -c -I `root-config --incdir` -std=c++11 EventPick.cpp 

Histogrammer.o: Histogrammer.cpp Histogrammer.h
	g++ -c -I `root-config --incdir` -std=c++11 Histogrammer.cpp

HistCollect.o: HistCollect.cpp HistCollect.h
	g++ -c -I`root-config --incdir` -std=c++11 HistCollect.cpp

clean:
	rm EventTree.o Selector.o Histogrammer.o HistCollect.o EventPick.o
