#include"PUReweight.h"
#include<iostream>
#include<TCanvas.h> 
PUReweight::PUReweight(int nFiles, char** fileNames, std::string PUfilename){
	PUweightSum = 0.0;
	events = 0;
	TFile* pileupFile = new TFile(PUfilename.c_str(),"READ");
	PUweightHist = (TH1D*)pileupFile->Get("pileup");
	PUweightHist->SetDirectory(0);
	pileupFile->Close();
	TH1D* PUbackup;
	if(PUweightHist->GetNbinsX() != 200){
		std::cout << "Wrong number of bins in the pileup histogram" << std::endl;
		PUbackup = new TH1D("pileup_new","pileup_new",200,0,200);
		for(int ibin=1; ibin <= PUweightHist->GetNbinsX(); ibin++){
			PUbackup->SetBinContent(ibin, PUweightHist->GetBinContent(ibin));
			// assuming the same scale
		}
		PUweightHist = PUbackup;
	}
	double PUweightInt = PUweightHist->Integral();
	TH1F* mcPU = NULL;
	for(int nmcfile = 0; nmcfile<nFiles; nmcfile++){
		std::cout << "reading file " << std::string(fileNames[nmcfile]) << std::endl;
		TFile* mcFile = TFile::Open(fileNames[nmcfile],"READ");
		if(!(mcFile->Get("ggNtuplizer/hPUTrue"))) {
			std::cout << "no hPU histogram here!" << std::endl;
			delete PUweightHist;
			PUweightHist = NULL;
			return;
		}
		if( mcPU==NULL) mcPU = (TH1F*)mcFile->Get("ggNtuplizer/hPUTrue");
		else mcPU->Add((TH1F*)mcFile->Get("ggNtuplizer/hPUTrue"));
		mcPU->SetDirectory(0);
		mcFile->Close();
	}
	if (mcPU->GetNbinsX()!=PUweightHist->GetNbinsX()){
		if (mcPU->GetNbinsX()>PUweightHist->GetNbinsX()){
			mcPU->Rebin(mcPU->GetNbinsX()/PUweightHist->GetNbinsX());
		}
	}
	mcPU->Scale(1.0/mcPU->Integral());
	PUweightHist->Divide(mcPU);
	PUweightHist->Scale(1.0/PUweightInt);
	PUweightHist->Draw();
	delete mcPU;
}

PUReweight::~PUReweight(){
	delete PUweightHist;
}

double PUReweight::getWeight(int nPUInfo, std::vector<int> *puBX, std::vector<float> *puTrue){
	double PUweight=0.0;
	if(!PUweightHist) {
	  std::cout << "you are calling getWeight by mistake" << std::endl; 
	  return 1.0;}
	
	for(int puInd=0; puInd<nPUInfo; ++puInd){
		if( puBX->at(puInd) == 0 ){
			PUweight = PUweightHist->GetBinContent(PUweightHist->GetXaxis()->FindBin(puTrue->at(puInd)));
			break;
		}
	}
	events++;
	PUweightSum+=PUweight;
	return PUweight;
}

double PUReweight::getAvgWeight(){
	if(events!=0) return PUweightSum/events;
	else return -1.0;
}
