#include<iostream>
#include<string>
//#include"JECvariation.h"
#include"EventTree.h"
#include"Selector.h"
#include"EventPick.h"
#include<TFile.h>
#include<TTree.h>
#include<TDirectory.h>
#include<TObject.h>
#include<TH1F.h>
#include<TCanvas.h>

int main(int ac, char** av){
	if(ac < 3){
		std::cout << "usage: ./makeSkim outputFileName inputFile[s]" << std::endl;
		return -1;
	}
	// input: dealing with TTree first
	bool isMC = true;
	EventTree* tree = new EventTree(ac-2, av+2);
	Selector* selector = new Selector();
	selector->jet_Pt_cut = 30;
	EventPick* evtPick = new EventPick("nominal");
        evtPick->MET_cut = -1.0;	
//	JECvariation* JEC = new JECvariation("./Spring16_25nsV6_MC/Spring16_25nsV6", isMC);
	std::string outDirName(av[1]);
	// antiselection for QCD fit
	if( outDirName.find("QCD") != std::string::npos){
		std::cout << "muon antiselection is on" << std::endl;
		selector->mu_RelIso_range[0] = 0.25; 
		selector->mu_RelIso_range[1] = 1.;
		selector->mu_Iso_invert = true;
	}
	

	evtPick->NBjet_ge = 1;

	TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",1000,500);
        c1->SetFillColor(42);
        c1->SetGrid();

        TCanvas *c2 = new TCanvas("c2","A Simple Graph Example",1000,500);
        c2->SetFillColor(42);
        c2->SetGrid();


        TCanvas *c3 = new TCanvas("c3","A Simple Graph Example",1000,500);
        c3->SetFillColor(42);
        c3->SetGrid();

        TH1F *h1 = new TH1F("leading jet pt", "jet_pt", 20, 30, 600);
        TH1F *h2 = new TH1F("second jet pt", "jet_pt", 20, 30, 600);
        TH1F *h3 = new TH1F("third jet pt", "jet_pt", 20, 30, 600);
	
	
	//TFile *theFile = TFile::Open("root://cmsxrootd.fnal.gov//store/user/troy2012/rootFile.root");
	TFile* outFile = TFile::Open( av[1] ,"RECREATE" );
	TDirectory* ggDir = outFile->mkdir("ggNtuplizer","ggNtuplizer");
	ggDir->cd();
	TTree* newTree = tree->chain->CloneTree(0);
	
	Long64_t nEntr = tree->GetEntries();
	for(Long64_t entry= 0; entry < nEntr; entry++){
//	for(Long64_t entry= 0; entry < 210000; entry++){ 	
		if(entry%1000 == 0) {
			std::cout << "processing entry " << entry << " out of " << nEntr << std::endl;
		}
		tree->GetEntry(entry);
		//int check = tree->GetEntry(entry);
		//std::cout <<"tree->entry : " << check <<std::endl;
		//if (check == 0) continue;
		selector->process_objects(tree);
		
		//JEC->applyJEC(tree,1);
		evtPick->process_event(tree,selector);
		// make selection here
		if( evtPick->passSkim ){
			newTree->Fill();
		
                }
		
	}

	newTree->Write();
	evtPick->print_cutflow();
	std::map<std::string, TH1F*> histMap;
	// copy histograms
	for(int fileInd = 2; fileInd < ac; ++fileInd){
		TFile* tempFile = TFile::Open(av[fileInd], "READ");
		TIter next(((TDirectory*)tempFile->Get("ggNtuplizer"))->GetListOfKeys());
		TObject* obj;
		while ((obj = next())){
			std::string objName(obj->GetName());
			if( objName != "EventTree"){
				TH1F* hist = (TH1F*)tempFile->Get(("ggNtuplizer/"+objName).c_str());
				if( histMap.find(objName) != histMap.end() ){
					histMap[objName]->Add(hist);
				}
				else {
					hist->SetDirectory(0);
					histMap[objName] = hist;
				}
			}
		}
		tempFile->Close();
	}
	
	ggDir->cd();
	for(std::map<std::string, TH1F*>::iterator it = histMap.begin(); it!= histMap.end(); ++it){
		it->second->SetDirectory(ggDir);
		it->second->Write();
	}
	outFile->Close();
	
	return 0;
}
