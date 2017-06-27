#include<iostream>
#include<string>

#include"EventTree.h"
#include"Selector.h"
#include"EventPick.h"
#include"Histogrammer.h"
#include"HistCollect.h"
#include"PUReweight.h"

#include"TRandom3.h"
// temporary solution
//#include"JECvariation.h"
//#include"JECvariationNew.cpp"

//#include"/uscms_data/d2/dnoonan/LHAPDF/local/include/LHAPDF/LHAPDF.h"


//#include"LHAPDF/LHAPDF.h"

// using namespace LHAPDF;
// using namespace std;

int jecvar012_g = 1; // 0:down, 1:norm, 2:up
int jervar012_g = 1; // 0:down, 1:norm, 2:up
int mueff012_g = 1; // 0:down, 1:norm, 2:up
int eleeff012_g = 1;
int btagvar012_g = 1; // 0:down, 1:norm, 2:up
int phosmear012_g = 1; // 0:down, 1:norm, 2:up 
int musmear012_g = 1; // 0:down, 1:norm, 2: up

int toppt012_g = 1; // 0:down, 1:norm, 2: up

int top_sample_g = 0; // 0: no ttbar, 1: ttjets_1l, 2: ttjets_2l, 3: ttjets_had

int pdfweight_g = 1; //0:down, 1:norm, 2: up

//vector<LHAPDF::PDF*> pdfs;

double topPtWeight(EventTree* tree);
double getMuEff(EventTree* tree, EventPick* evt );
double getEleEff(EventTree* tree, EventPick* evt);
double getBtagSF(EventTree* tree, EventPick* evt);
void doMuSmearing(EventTree* tree);
void doPhoSmearing(EventTree* tree);
void doJER(EventTree* tree);
double JERcorrection(double eta);
bool overlapWHIZARD(EventTree* tree);
bool overlapMadGraph(EventTree* tree);
bool overlapISRFSR(EventTree* tree);
//double pdfWeight(EventTree* tree, int pdfNum);
double WjetsBRreweight(EventTree* tree);

int main(int ac, char** av){
	if(ac < 4){
		std::cout << "usage: ./makeTemplates sampleName outputDir inputFile[s]" << std::endl;
		return -1;
	}
	std::string PUfilename = "MyDataPileupHistogram.root";
	bool systematics = false;
	
	std::string inpFileName(av[3]);
	if( inpFileName.find("TTbar") != std::string::npos) top_sample_g = 1;
	//if( inpFileName.find("ttjets_2l") != std::string::npos) top_sample_g = 2;
	//if( inpFileName.find("ttjets_had") != std::string::npos) top_sample_g = 3;
	std::cout << "top_sample: " << top_sample_g << std::endl;
	
	std::string outDirName(av[2]);
	if( outDirName.find("JEC_up") != std::string::npos) {systematics=true; jecvar012_g = 2;}
	if( outDirName.find("JEC_down") != std::string::npos) {systematics=true; jecvar012_g = 0;}
	if( outDirName.find("JER_up") != std::string::npos) {systematics=true; jervar012_g = 2;}
	if( outDirName.find("JER_down") != std::string::npos) {systematics=true; jervar012_g = 0;}
	if( outDirName.find("EleEff_up") != std::string::npos) {systematics=true; eleeff012_g = 2;}
       	if( outDirName.find("EleEff_down") != std::string::npos) {systematics=true; eleeff012_g = 0;}
	if( outDirName.find("MuEff_up") != std::string::npos) {systematics=true; mueff012_g = 2;}
        if( outDirName.find("MuEff_down") != std::string::npos) {systematics=true; mueff012_g = 0;}
	if( outDirName.find("Btag_up") != std::string::npos) {systematics=true; btagvar012_g = 2;}
	if( outDirName.find("Btag_down") != std::string::npos) {systematics=true; btagvar012_g = 0;}
	if( outDirName.find("pho_up") != std::string::npos) {systematics=true; phosmear012_g = 2;}
	if( outDirName.find("pho_down") != std::string::npos) {systematics=true; phosmear012_g = 0;}
	if( outDirName.find("musmear_up") != std::string::npos) {systematics=true; musmear012_g = 2;}
	if( outDirName.find("musmear_down") != std::string::npos) {systematics=true; musmear012_g = 0;}
	if( outDirName.find("PU_up") != std::string::npos) {systematics=true; PUfilename = "MyDataPileupHistogram_up.root";}
	if( outDirName.find("PU_down") != std::string::npos) {systematics=true; PUfilename = "MyDataPileupHistogram_down.root";}
	if( outDirName.find("toppt_up") != std::string::npos) {systematics=true; toppt012_g = 2;}
	if( outDirName.find("toppt_down") != std::string::npos) {systematics=true; toppt012_g = 0;}	
	if( outDirName.find("PDF") != std::string::npos) {systematics=true; pdfweight_g=2;}

	int pdfNum = 0;
	if( pdfweight_g==2 ){
	  string tempNum = "";
	  for (int i = outDirName.find("PDF_")+4;i <outDirName.length(); i++){
	    if (outDirName[i]=='_' || outDirName[i]=='/'){i+outDirName.length();}
	    tempNum += outDirName[i];
	  }
	  pdfNum = atoi(tempNum.c_str());
	  std::cout << "PDF number: " << pdfNum << endl;
	}

	std::cout << "JEC: " << jecvar012_g << "  JER: " << jervar012_g << "  MuEff: " << mueff012_g << "  BtagVar: " << btagvar012_g << "  ";
	std::cout << "  PhoSmear: " << phosmear012_g << "  muSmear: " << musmear012_g << "  pileup: " << PUfilename << "  ";
	std::cout << "  topPt: " << toppt012_g << std::endl;
	// book HistCollect
	HistCollect* looseCollect = new HistCollect("PreSel",std::string("top_")+av[1]);
	//	looseCollect->fillEndcap = false;
	looseCollect->fillEndcap = true;
	HistCollect* looseCollectNoMET = new HistCollect("PreSelNoMET",std::string("top_")+av[1]);
	//	looseCollectNoMET->fillEndcap = false;	
	looseCollectNoMET->fillEndcap = true;	
	//HistCollect* fourjCollect = new HistCollect("1pho4j",std::string("top4j_")+av[1]);
	// HistCollect for tight Photon ID
	//HistCollect* tightCollect = new HistCollect("1photight",std::string("top_")+av[1]);
	
	// object selectors
	Selector* selectorLoose = new Selector();
	bool isQCD = false;

	// create event selectors here
	EventPick* evtPickLoose = new EventPick("LoosePhotonID");
	EventPick* evtPickLooseNoMET = new EventPick("LoosePhotonID");

	std::cout << av[2] << std::endl;
	if(std::string(av[2]).find("QCD") != std::string::npos){
	        std::cout << "IS QCD" << std::endl;
		isQCD = true;
		//selectorLoose->mu_MVA_range[0] = -1.0;
		//selectorLoose->mu_MVA_range[1] = -0.1;
		selectorLoose->mu_RelIso_range[0] = 0.25;
		selectorLoose->mu_RelIso_range[1] = 1.0;
		//selectorLoose->ele_Iso_MVA_invert = true;

		// no MC categories for data-driven QCD needed
		looseCollect->fillRS = false;
		looseCollect->fillFE = false;
		looseCollect->fillFJRB = false;

		looseCollectNoMET->fillRS = false;
		looseCollectNoMET->fillFE = false;
		looseCollectNoMET->fillFJRB = false;
		
		evtPickLoose->NlooseMuVeto_le = 99.;
		evtPickLooseNoMET->NlooseMuVeto_le = 99.;

	}
	//Selector* selectorTight = new Selector();
	// set up the parameters for object selectors here
	//selectorTight->pho_ID_ind = 2; // tight ID
	//std::cout << selectorLoose->mu_RelIso_range[0] << std::endl;
	//std::cout << selectorLoose->mu_RelIso_range[1] << std::endl;
	evtPickLoose->MET_cut = 20.0;
	evtPickLooseNoMET->MET_cut = -1.0;
	//evtPickLoose->veto_pho_jet_dR = 0.05;
	//evtPickLoose->Njet_ge = 4;
	//evtPickLoose->NBjet_ge = 2;
	// evtPickLoose->NlooseMuVeto_le = 0;
	// evtPickLoose->NlooseEleVeto_le = 0;

	if( outDirName.find("zeroB") != std::string::npos){
		evtPickLoose->NBjet_ge = 0;
		evtPickLooseNoMET->NBjet_ge = 0;
	}
	
	if( outDirName.find("twoMu") != std::string::npos){
		evtPickLoose->Nmu_eq = 2;
		evtPickLooseNoMET->Nmu_eq = 2;
	}

	if( outDirName.find("twoEle") != std::string::npos){
		evtPickLoose->Nele_eq = 2;
		evtPickLooseNoMET->Nele_eq = 2;
	}

	bool WHIZARD = false;
	if( std::string(av[1]).find("WHIZARD") != std::string::npos) WHIZARD = true;

	bool MGttgamma = false;
	if( std::string(av[1]).find("ttgamma") != std::string::npos) MGttgamma = true;
	if (MGttgamma) std::cout << "THIS IS TTGAMMA" << std::endl;

	bool doOverlapRemoval = false;
	bool doOverlapRemovalWZ = false;
	bool doInvertedOverlapRemoval = false;
	bool skipOverlap = false;
	if( std::string(av[1]).find("TTbar") != std::string::npos) doOverlapRemoval = true;
	if( std::string(av[1]).find("DYJets") != std::string::npos) doOverlapRemovalWZ = true;
	if( std::string(av[1]).find("W1jets") != std::string::npos) doOverlapRemovalWZ = true;
	if( std::string(av[1]).find("W2jets") != std::string::npos) doOverlapRemovalWZ = true;
	if( std::string(av[1]).find("W3jets") != std::string::npos) doOverlapRemovalWZ = true;
	if( std::string(av[1]).find("W4jets") != std::string::npos) doOverlapRemovalWZ = true;
	if( std::string(av[1]).find("Invert") != std::string::npos) doInvertedOverlapRemoval = true;
	if( std::string(av[1]).find("SkipOverlap") != std::string::npos) skipOverlap = true;

	if(skipOverlap){
	  doOverlapRemoval = false;
	  doOverlapRemovalWZ = false;
	  doInvertedOverlapRemoval = false;
	}
	if(doOverlapRemoval) std::cout << "########## Will apply overlap removal ###########" << std::endl;
	if(doOverlapRemovalWZ) std::cout << "########## Will apply WZ overlap removal ###########" << std::endl;
	if(doInvertedOverlapRemoval) std::cout << "########## Will invert overlap removal ###########" << std::endl;

	EventTree* tree = new EventTree(ac-3, av+3);
	double PUweight = 1.0;
	bool isMC = false;
	
	// initialize PU reweighting here
	PUReweight* PUweighter = new PUReweight(ac-3, av+3, PUfilename);
	
	tree->GetEntry(0);
	isMC = !(tree->isData_);
	//JECvariation* jecvar;
	//jecvar = new JECvariation("./Spring16_25nsV6_MC", isMC);

	// we don't need systematics variations for Data
	if(!isMC && systematics) {
		std::cout << "Systematics on Data is not applied" << std::endl;
		delete tree;
		return 0;
	}
	

	//	LHAPDF::PDFSet set("CT10nnlo");
	//	pdfs = set.mkPDFs();

	// initPDFSet(1, "cteq66.LHgrid");
//	TTree* newTree_pre = tree->chain->CloneTree(0);
	Long64_t nEntr = tree->GetEntries();
	for(Long64_t entry=0; entry<nEntr; entry++){
		if(entry%10000 == 0) std::cout << "processing entry " << entry << " out of " << nEntr << std::endl;
		//if(entry==5000000) break;
		tree->GetEntry(entry);
	//	std::cout<<"tree"<< tree <<std::endl;
	//	isMC = false;
		isMC = !(tree->isData_);
	//	std::cout<< "apply PU reweighting"<<std::endl;	
		//apply PU reweighting
		if(isMC) PUweight = PUweighter->getWeight(tree->nPUInfo_, tree->puBX_, tree->puTrue_);
	//	std::cout << "PUweight: " << PUweight <<std::endl;
		
		if(isMC && !isQCD){
			// JEC
		//	jecvar->applyJEC(tree, jecvar012_g); // 0:down, 1:norm, 2:up
			//JER smearing 
			doJER(tree);
		//	std::cout << "doing JER"<<std::endl;		
			// photon energy smearing
			doPhoSmearing(tree);
			// electron energy smearing
			//doEleSmearing(tree);
			doMuSmearing(tree);
		}
		 //do overlap removal here: overlapMadGraph(tree) or overlapWHIZARD(tree)
	        if( isMC && doOverlapRemoval){
		  if (!doInvertedOverlapRemoval){
		    if (overlapMadGraph(tree)){
		      continue;
		    }
		  }
		  if (doInvertedOverlapRemoval){
		    if (!overlapMadGraph(tree)){
		      continue;
		    }
		  }
	}
		    
		// if( isMC && doOverlapRemoval && overlapMadGraph(tree)){
		// 	//std::cout << "overlap!" << std::endl;
		// 	// overlapping part, not needed
		// 	continue;
		// }
         	 if( isMC && doOverlapRemovalWZ && overlapISRFSR(tree)){
	//		continue;
		}
		// this part cuts MadGraph ttgamma phase space to match WHIZARD, for comparison
		//if( isMC && MGttgamma && doOverlapRemoval && !overlapWHIZARD(tree)){
		//	// do not need MadGraph ttgamma events that are outside WHIZARD phase space
		//	continue;
		//}		

		selectorLoose->process_objects(tree);
		//selectorTight->process_objects(tree);
		evtPickLoose->process_event(tree, selectorLoose, PUweight);
		evtPickLooseNoMET->process_event(tree, selectorLoose, PUweight);
		//evtPickLoose4j->process_event(tree, selectorLoose, PUweight);
		//evtPickTight->process_event(tree, selectorTight, PUweight);
		double evtWeight = PUweight;
		if(isMC && !isQCD){
			//electron trigger efficiency reweighting
			evtWeight *= getMuEff(tree, evtPickLoose);
			// b-tag SF reweighting
			evtWeight *= getBtagSF(tree, evtPickLoose);
		}
		//double evtWeight = PUweight;
		if(isMC){
			//electron trigger efficiency reweighting
		//	evtWeight *= getEleEff(tree, evtPickLoose);
			// b-tag SF reweighting
			evtWeight *= getBtagSF(tree, evtPickLoose);
		}
		// top pt reweighting
		if(isMC){
			evtWeight *= topPtWeight(tree);
		}

		//if(isMC && pdfweight_g!=1){
		  
		   //double tempWeight = pdfWeight(tree,pdfNum);
		  // evtWeight *= tempWeight;
		//   // cout << tempWeight << endl;
		//}

	//	if(isMC && MGttgamma){
	//	 double tempWeight = WjetsBRreweight(tree);
	//	 evtWeight *= tempWeight;
	//	}

		// fill the histograms
		//std::cout << "fill, weight " << evtWeight << "  passPresel " << evtPickLoose->passPreSel << std::endl;
	//	if (evtPickLoose->passPreSel){
				//newTree_pre->Fill();
	//	double evtWeight = 1.0;
		//std::cout << "PUweight" << evtWeight <<std::endl;
	///	double PUweight = 1.0;		
	//	std::cout<<"Stage before filling histogarams"<<std::endl;
		looseCollectNoMET->fill_histograms(selectorLoose, evtPickLooseNoMET, tree, isMC, evtWeight);
	//	std::cout<<"Stage after NoMET"<<std::endl; 
		looseCollect->fill_histograms(selectorLoose, evtPickLoose, tree, isMC, evtWeight);
	//	}
	//	std::cout<<"Stage after filling histogarams"<<std::endl;

		//fourjCollect->fill_histograms(selectorLoose, evtPickLoose4j, tree, isMC, evtWweight);
	}
	
	looseCollect->write_histograms(evtPickLoose, isMC, av[2]);
	looseCollectNoMET->write_histograms(evtPickLooseNoMET, isMC, av[2]);

	//fourjCollect->write_histograms(evtPickLoose4j, isMC, av[2]);

//	std::cout << "Average PU weight " << PUweighter->getAvgWeight() << std::endl;
	evtPickLoose->print_cutflow();
	
	delete tree;
	return 0;
}

//newTree_pre->Write();
double muTrigSF(double pt, double eta){
	static double trigEffSF_PtEta[3][3] = { {0.984, 0.967, 0.991}, {0.999, 0.983, 1.018}, {0.999, 0.988, 0.977} };
	static double trigEffSFerr_PtEta[3][3] = { {0.002, 0.002, 0.007}, {0.003, 0.002, 0.012}, {0.002, 0.003, 0.015} };
	
	int etaRegion = 0;
	if( eta > 0.80) etaRegion++;
	if( eta > 1.48) etaRegion++;

	int ptRegion = 0;
	if( pt > 40 ) ptRegion++;
	if( pt > 50 ) ptRegion++;

	if(mueff012_g == 1) return trigEffSF_PtEta[ptRegion][etaRegion];
	if(mueff012_g == 0) return trigEffSF_PtEta[ptRegion][etaRegion] - trigEffSFerr_PtEta[ptRegion][etaRegion];
	if(mueff012_g == 2) return trigEffSF_PtEta[ptRegion][etaRegion] + trigEffSFerr_PtEta[ptRegion][etaRegion];
	return 1.0;	
}

double muIDSF(double pt, double eta){
	static double idEffSF_PtEta[3][3] = { {0.950, 0.957, 0.922}, {0.966, 0.961, 0.941}, {0.961, 0.963, 0.971} };
	static double idEffSFerr_PtEta[3][3] = { {0.003, 0.002, 0.004}, {0.001, 0.002, 0.007}, {0.002, 0.003, 0.0} };

	int etaRegion = 0;
	if( eta > 0.80) etaRegion++;
	if( eta > 1.48) etaRegion++;

	int ptRegion = 0;
	if( pt > 40 ) ptRegion++;
	if( pt > 50 ) ptRegion++;

	if(mueff012_g == 1) return idEffSF_PtEta[ptRegion][etaRegion];
	if(mueff012_g == 0) return idEffSF_PtEta[ptRegion][etaRegion] - idEffSFerr_PtEta[ptRegion][etaRegion];
	if(mueff012_g == 2) return idEffSF_PtEta[ptRegion][etaRegion] + idEffSFerr_PtEta[ptRegion][etaRegion];	
	return 1.0;
}

// AN2012_438_v10 page9
double getMuEff(EventTree* tree, EventPick* evt){
	if( evt->Muons.size() < 1 ) return 1.0; // no electrons, no weight
	int muInd = evt->Muons[0];
	double pt = tree->muPt_->at(muInd);
	double eta = TMath::Abs(tree->muEta_->at(muInd));
	
	double trigSF = muTrigSF(pt, eta);
	double idSF = muIDSF(pt, eta);
	if(evt->Muons.size() == 1) return trigSF*idSF;
	if(evt->Muons.size() == 2){
		int muInd2 = evt->Muons[1];	
		double pt2 = tree->muPt_->at(muInd2);
		double eta2 = TMath::Abs(tree->muEta_->at(muInd2));
		double trigSF2 = muTrigSF(pt2, eta2);
		double idSF2 = muIDSF(pt2, eta2);
		return ( 1.0 - (1.0-trigSF)*(1.0-trigSF2) )*idSF*idSF2;
	}
	return 1.0;
}



double eleTrigSF(double pt, double eta){
	static double trigEffSF_PtEta[3][3] = { {0.984, 0.967, 0.991}, {0.999, 0.983, 1.018}, {0.999, 0.988, 0.977} };
	static double trigEffSFerr_PtEta[3][3] = { {0.002, 0.002, 0.007}, {0.003, 0.002, 0.012}, {0.002, 0.003, 0.015} };
	
	int etaRegion = 0;
	if( eta > 0.80) etaRegion++;
	if( eta > 1.48) etaRegion++;

	int ptRegion = 0;
	if( pt > 40 ) ptRegion++;
	if( pt > 50 ) ptRegion++;

	if(eleeff012_g == 1) return trigEffSF_PtEta[ptRegion][etaRegion];
	if(eleeff012_g == 0) return trigEffSF_PtEta[ptRegion][etaRegion] - trigEffSFerr_PtEta[ptRegion][etaRegion];
	if(eleeff012_g == 2) return trigEffSF_PtEta[ptRegion][etaRegion] + trigEffSFerr_PtEta[ptRegion][etaRegion];
	return 1.0;	
}

double eleIDSF(double pt, double eta){
	static double idEffSF_PtEta[3][3] = { {0.950, 0.957, 0.922}, {0.966, 0.961, 0.941}, {0.961, 0.963, 0.971} };
	static double idEffSFerr_PtEta[3][3] = { {0.003, 0.002, 0.004}, {0.001, 0.002, 0.007}, {0.002, 0.003, 0.0} };

	int etaRegion = 0;
	if( eta > 0.80) etaRegion++;
	if( eta > 1.48) etaRegion++;

	int ptRegion = 0;
	if( pt > 40 ) ptRegion++;
	if( pt > 50 ) ptRegion++;

	if(eleeff012_g == 1) return idEffSF_PtEta[ptRegion][etaRegion];
	if(eleeff012_g == 0) return idEffSF_PtEta[ptRegion][etaRegion] - idEffSFerr_PtEta[ptRegion][etaRegion];
	if(eleeff012_g == 2) return idEffSF_PtEta[ptRegion][etaRegion] + idEffSFerr_PtEta[ptRegion][etaRegion];	
	return 1.0;
}

// AN2012_438_v10 page9
double getEleEff(EventTree* tree, EventPick* evt){
	if( evt->Electrons.size() < 1 ) return 1.0; // no electrons, no weight
	int eleInd = evt->Electrons[0];
	double pt = tree->elePt_->at(eleInd);
	double eta = TMath::Abs(tree->eleSCEta_->at(eleInd));
	
	double trigSF = eleTrigSF(pt, eta);
	double idSF = eleIDSF(pt, eta);
	if(evt->Electrons.size() == 1) return trigSF*idSF;
	if(evt->Electrons.size() == 2){
		int eleInd2 = evt->Electrons[1];	
		double pt2 = tree->elePt_->at(eleInd2);
		double eta2 = TMath::Abs(tree->eleSCEta_->at(eleInd2));
		double trigSF2 = eleTrigSF(pt2, eta2);
		double idSF2 = eleIDSF(pt2, eta2);
		return ( 1.0 - (1.0-trigSF)*(1.0-trigSF2) )*idSF*idSF2;
	}
	return 1.0;
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
double SFtop(double pt){
	if(top_sample_g==1) return exp(0.159 - 0.00141*pt);
	if(top_sample_g==2) return exp(0.148 - 0.00129*pt);
	return 1.0;
}

double topPtWeight(EventTree* tree){
	double toppt=0.0;
	double antitoppt=0.0;
	double weight = 1.0;
	for(int mcInd=0; mcInd<tree->nMC_; ++mcInd){
		if(tree->mcPID->at(mcInd)==6) toppt = tree->mcPt->at(mcInd);
		if(tree->mcPID->at(mcInd)==-6) antitoppt = tree->mcPt->at(mcInd);
	}
	if(toppt > 0.001 && antitoppt > 0.001)
		weight = sqrt( SFtop(toppt) * SFtop(antitoppt) );
	
	// if(toppt012_g == 1) return weight;
	// if(toppt012_g == 0) return 1.0;
	// if(toppt012_g == 2) return weight*weight;

	//This has been changed, the new prescription is to not use the top pt reweighting, and the syst is using it
	if(toppt012_g == 1) return 1.0;
	if(toppt012_g == 0) return weight;
	if(toppt012_g == 2) return weight;


	// should not get here
	return 1.0;
}

void doMuSmearing(EventTree* tree){
	static TRandom3 rand;
	if(musmear012_g == 1) return;
	for(int muInd = 0; muInd < tree->nMu_; ++muInd){
		if(tree->muPt_->at(muInd) < 15) continue;
		//std::cout << "electron Pt before " << tree->elePt_->at(eleInd) << "   "; 
		double factor = 1.0;
		if(musmear012_g == 0) factor = 0.99;
		if(musmear012_g == 2) factor = 1.01;
		//std::cout << "factor " << factor << "  ";
		tree->muPt_->at(muInd) *= factor;
		//std::cout << "electron Pt after " << tree->elePt_->at(eleInd) << std::endl;
	}
}

void doPhoSmearing(EventTree* tree){
	if(phosmear012_g == 1) return;
	for(int phoInd = 0; phoInd < tree->nPho_; ++phoInd){
		if(tree->phoEt_->at(phoInd) < 15) continue;
		//std::cout << "photon Et before " << tree->phoEt_->at(phoInd) << "   "; 
		double factor = 1.0;
		if(phosmear012_g == 0) factor = 0.99;
		if(phosmear012_g == 2) factor = 1.01;
		//std::cout << "factor " << factor << "  ";
		tree->phoEt_->at(phoInd) *= factor;
	}
}

// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
void doJER(EventTree* tree){
	// correct MET when jets are smeared
	TLorentzVector tMET;
	tMET.SetPtEtaPhiM(tree->pfMET_,0.0,tree->pfMETPhi_,0.0);
	//std::cout << "before correction MET " << tree->pfMET_ << " " << tree->pfMETPhi_ << "    ";
	// scale jets
	for(int jetInd = 0; jetInd < tree->nJet_ ; ++jetInd){
	        // This is changed, previously was 20 GeV (was this from 2011 recommendation?) 
		if(tree->jetPt_->at(jetInd) < 10) continue;   
		//if(tree->jetGenJetIndex_->at(jetInd)>0){
		//	TLorentzVector tjet;
		//	tjet.SetPtEtaPhiM(tree->jetPt_->at(jetInd), tree->jetEta_->at(jetInd), tree->jetPhi_->at(jetInd), 0.0);
		//	tMET+=tjet;
		//	double oldPt = tree->jetPt_->at(jetInd);
		//	double genPt = tree->jetGenJetPt_->at(jetInd);
		//	double eta = tree->jetEta_->at(jetInd);
		//	tree->jetPt_->at(jetInd) = std::max(0.0, genPt + JERcorrection(eta)*(oldPt-genPt));
		//	tjet.SetPtEtaPhiM(tree->jetPt_->at(jetInd), tree->jetEta_->at(jetInd), tree->jetPhi_->at(jetInd), 0.0);			
		//	tMET-=tjet;
			//std::cout << "old " << oldPt << "  new " << tree->jetPt_->at(jetInd) << std::endl;
		//}
	}
	// save updated MET values
	tree->pfMET_ = tMET.Pt();
	tree->pfMETPhi_ = tMET.Phi();
	//std::cout << "after corrections " << tree->pfMET_ << " " << tree->pfMETPhi_ << std::endl;
}


double JERcorrection(double JetEta){
	double eta = TMath::Abs(JetEta);
	// static const double corr[5] = {1.052, 1.057, 1.096, 1.134, 1.288};
	// static const double corrDown[5] = {0.990, 1.001, 1.032, 1.042, 1.089};
	// static const double corrUp[5] = {1.115, 1.114, 1.161, 1.228, 1.488};
	static const double corr[7] = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
	static const double corrDown[7] = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};
	static const double corrUp[7] = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};

	int region = 0;
	if( eta >= 0.5 ) region++;
	if( eta >= 1.1 ) region++;
	if( eta >= 1.7 ) region++;
	if( eta >= 2.3 ) region++;
	if( eta >= 2.8 ) region++;
	if( eta >= 3.2 ) region++;
	if(jervar012_g == 0) return corrDown[region];
	if(jervar012_g == 1) return corr[region];
	if(jervar012_g == 2) return corrUp[region];
	
	// should not get here
	return 1.0;
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods
// https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_NOttbar_payload_EPS13.txt
// weight for >=1 btag :  1 - prod(1-SFi) over all b-tagged jets

double lfJetCSVM(double x, double jeteta){
	if(jeteta < 0.8) return ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
	if(jeteta < 1.6) return ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
	return ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x))); 
}

double lfJetCSVMmin(double x, double jeteta){
	if(jeteta < 0.8) return ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
	if(jeteta < 1.6) return ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
	return ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
}

double lfJetCSVMmax(double x, double jeteta){
	if(jeteta < 0.8) return ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)));
	if(jeteta < 1.6) return ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)));
	return ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)));
}

double lfJetSF(double jetpt, double jeteta){
	bool pthigh = false;
	if(jeteta < 1.6 && jetpt > 1000.0) {jetpt = 1000.0; pthigh = true;}
	if(jeteta >= 1.6 && jetpt > 650.0) {jetpt = 650.0; pthigh = true;}
	
	if(!pthigh){
		if( btagvar012_g == 1 ) return lfJetCSVM(jetpt, jeteta);
		if( btagvar012_g == 0 ) return lfJetCSVMmin(jetpt, jeteta);
		if( btagvar012_g == 2 ) return lfJetCSVMmax(jetpt, jeteta);
	} else{ // in case jetpt is too high:
		if( btagvar012_g == 1 ) return lfJetCSVM(jetpt, jeteta);
		if( btagvar012_g == 0 ) return 2.0*lfJetCSVMmin(jetpt, jeteta) - lfJetCSVM(jetpt, jeteta);
		if( btagvar012_g == 2 ) return 2.0*lfJetCSVMmax(jetpt, jeteta) - lfJetCSVM(jetpt, jeteta);
	}
	std::cout << "should not be here!" << std::endl;
	return 1.0;
}

double bSFerr(double jetpt){
	//Tagger: CSVM within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
	static double ptmax[16] = {30, 40, 50, 60, 
				70, 80, 100, 120, 
				160, 210, 260, 320, 
				400, 500, 600, 800};
	static double SFb_error[17] = {0.0415694, 0.023429, 0.0261074, 0.0239251, 
				0.0232416, 0.0197251, 0.0217319, 0.0198108, 
				0.0193, 0.0276144, 0.020583, 0.026915, 
				0.0312739, 0.0415054, 0.074056, 0.0598311, 
				0.0598311*2};

	int ptInd = 0;
	if( btagvar012_g == 1) return 0.0;

	for(int ipt=0; ipt<16; ipt++) 
		if(jetpt > ptmax[ipt]) ptInd++;
	
	if( btagvar012_g == 0 ) return -1.0*SFb_error[ptInd];
	if( btagvar012_g == 2 ) return  SFb_error[ptInd];

	return 0.0;
}

double bJetSF(double jetpt){
	return (0.939158+(0.000158694*jetpt))+(-2.53962e-07*(jetpt*jetpt)) + bSFerr(jetpt);
}

double cjetSF(double jetpt){
	return (0.939158+(0.000158694*jetpt))+(-2.53962e-07*(jetpt*jetpt)) + 2.0*bSFerr(jetpt);
}

double getBtagSF(EventTree* tree, EventPick* evt){
	
	double prod = 1.0;
	double jetpt;
	double jeteta;
	int jetflavor;
	double SFb;

	if(evt->bJets.size() == 0) return 1.0;

	for(std::vector<int>::const_iterator bjetInd = evt->bJets.begin(); bjetInd != evt->bJets.end(); bjetInd++){
		jetpt = tree->jetPt_->at(*bjetInd);
		jeteta = fabs(tree->jetEta_->at(*bjetInd));
		jetflavor = abs(tree->jetPartonID_->at(*bjetInd));
		
		if(jetflavor == 5) SFb = bJetSF(jetpt);
		else if( jetflavor == 4) SFb = cjetSF(jetpt);
		else SFb = lfJetSF(jetpt, jeteta);
	
		prod *= 1.0 - SFb;
	}
	return 1.0 - prod;
}

// double pdfWeight(EventTree* tree, int pdfNum){
//   if( pdfweight_g==1){ return 1.;  } 
//   if (pdfNum >= pdfs.size() || pdfNum < 1){return 1.;}

//   double id1 = tree->pdf_[0];
//   double id2 = tree->pdf_[1];
//   double x1  = tree->pdf_[2];
//   double x2  = tree->pdf_[3];
//   double q   = tree->pdf_[6];

//   double pdfweightMin = 999.;
//   double pdfweightMax = 0.;

//   double xpdf1 = pdfs[0]->xfxQ(id1, x1, q);
//   double xpdf2 = pdfs[0]->xfxQ(id2, x2, q);

//   // LHAPDF::usePDFMember(1,0);
//   // double xpdf1 = LHAPDF::xfx(1, x1, q, id1);
//   // double xpdf2 = LHAPDF::xfx(1, x2, q, id2);
//   double w0 = xpdf1 * xpdf2;

//   double xpdf1_new = pdfs[pdfNum]->xfxQ(id1, x1, q);
//   double xpdf2_new = pdfs[pdfNum]->xfxQ(id2, x2, q);
//   // LHAPDF::usePDFMember(1,0);
//   // double xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
//   // double xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
//   double weight = xpdf1_new * xpdf2_new / w0;
//   return weight;

//   // for(int i=1; i < pdfs.size(); ++i){
//   //   double xpdf1_new = pdfs[i]->xfxQ(id1, x1, q);
//   //   double xpdf2_new = pdfs[i]->xfxQ(id2, x2, q);
//   //   // LHAPDF::usePDFMember(1,0);
//   //   // double xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
//   //   // double xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
//   //   double weight = xpdf1_new * xpdf2_new / w0;
//   //   if(weight > pdfweightMax) pdfweightMax=weight;
//   //   if(weight < pdfweightMin) pdfweightMin=weight;
//   // }

//   // if( pdfweight_g==0){ return pdfweightMin;}
//   // if( pdfweight_g==2){ return pdfweightMax;}
// }


double WjetsBRreweight(EventTree* tree){

  int countLeps = 0;
  //Need to try to avoid double counting of tau's (check if momPt is the same before counting, since status flag isn't there)
  double tauMomPt1 = -1.;
  double tauMomPt2 = -1.;
  //  std::cout << "-------------------" << std::endl;
  for(int mcInd=0; mcInd<tree->nMC_; ++mcInd){
    if( (TMath::Abs(tree->mcPID->at(mcInd)) == 11 ||
	 TMath::Abs(tree->mcPID->at(mcInd)) == 13 ||
	 TMath::Abs(tree->mcPID->at(mcInd)) == 15 ) &&
	TMath::Abs(tree->mcMomPID->at(mcInd)) == 24 &&
	TMath::Abs(tree->mcGMomPID->at(mcInd)) == 6) {
      if (TMath::Abs(tree->mcPID->at(mcInd)) == 15){
	if (tauMomPt1==-1.){
	  tauMomPt1 = tree->mcMomPt->at(mcInd);
	  countLeps += 1;
	}
	else if (tauMomPt2==-1.){
	  if (tree->mcMomPt->at(mcInd)!=tauMomPt1) {
	    tauMomPt2 = tree->mcMomPt->at(mcInd);
	    countLeps += 1;
	  }
	}
	else{
	  if (tree->mcMomPt->at(mcInd)!=tauMomPt1 && tree->mcMomPt->at(mcInd)!=tauMomPt2) {
	    countLeps += 1;
	  }
	}
      }
      else {
	countLeps += 1;
      }
      //      std::cout << tree->mcPID->at(mcInd) << "\t" << tree->mcMomPID->at(mcInd) << "\t" << tree->mcGMomPID->at(mcInd) << "\t" << tree->mcMomPt->at(mcInd) << "\t" << tree->mcIndex->at(mcInd) << std::endl;
    }
  }
  double reweight=1.;
  if (countLeps==0){reweight = .6741*.6741*9./4.;}
  else if(countLeps==1){reweight = .6741*.3259*2*9./4.;}
  else if(countLeps==2){reweight = .3259*.3259*9.;}
  else {
    std::cout << "MORE THAN TWO LEPTONS???????" << std::endl;
    std::cout << countLeps << std::endl;
  }

  //  std::cout << reweight << std::endl;
  return reweight;

}

