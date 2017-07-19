#include<iostream>
#include<string>
#include"EventTree.h"
#include"Selector.h"
#include"EventPick.h"
#include"Histogrammer.h"
#include"HistCollect.h"
#include"PUReweight.h"

#include"TRandom3.h"

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

#include "BTagCalibrationStandalone.h"
double topPtWeight(EventTree* tree); 
double getMuEff(EventTree* tree, EventPick* evt, int systlevel );
double getEleEff(EventTree* tree, EventPick* evt);
double getBtagSF(EventTree* tree, EventPick* evt,string sysType, BTagCalibrationReader reader);
void doMuSmearing(EventTree* tree);
void doPhoSmearing(EventTree* tree);
void doJER(EventTree* tree);
double JERcorrection(double eta);
bool overlapWHIZARD(EventTree* tree);
bool overlapMadGraph(EventTree* tree);
bool overlapISRFSR(EventTree* tree);
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

	BTagCalibration calib("csvv2", "CSVv2_Moriond17_B_H.csv");

	BTagCalibrationReader reader(BTagEntry::OP_MEDIUM,"central", {"up", "down"});	

	reader.load(calib, BTagEntry::FLAV_B, "comb");               // measurement type

	reader.load(calib, BTagEntry::FLAV_C, "comb");               // measurement type

	reader.load(calib, BTagEntry::FLAV_UDSG, "incl");               // measurement type

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
	

	int count_overlapTTbar=0;
        int count_overlapVJets=0;
	Long64_t nEntr = tree->GetEntries();
	for(Long64_t entry=0; entry<nEntr; entry++){
		if(entry%10000 == 0) std::cout << "processing entry " << entry << " out of " << nEntr << std::endl;
		tree->GetEntry(entry);
	//	std::cout << "entry is " << entry << std::endl;
	//	std::cout << "tree is " << tree <<std::endl;
		isMC = !(tree->isData_);
		if(isMC) PUweight = PUweighter->getWeight(tree->nPUInfo_, tree->puBX_, tree->puTrue_);
		
		if(isMC && !isQCD){
			doJER(tree);
			doPhoSmearing(tree);
			doMuSmearing(tree);
		}
		 //do overlap removal here: overlapMadGraph(tree) or overlapWHIZARD(tree)
	        if( isMC && doOverlapRemoval){
		
		    if (overlapMadGraph(tree)){
		      count_overlapTTbar++;
		      continue;
		    }
		  
		}
		    
         	 if( isMC && doOverlapRemovalWZ){ 
		     if(overlapISRFSR(tree)){
			count_overlapVJets++;
			continue;
			}
		 }

		selectorLoose->process_objects(tree);
		evtPickLoose->process_event(tree, selectorLoose, PUweight);
		evtPickLooseNoMET->process_event(tree, selectorLoose, PUweight);
		double evtWeight = PUweight;
	//	std::cout << "getting BTag and MuSF "<< std::endl;
		if(isMC && !isQCD){
			//Muon trigger efficiency reweighting
			evtWeight *= getMuEff(tree, evtPickLoose, 1);
			// b-tag SF reweighting
			evtWeight *= getBtagSF(tree, evtPickLoose,"central", reader);
	//		std::cout << "After BTag:"<<  evtWeight <<std::endl;
		}
	//	std::cout << "after BTag"<< std::endl;
		// top pt reweighting
//		if(isMC){
//			evtWeight *= topPtWeight(tree);
//		}


//		if(isMC && MGttgamma){
//		 double tempWeight = WjetsBRreweight(tree);
//		 evtWeight *= tempWeight;
/////		}
	//	std::cout << "start filling histo"<< std::endl;
		looseCollectNoMET->fill_histograms(selectorLoose, evtPickLooseNoMET, tree, isMC, evtWeight);
		looseCollect->fill_histograms(selectorLoose, evtPickLoose, tree, isMC, evtWeight);
	}
	std::cout << "Total number of events removed from TTbar:"<< count_overlapTTbar <<std::endl;
        std::cout << "Total number of events removed from W/ZJets:"<< count_overlapVJets <<std::endl;
	
	looseCollect->write_histograms(evtPickLoose, isMC, av[2]);
	looseCollectNoMET->write_histograms(evtPickLooseNoMET, isMC, av[2]);
	//std::cout << "after filling histo"<< std::endl;
	evtPickLoose->print_cutflow();
	
	delete tree;
	return 0;
}

//https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
double muTrigSF[7][4][3]={{{0.979945,0.980956,0.981966}, {0.954675,0.956924,0.959174}, {0.980779,0.982631,0.984484}, {0.903027,0.906908,0.910789}, },
  			      	{ {0.984395,0.984724,0.985052}, {0.965922,0.966496,0.967069}, {0.994957,0.995577,0.996197}, {0.943534,0.944924,0.946313}, },
			     	{ {0.985344,0.985599,0.985853}, {0.968089,0.968468,0.968846}, {0.999035,0.999465,0.999894}, {0.956917,0.958023,0.959128}, },
			     	{ {0.985190,0.985745,0.986299}, {0.968046,0.968827,0.969609}, {0.998521,0.999418,1.000314}, {0.958654,0.961000,0.963346}, },
			     	{ {0.983833,0.984741,0.985649}, {0.964031,0.965445,0.966858}, {0.997885,0.999318,1.000752}, {0.949569,0.953347,0.957124}, },
			     	{ {0.971605,0.976515,0.981424}, {0.941747,0.947421,0.953094}, {0.997657,1.005482,1.013308}, {0.946883,0.977561,1.008239}, },
			     	{ {0.975026,0.984242,0.993459}, {0.934809,0.952092,0.969374}, {0.959526,0.984215,1.008903}, {0.751047,0.919715,1.088383}, } };


double muIdIsoSF[6][4][3] = {{{0.982175,0.984434,0.986693}, {0.991846,0.995412,0.998978}, {0.990022,0.991702,0.993382}, {0.981885,0.984787,0.987688}, },
					{ {0.992057,0.993262,0.994468}, {0.998298,1.000365,1.002432}, {0.995093,0.996092,0.997091}, {0.990776,0.992530,0.994284}, },
					{ {0.993362,0.993752,0.994143}, {0.998745,0.999447,1.000149}, {0.997677,0.998057,0.998437}, {0.995883,0.996573,0.997263}, },
					{ {0.995116,0.995289,0.995461}, {0.997599,0.997752,0.997904}, {0.997983,0.998067,0.998150}, {0.998147,0.998495,0.998844}, },
					{ {0.996429,0.996797,0.997165}, {0.998515,0.999121,0.999726}, {0.998046,0.998386,0.998725}, {0.997953,0.998730,0.999506}, },
					{ {0.998338,0.998806,0.999274}, {0.998320,0.999103,0.999885}, {0.998777,0.999235,0.999693}, {1.000398,1.001523,1.002648}, } };



double muTrackingSF[12][3] = { {0.996924,0.996996,0.997069},
								{0.997629,0.997712,0.997794},
								{0.998007,0.998078,0.998149},
								{0.997729,0.997804,0.997878},
								{0.997863,0.997971,0.998077},
								{0.996962,0.997148,0.997334},
								{0.996047,0.996227,0.996409},
								{0.995308,0.995479,0.995649},
								{0.995606,0.995781,0.995958},
								{0.993657,0.993892,0.994127},
								{0.992617,0.992943,0.993273},
								{0.986461,0.987313,0.988173}};




// AN2012_438_v10 page9
double getMuEff(EventTree* tree, EventPick* evt, int systLevel){
	if( evt->Muons.size() < 1 ) return 1.0; // no electrons, no weight
	int muInd = evt->Muons[0];
	double pt = tree->muPt_->at(muInd);
	double eta = TMath::Abs(tree->muEta_->at(muInd));
	
	int muTrackEtaRegion = int(eta/0.2);
	int muEtaRegion = -1;
	if (eta < 0.9) {muEtaRegion = 0;}
	else if (eta < 1.2) {muEtaRegion = 1;}
	else if (eta < 2.1) {muEtaRegion = 2;}
	else {muEtaRegion = 3;}

	int muPtRegion_Trigger = -1;
	if (pt < 30){muPtRegion_Trigger = 0;}
	else if (pt < 40){muPtRegion_Trigger = 1;}
	else if (pt < 50){muPtRegion_Trigger = 2;}
	else if (pt < 60){muPtRegion_Trigger = 3;}
	else if (pt < 120){muPtRegion_Trigger = 4;}
	else if (pt < 200){muPtRegion_Trigger = 5;}
	else {muPtRegion_Trigger = 6;}

	int muPtRegion_IDIso = -1;
	if (pt < 25){muPtRegion_IDIso = 0;}
	else if (pt < 30){muPtRegion_IDIso = 1;}
	else if (pt < 40){muPtRegion_IDIso = 2;}
	else if (pt < 50){muPtRegion_IDIso = 3;}
	else if (pt < 60){muPtRegion_IDIso = 4;}
	else {muPtRegion_IDIso = 5;}
	
	double muEffSF = muTrackingSF[muTrackEtaRegion][systLevel] * muIdIsoSF[muPtRegion_IDIso][muEtaRegion][systLevel] * muTrigSF[muPtRegion_Trigger][muEtaRegion][systLevel];
	
	return muEffSF;
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
	


	return  weight;
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
double getBtagSF(EventTree* tree, EventPick* evt,string sysType, BTagCalibrationReader reader){
        double prod = 1.0;
        double jetpt;
        double jeteta;
        int jetflavor;
        double SFb;
        double prod1 = 1.0;
        double jetpt1;
        double jeteta1;
        int jetflavor1;
        double SFb1;
        double weight0=1.0;
        double weight1=0.;
        double weight=0.;
        double weight2=0.;
        std::vector<float> weights;
        std::vector<float> ScaleFactors;
        if ( evt->NBjet_ge == 0) return 1.0;
	if (evt-> NBjet_ge == 1){
		for(std::vector<int>::const_iterator bjetInd = evt->bJets.begin(); bjetInd != evt->bJets.end(); bjetInd++){
                jetpt = tree->jetPt_->at(*bjetInd);
                jeteta = fabs(tree->jetEta_->at(*bjetInd));
                jetflavor = abs(tree->jetPartonID_->at(*bjetInd));

                if(jetflavor == 5) {
                        ScaleFactors.push_back(reader.eval_auto_bounds(sysType, BTagEntry::FLAV_B, jeteta, jetpt));

                } else if( jetflavor == 4) {
                } else {
                        ScaleFactors.push_back(reader.eval_auto_bounds(sysType, BTagEntry::FLAV_UDSG, jeteta, jetpt));
                }
		if (ScaleFactors.size() < 2) return 0;
		for (int j = 0; j < ScaleFactors.size();++j){

                        prod=ScaleFactors[j];
                        for (int i =0; i < ScaleFactors.size(); ++i){
                                if (i==j) continue;
                                prod*=(1-ScaleFactors[i]);
                                }
                        weight1 += prod;

                }
		return weight1;
	}				
 	}
        if (evt-> NBjet_ge >= 1) {
        for(std::vector<int>::const_iterator bjetInd = evt->bJets.begin(); bjetInd != evt->bJets.end(); bjetInd++){
                jetpt = tree->jetPt_->at(*bjetInd);
                jeteta = fabs(tree->jetEta_->at(*bjetInd));
                jetflavor = abs(tree->jetPartonID_->at(*bjetInd));

                if(jetflavor == 5) SFb = reader.eval_auto_bounds(sysType, BTagEntry::FLAV_B, jeteta, jetpt);
                else if( jetflavor == 4) SFb = reader.eval_auto_bounds(sysType, BTagEntry::FLAV_C, jeteta, jetpt);
                else SFb = reader.eval_auto_bounds(sysType, BTagEntry::FLAV_UDSG, jeteta, jetpt);

                weight0*= 1.0 - SFb;
        }
        return 1.0 - weight0;
       }
	if (evt->NBjet_ge >= 2) {
                for(std::vector<int>::const_iterator bjetInd = evt->bJets.begin(); bjetInd != evt->bJets.end(); bjetInd++){
                jetpt = tree->jetPt_->at(*bjetInd);
                jeteta = fabs(tree->jetEta_->at(*bjetInd));
                jetflavor = abs(tree->jetPartonID_->at(*bjetInd));

                if(jetflavor == 5) {
                        SFb = reader.eval_auto_bounds(sysType, BTagEntry::FLAV_B, jeteta, jetpt);
                        ScaleFactors.push_back(reader.eval_auto_bounds(sysType, BTagEntry::FLAV_B, jeteta, jetpt));

                } else if( jetflavor == 4) {
                        SFb = reader.eval_auto_bounds(sysType, BTagEntry::FLAV_C, jeteta, jetpt);
                } else {
                        SFb = reader.eval_auto_bounds(sysType, BTagEntry::FLAV_UDSG, jeteta, jetpt);
                        ScaleFactors.push_back(reader.eval_auto_bounds(sysType, BTagEntry::FLAV_UDSG, jeteta, jetpt));
                }

                weight0*= 1.0 - SFb;
                weights.push_back(weight0);
                }
                if (ScaleFactors.size() < 2) return 0;

                for (int j = 0; j < ScaleFactors.size();++j){

                        prod=ScaleFactors[j];
                        for (int i =0; i < ScaleFactors.size(); ++i){
                                if (i==j) continue;
                                prod*=(1-ScaleFactors[i]);
                                }
                        weight1 += prod;

                }
                weights.push_back(weight1);
                weight2 = 1 - weight0 - weight1;
                weights.push_back(weight2);

                return weight2;
        }
}




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

