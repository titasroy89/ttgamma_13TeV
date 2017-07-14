#include"Selector.h"


double dR(double eta1, double phi1, double eta2, double phi2){
	double dphi = phi2 - phi1;
	double deta = eta2 - eta1;
	static const double pi = TMath::Pi();
	dphi = TMath::Abs( TMath::Abs(dphi) - pi ) - pi;
	return TMath::Sqrt( dphi*dphi + deta*deta );
}

Selector::Selector(){
	// jets
	jet_Pt_cut = 30;
	btag_cut = 0.8484;

	// electrons
	ele_Pt_cut = 35.0;
	ele_Eta_tight = 2.1;
	ele_Eta_loose = 2.5;
	ele_PtLoose_cut = 15.0;
	
	ele_Ptmedium_cut = 20.0;
	ele_RelIso_range[0] = 0.0;
	ele_RelIso_range[1] = 0.1;
	ele_RelIsoLoose_cut = 0.2;
	ele_MVA_range[0] = 0.9;
	ele_MVA_range[1] = 1.0;
	ele_MVALoose_cut = 0.0;
        ele_cutbased_range[0] =0.9 ;
	ele_cutbased_range[1] =1.0 ;
	ele_Dxy_cut = 0.02;
	ele_MissInnHit_cut = 0;
	ele_Iso_MVA_invert = false;

	// photons
	pho_Et_cut = 15.0; //check if this works 
	pho_ID_ind = 0; // 0 - Loose, 1 - Medium, 2 - Tight
	pho_noPixelSeed_cut = false;
	pho_noEleVeto_cut = false;
	
	// muons :https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
	mu_PtLoose_cut = 15.0;
	mu_RelIsoLoose_cut = 0.25;
        mu_RelIso_range[0] = 0.0;
        mu_RelIso_range[1] = 0.15;
	mu_Iso_invert = false;
	mu_Eta_tight = 2.4;
	mu_Eta_loose = 2.4;
	mu_Pt_cut = 30;
//	bool passPhoMediumID;
}

void Selector::process_objects(const EventTree* inp_tree){
//	std::cout << "starting selector" << std::endl;
	tree = inp_tree;
//	std::cout << tree << std::endl;
	clear_vectors();
//	std::cout << "vectors done " << std::endl;
	filter_photons();
//	std::cout << "photons done " << std::endl;
        filter_muons();
  //      std::cout << "muons done" << std::endl;

 	filter_electrons();
//	std::cout << " electrons done" << std::endl;
	
	filter_jets();
  //      std::cout << "jets done" << std::endl;

}

void Selector::clear_vectors(){
	PhotonsPresel.clear();
	PhoPassChHadIso.clear();
	PhoPassPhoIso.clear();
	PhoPassSih.clear();
	
	Electrons.clear();
	ElectronsLoose.clear();
	ElectronsMedium.clear();
	Muons.clear();
	MuonsLoose.clear();
	Jets.clear();
	bJets.clear();
	
	Ele03RelIso.clear();
	Mu04RelIso.clear();
	Pho03ChHadIso.clear();
	Pho03ChHadSCRIso.clear();
	Pho03NeuHadIso.clear();
	Pho03PhoIso.clear();
	//Pho03PhoSCRIso.clear();
//	Pho03RandPhoIso.clear();
//	Pho03RandChHadIso.clear();
}

void Selector::filter_photons(){
//	std::cout << "starting to select photons" << std::endl;

	for(int phoInd = 0; phoInd < tree->nPho_; ++phoInd){
		double eta = tree->phoEta_->at(phoInd);
		double et = tree->phoEt_->at(phoInd);
	//	 std::cout << "starting to select photons" << std::endl;
		// uint photonIDbit = tree->phoIDbit_->at(phoInd);
		// bool passLoosePhotonID  = photonIDbit >> 0 & 1;
		// bool passMediumPhotonID = photonIDbit >> 1 & 1;
		// bool passTightPhotonID  = photonIDbit >> 2 & 1;
		Pho03ChHadIso.push_back( max(0.,tree->phoPFChIso_->at(phoInd) - phoEffArea03ChHad(eta)*tree->rho_));
		Pho03NeuHadIso.push_back(max(0.,tree->phoPFNeuIso_->at(phoInd) - phoEffArea03NeuHad(eta)*tree->rho_));
		Pho03PhoIso.push_back(max(0.,tree->phoPFPhoIso_->at(phoInd) - phoEffArea03Pho(eta)*tree->rho_));
	//	bool passMediumPhotonID = passPhoMediumID(phoInd);

		bool hasPixelSeed = tree->phohasPixelSeed_->at(phoInd);
		bool passMediumPhotonID = false;
		//implementing photon ID bit for medium cut
		passMediumPhotonID = ( tree->phoIDbit_->at(phoInd) >> 1 & 1)  ;
		int region = 0; //barrel
		if(TMath::Abs( eta )>1.5) region = 1; //endcap
		bool phoPresel = (fidEtaPass(eta) &&
				  et > pho_Et_cut && 
				  passMediumPhotonID && 
				  !hasPixelSeed);
		if(phoPresel){
			PhotonsPresel.push_back(phoInd);
			PhoPassSih.push_back( tree->phoSigmaIEtaIEtaFull5x5_->at(phoInd) < photonID_SigmaIEtaIEta[region][pho_ID_ind] );
		}
	}
}

void Selector::filter_electrons(){
	for(int eleInd = 0; eleInd < tree->nEle_; ++eleInd){
		double eta = tree->eleSCEta_->at(eleInd);
		double pt = tree->elePt_->at(eleInd);
		double rho_zero = std::max(0.0, (double)tree->rho_);
		Ele03RelIso.push_back( 
			(tree->elePFChIso_->at(eleInd) + 
			 std::max(0.0, tree->elePFNeuIso_->at(eleInd) + tree->elePFPhoIso_->at(eleInd) - rho_zero * eleEffArea03(eta))
			) / pt );
		
		
		
	       bool eleSel = false;
	       bool looseSel = false;			
               bool cutbasedtightID = false;
	       bool cutbasedmediumID = false;
	       bool cutbasedlooseID = false;
               //ID cut for electron
               cutbasedtightID = ( tree->eleIDbit_->at(eleInd) >> 3 & 1)  ;//tight cut
	       cutbasedmediumID = ( tree->eleIDbit_->at(eleInd) >> 2 & 1)  ;//medium cut	
	       cutbasedlooseID = ( tree->eleIDbit_->at(eleInd) >> 1 & 1)  ;//loose cut
	       if (fabs( eta )>1.5){	
	       		eleSel = fabs( eta )  < ele_Eta_tight &&
                                                pt > ele_Pt_cut && tree->eleD0_->at(eleInd) <0.10 && tree->eleDz_->at(eleInd) <0.20 ;
			looseSel = fabs(eta) < ele_Eta_loose &&
                                                pt > ele_PtLoose_cut && tree->eleD0_->at(eleInd) <0.10 && tree->eleDz_->at(eleInd) <0.20 ;
			}
	       else{
			eleSel = fabs( eta )  < ele_Eta_tight &&
                                                pt > ele_Pt_cut && tree->eleD0_->at(eleInd) <0.05 && tree->eleDz_->at(eleInd) <0.10 ;
			looseSel = fabs(eta) < ele_Eta_loose && pt > ele_PtLoose_cut && tree->eleD0_->at(eleInd) < 0.05 && tree->eleDz_->at(eleInd) <0.10;
	
		}

	
	       bool passEtaEBEEGap = (fabs( eta ) < 1.4442) || (fabs( eta ) > 1.566);		
		
	       if( eleSel && cutbasedtightID && passEtaEBEEGap){
			Electrons.push_back(eleInd);
		}
		else if( looseSel && cutbasedlooseID && passEtaEBEEGap){ 
			ElectronsLoose.push_back(eleInd);
		}
	}
}
void Selector::filter_muons(){
	for(int muInd = 0; muInd < tree->nMu_; ++muInd){
		const unsigned int GlobalMuon     =  1<<1;
		const unsigned int TrackerMuon    =  1<<2;
		const unsigned int PFMuon =  1<<5;
		bool isGlobalMuon  = tree->muType_->at(muInd) & GlobalMuon;
		bool isTrackerMuon = tree->muType_->at(muInd) & TrackerMuon;
		bool isPFMuon      = tree->muType_->at(muInd) & PFMuon;
		double eta = tree->muEta_->at(muInd);
		double pt = tree->muPt_->at(muInd);
		double frelIsocorr = ( tree->muPFChIso_->at(muInd) + 
					fmax(0.0, tree->muPFNeuIso_->at(muInd) + 
						tree->muPFPhoIso_->at(muInd) -
						0.5*tree->muPFPUIso_->at(muInd)
					) 
				     ) / pt;
		Mu04RelIso.push_back( frelIsocorr );
		//muon ID cut
	        bool MuontightID = false;
		bool MuonlooseID = false;
		MuontightID = ( tree->muIDbit_->at(muInd) >> 2 & 1); 
		MuonlooseID = ( tree->muIDbit_->at(muInd) >> 0 & 1);

		bool IsoPass = frelIsocorr >= mu_RelIso_range[0] && frelIsocorr <= mu_RelIso_range[1];		

	//	if (mu_Iso_invert) IsoPass = !IsoPass;

		bool passLoose = pt > mu_PtLoose_cut && TMath::Abs(eta) < 2.4 && frelIsocorr < 0.25 &&  MuonlooseID;//isPFMuon && ( isGlobalMuon || isTrackerMuon);
	
		bool passTight = pt > mu_Pt_cut && TMath::Abs(eta) < 2.4 && 
				frelIsocorr < 0.15 &&  MuontightID;
	//			tree->muChi2NDF_->at(muInd) < 10 &&
	//			tree->muTrkLayers_->at(muInd) > 5 &&
	//			tree->muMuonHits_->at(muInd) > 0 &&
	//			fabs(tree->muD0_->at(muInd)) < 0.2 && // use absolute values
	//			fabs( tree->muDz_->at(muInd) ) < 0.5 &&
	//			tree->muPixelHits_->at(muInd) > 0 &&
	//			tree->muStations_->at(muInd) > 1 &&
	//	                IsoPass &&
	//			isPFMuon && isGlobalMuon && isTrackerMuon;
		

		if(passTight){
		 	Muons.push_back(muInd);
		 }
		else if (passLoose){
		 	MuonsLoose.push_back(muInd);
		 }

	}
}

void Selector::filter_jets(){
	for(int jetInd = 0; jetInd < tree->nJet_; ++jetInd){
		bool jetID_pass = false;
		//if ( tree->jetPt_->size() == tree->jetPFLooseID_->size() ) {
		//	 jetID_pass = ( tree->jetPFLooseID_->at(jetInd) == 1) ; 
	//	}
		jetID_pass = ( tree->jetID_->at(jetInd) >> 2 & 1) ;//tight jet ID
		bool jetPresel = TMath::Abs(tree->jetEta_->at(jetInd)) < 2.4 &&
						 tree->jetPt_->at(jetInd) > 30.0 &&
						 jetID_pass ;
		if( jetPresel){
			Jets.push_back(jetInd);
			if(tree->jetCSV2BJetTags_->at(jetInd) > btag_cut) bJets.push_back(jetInd);
		}
		
	}
}

bool Selector::fidEtaPass(double Eta){
	double fabsEta = TMath::Abs(Eta);
	if( fabsEta > 2.5) return false;
	if( 1.4442 < fabsEta && fabsEta < 1.566) return false;
	return true;
}

// kEleGammaAndNeutralHadronIso03
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h
double Selector::eleEffArea03(double SCEta){
	double eta = TMath::Abs(SCEta);
	static const double area[7] = {0.130, 0.137, 0.067, 0.089, 0.107, 0.110, 0.138};
	int region = 0;
	if( eta >= 1.0 )   region++;
	if( eta >= 1.479 ) region++;
	if( eta >= 2.0 )   region++;
	if( eta >= 2.2 )   region++;
	if( eta >= 2.3 )   region++;
	if( eta >= 2.4 )   region++;
	return area[region];
}

// kMuGammaAndNeutralHadronIso04
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/sixie/Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h?revision=1.7&view=markup
double Selector::muEffArea04(double muEta){
	double eta = TMath::Abs(muEta);
	static const double area[6] = {0.674, 0.565, 0.442, 0.515, 0.821, 0.660};
	int region = 0;
	if( eta >= 1.0 )   region++;
	if( eta >= 1.479 ) region++;
	if( eta >= 2.0 )   region++;
	if( eta >= 2.2 )   region++;
	if( eta >= 2.3 )   region++;
	return area[region];
}

//https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
//https://indico.cern.ch/event/491548/contributions/2384977/attachments/1377936/2117789/CutBasedPhotonID_25-11-2016.pdf
int Selector::phoRegion(double absEta){
	int region = 0;
	if( absEta >= 1.0  ) region++;
	if( absEta >= 1.479) region++;
	if( absEta >= 2.0  ) region++;
	if( absEta >= 2.2  ) region++;
	if( absEta >= 2.3  ) region++;
	if( absEta >= 2.4  ) region++;
	return region;
}
double Selector::phoEffArea03ChHad(double phoEta){
	double eta = TMath::Abs(phoEta);
	static const double area[7] = {0.0360, 0.0377, 0.0306, 0.0283, 0.0254, 0.0217, 0.0167};
	return area[phoRegion(eta)];
}
double Selector::phoEffArea03NeuHad(double phoEta){
	double eta = TMath::Abs(phoEta);
	static const double area[7] = {0.0597, 0.0807, 0.0629, 0.0197,0.0184, 0.0284, 0.0591};
	return area[phoRegion(eta)];
}
double Selector::phoEffArea03Pho(double phoEta){
	double eta = TMath::Abs(phoEta);
	static const double area[7] = {0.1210, 0.1107, 0.0699, 0.1056, 0.1457, 0.1719, 0.1998};
	return area[phoRegion(eta)];
}

bool Selector::passPhoMediumID(int phoInd){
  double pt = tree->phoEt_->at(phoInd);  
  double eta = TMath::Abs(tree->phoEta_->at(phoInd));
  bool passMediumID = false;

  double rhoCorrPFChIso = tree->phoPFChIso_->at(phoInd) - phoEffArea03ChHad(eta)*tree->rho_;
  double rhoCorrPFNeuIso = tree->phoPFNeuIso_->at(phoInd) - phoEffArea03NeuHad(eta)*tree->rho_;
  double rhoCorrPFPhoIso = tree->phoPFPhoIso_->at(phoInd) - phoEffArea03Pho(eta)*tree->rho_;
  

  if (eta < 1.47){
    if (tree->phoHoverE_->at(phoInd)                < 0.0396  &&
	tree->phoSigmaIEtaIEtaFull5x5_->at(phoInd)  < 0.01022 &&
	rhoCorrPFChIso                             < 0.441 &&
	rhoCorrPFNeuIso                             < 2.725+0.0148*pt+0.000017*pt*pt &&
	rhoCorrPFPhoIso                             < 2.571+0.0047*pt){
      passMediumID = true;
    }
  } else {
    if (tree->phoHoverE_->at(phoInd)                < 0.0219  &&
	tree->phoSigmaIEtaIEtaFull5x5_->at(phoInd)  < 0.03001 &&
	rhoCorrPFChIso                             < 0.442 &&
	rhoCorrPFNeuIso                             < 1.715+0.0163*pt+0.000014*pt*pt &&
	rhoCorrPFPhoIso                             < 3.863+0.0034*pt){
      passMediumID = true;
    }
  }
  return passMediumID;
}


Selector::~Selector(){
}
