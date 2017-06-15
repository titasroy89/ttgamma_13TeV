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
	ele_Pt_cut = 30.0;
	ele_Eta = 2.4;
	ele_PtLoose_cut = 10.0;
	
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
	pho_Et_cut = 25.0; 
	pho_ID_ind = 0; // 0 - Loose, 1 - Medium, 2 - Tight
	pho_noPixelSeed_cut = false;
	pho_noEleVeto_cut = false;
	
	// muons :https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
	mu_PtLoose_cut = 10.0;
	mu_RelIsoLoose_cut = 0.25;
        mu_RelIso_range[0] = 0.0;
        mu_RelIso_range[1] = 0.15;
	mu_Iso_invert = false;
	mu_Eta_tight = 2.1;
	mu_Eta_loose = 2.4;
	mu_Pt_cut = 26;
//	bool passPhoMediumID;
}

void Selector::process_objects(const EventTree* inp_tree){
	tree = inp_tree;
	//std::cout << "tree done" << std::endl;
	clear_vectors();
	//std::cout << "vectors done " << std::endl;
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
	Pho03PhoSCRIso.clear();
	Pho03RandPhoIso.clear();
	Pho03RandChHadIso.clear();
}

void Selector::filter_photons(){
	for(int phoInd = 0; phoInd < tree->nPho_; ++phoInd){
		double eta = tree->phoEta_->at(phoInd);
		double et = tree->phoEt_->at(phoInd);

		// uint photonIDbit = tree->phoIDbit_->at(phoInd);
		// bool passLoosePhotonID  = photonIDbit >> 0 & 1;
		// bool passMediumPhotonID = photonIDbit >> 1 & 1;
		// bool passTightPhotonID  = photonIDbit >> 2 & 1;
		Pho03ChHadIso.push_back( max(0.,tree->phoPFChIso_->at(phoInd) - phoEffArea03ChHad(eta)*tree->rho_));
		bool passMediumPhotonID = passPhoMediumID(phoInd);

		bool hasPixelSeed = tree->phohasPixelSeed_->at(phoInd);
		

		int region = 0; //barrel
		if(TMath::Abs( eta )>1.5) region = 1; //endcap
		bool phoPresel = (fidEtaPass(eta) &&
				  et > pho_Et_cut && 
				  passMediumPhotonID && 
				  !hasPixelSeed);// &&
					//	( pho_noPixelSeed_cut || tree->phohasPixelSeed_->at(phoInd) == 0 ) &&
					//	( pho_noEleVeto_cut || tree->phoEleVeto_->at(phoInd) == 0 ) &&
					//	tree->phoSeedBCE_->at(phoInd) == photonID_IsConv[region][pho_ID_ind] &&
					//	tree->phoHoverE_->at(phoInd) < photonID_HoverE[region][pho_ID_ind] &&
					//	Pho03NeuHadIso[phoInd] < 
//			(photonID_RhoCorrR03NeuHadIso_0[region][pho_ID_ind] + et*photonID_RhoCorrR03NeuHadIso_1[region][pho_ID_ind]);
//		std::cout<< "the boolean is :" <<  phoPresel <<std::endl;
//		std::cout<<"passing fidEtaPass:"<< fidEtaPass(eta)<<std::endl;
		//bool phoPresel = true;	
		if(phoPresel){
			PhotonsPresel.push_back(phoInd);
//			std::cout<< "length of photons presel vector" << PhotonsPresel.size() <<std::endl;
			PhoPassSih.push_back( tree->phoSigmaIEtaIEtaFull5x5_->at(phoInd) < photonID_SigmaIEtaIEta[region][pho_ID_ind] );
			// substitute ChHadIso cut to loose SC footprint removed ChHadIso cut of 5 GeV in photon selection
			//PhoPassChHadIso.push_back( Pho03ChHadSCRIso[phoInd] < 5 /*photonID_RhoCorrR03ChHadIso[region][pho_ID_ind]*/ );
		//	PhoPassPhoIso.push_back( Pho03PhoIso[phoInd] < photonID_RhoCorrR03PhoIso_0[region][pho_ID_ind] + et * photonID_RhoCorrR03PhoIso_1[region][pho_ID_ind] );
		}
	}
}

void Selector::filter_electrons(){
	for(int eleInd = 0; eleInd < tree->nEle_; ++eleInd){
//		std::cout << "start to iterate over electron number " << eleInd << " out of " << tree->nEle_ << std::endl; 
		double eta = tree->eleSCEta_->at(eleInd);
		//double eta = tree->eleEta_->at(eleInd);
		double pt = tree->elePt_->at(eleInd);
		//std::cout << " pt of the electron " << pt << std::endl;
		double rho_zero = std::max(0.0, (double)tree->rho_);
		Ele03RelIso.push_back( 
			(tree->elePFChIso_->at(eleInd) + 
			 std::max(0.0, tree->elePFNeuIso_->at(eleInd) + tree->elePFPhoIso_->at(eleInd) - rho_zero * eleEffArea03(eta))
			) / pt );
		
		//bool Iso_pass = ele_RelIso_range[0] <= Ele03RelIso[eleInd] &&
			//	Ele03RelIso[eleInd] < ele_RelIso_range[1];
		//bool MVA_pass = ele_MVA_range[0] < tree->eleIDMVATrig_->at(eleInd) &&
		//		tree->eleIDMVATrig_->at(eleInd) <= ele_MVA_range[1];
		//
	//	std::cout << "eleIDbit " << tree->eleIDbit_->at(eleInd)  << std::endl;
	

		
	//	bool cutbasedID = ( tree->eleIDbit_ & 4) >> 2  ; 
	//	std::cout << "after defining cutbasedID bool" << std::endl;			
		
			
                
		//bool Iso_MVA_pass;
		//if(ele_Iso_MVA_invert) Iso_MVA_pass = !Iso_pass || !MVA_pass;
		//else Iso_MVA_pass = Iso_pass && MVA_pass;
      		//bool eleSel = fidEtaPass( eta ) &&
		//				pt > ele_Pt_cut &&
		//				
		//				tree->eleEtaseedAtVtx_->at(eleInd) == 0 &&
		//				TMath::Abs(tree->eleD0_->at(eleInd)) < ele_Dxy_cut &&
		//				tree->eleMissHits_->at(eleInd) <= ele_MissInnHit_cut;
		
	//	std::cout << "len eta " << tree->eleSCEta_->size() << std::endl;
             	//std::cout << "len eleID " << tree->eleIDbit_->size() << std::endl;
            //    std::cout << "len mvaID " << tree->eleIDMVATrg_->size() << std::endl;
              //  std::cout << "len mvaIDNT " << tree->eleIDMVANonTrg_->size() << std::endl;
                //std::cout << "len pt " << tree->elePt_->size() << std::endl;
		//if (tree->eleIDbit_->size() > 0 ){
		//	std::cout << "what is inside eleIDbit_ actually? " << tree->eleIDbit_
		if (tree->eleIDbit_->size() != tree->elePt_->size()){	
			std::cout << " vectors are of different length skip to next electron " << std::endl;
			continue;
		}
		
		
               bool cutbasedID = false;
               //if (tree->eleIDbit_->size() == tree->elePt_->size()){
                  //std::cout << "cutbasedID bool" << std::endl;
  //             std::cout << " electron ID bit is " << tree->eleIDbit_->at(eleInd) << std::endl;
               cutbasedID = ( tree->eleIDbit_->at(eleInd) >> 2 & 1)  ;
	//       if (cutbasedID ) {
	//		std::cout << " It is passing the ID cut !! " << std::endl;
		   //continue;
         //        }
                //}
          //     std::cout << "after defining cutbasedID bool" << std::endl;
//	       std::cout << "eta is " << fabs(eta) << std::endl;
//	       std::cout << " electron pt is " << pt << std::endl;
//	       std::cout << " cutbasedID pass : " << cutbasedID << std::endl;
	       bool eleSel = fabs( eta )  < 2.4 &&
                                                pt > ele_Pt_cut ;// && 
                                                //cutbasedID  ;
                                                //tree->eleEtaseedAtVtx_->at(eleInd) == 0 &&
                                                //TMath::Abs(tree->eleD0_->at(eleInd)) < ele_Dxy_cut &&
                                                //tree->eleMissHits_->at(eleInd) <= ele_MissInnHit_cut;





	
		bool looseSel = fabs(eta) < 2.5 && 
						pt > ele_PtLoose_cut ;//&& 
						//Ele03RelIso[eleInd] < ele_RelIsoLoose_cut && 
						//tree->eleIDMVATrg_->at(eleInd) > ele_MVALoose_cut;
		
        //       std::cout << "after defining all cuts" << std::endl; 
		
		if( eleSel && cutbasedID ){
			Electrons.push_back(eleInd);
		}
		else if( looseSel ){ 
			ElectronsLoose.push_back(eleInd);
		}
	//	std::cout << " lenght of the vector Electrons " << ElectronsLoose.size() << std::endl;
	}
//	std::cout << " lenght of the vector Electrons " << Electrons.size() << std::endl;
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
	//	bool Iso_pass = mu_RelIso_range[0] <= Mu04RelIso[muInd] &&
         //                       Mu04RelIso[muInd] < mu_RelIso_range[1];

		double rho_zero = std::max(0.0, (double)tree->rho_);
		//Mu04RelIso.push_back( 
		//	(tree->muPFChIso_->at(muInd) +
		//	 std::max(0.0, tree->muPFNeuIso_->at(muInd) + tree->muPFPhoIso_->at(muInd) - rho_zero * muEffArea04(eta))
		//	) / pt );
		Mu04RelIso.push_back( frelIsocorr );
		// fill Muons vector with indices of muons passing selection

		bool IsoPass = frelIsocorr >= mu_RelIso_range[0] && frelIsocorr <= mu_RelIso_range[1];		

		if (mu_Iso_invert) IsoPass = !IsoPass;

		bool passLoose = pt > 10.0 && TMath::Abs(eta) < 2.4 && frelIsocorr < 0.25 && isPFMuon && ( isGlobalMuon || isTrackerMuon);
		bool passTight = pt > 30.0 && TMath::Abs(eta) < 2.1 && 
				frelIsocorr < 0.15 &&
				tree->muChi2NDF_->at(muInd) < 10 &&
				tree->muTrkLayers_->at(muInd) > 5 &&
				tree->muMuonHits_->at(muInd) > 0 &&
				tree->muD0_->at(muInd) < 0.2 &&
				fabs( tree->muDz_->at(muInd) ) < 0.5 && //check this
				tree->muPixelHits_->at(muInd) > 0 &&
				tree->muStations_->at(muInd) > 1 &&
		                IsoPass &&
				isPFMuon && isGlobalMuon && isTrackerMuon;
		
		//bool muSel = TMath::Abs(eta) < 2.5 &&
		//				pt > mu_PtLoose_cut &&
		//				Mu04RelIso[muInd] < mu_RelIsoLoose_cut;

		if(passTight){
		 	Muons.push_back(muInd);
		 }
		else if (passLoose){
		 	MuonsLoose.push_back(muInd);
		 }

//		if( passTight && passLoose ){
//			Muons.push_back(muInd);
//		}
//		if ( passTight && !passLoose){
//			Muons.push_back(muInd);
//		}
//		if( passLoose && !passTight ){
//			MuonsLoose.push_back(muInd);
//		}
	}
 //       std::cout << " lenght of the vector Muons " << Muons.size() << std::endl;	
}

// jet ID is not likely to be altered, so it is hardcoded
void Selector::filter_jets(){
	for(int jetInd = 0; jetInd < tree->nJet_; ++jetInd){
//		std::cout << "starts to iterate with jet number " << jetInd << " out of " << tree->nJet_ << std::endl;
		bool jetID_pass = false;
//		std::cout<< "leading jet Pt " << tree->jetPt_->at(0) << std::endl;		
//		std::cout << " length of jet Pt " << tree->jetPt_->size() << std::endl;
//		std::cout << " length of jetID " << tree->jetPFLooseID_->size() << std::endl;
		if ( tree->jetPt_->size() == tree->jetPFLooseID_->size() ) {
			 jetID_pass = ( tree->jetPFLooseID_->at(jetInd) == 1) ; 
		}
//		std::cout << " jet ID " <<  tree->jetPFLooseID_ << std::endl;
		bool jetPresel = TMath::Abs(tree->jetEta_->at(jetInd)) < 2.4 &&
						 tree->jetPt_->at(jetInd) > 30.0 &&
						 jetID_pass ;
						// tree->AK8JetCHF_->at(jetInd) > 0 &&
						 //tree->AK8JetNHF_->at(jetInd) < 0.99 &&
						//tree->AK8JetCEF_->at(jetInd) < 0.99 &&
						//tree->AK8JetNEF_->at(jetInd) < 0.99 &&
						//tree->AK8JetNCH_->at(jetInd) > 0 &&
						//tree->AK8Jetnconstituents_->at(jetInd) > 1;
//		std::cout << "crosses the AK8 stuff " << std::endl;	
		if( jetPresel){
			Jets.push_back(jetInd);
			if(tree->jetpfCombinedMVAV2BJetTags_->at(jetInd) > btag_cut) bJets.push_back(jetInd);
		}
		
	}
 //    std::cout << " Size of the Jets vector is : " << Jets.size() << std::endl;
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
	static const double area[7] = {0.012, 0.010, 0.014, 0.012, 0.016, 0.020, 0.012};
	return area[phoRegion(eta)];
}
double Selector::phoEffArea03NeuHad(double phoEta){
	double eta = TMath::Abs(phoEta);
	static const double area[7] = {0.030, 0.057, 0.039, 0.015, 0.024, 0.039, 0.072};
	return area[phoRegion(eta)];
}
double Selector::phoEffArea03Pho(double phoEta){
	double eta = TMath::Abs(phoEta);
	static const double area[7] = {0.148, 0.130, 0.112, 0.216, 0.262, 0.260, 0.266};
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
