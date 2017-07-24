#include"EventPick.h"
#include <iostream> 
#include <iomanip>


double secondMinDr(int myInd, const EventTree* tree);

EventPick::EventPick(std::string titleIn){
	title = titleIn;

	cutFlow = new TH1D("cut_flow","cut flow",14,-0.5,13.5);
	cutFlow->SetDirectory(0);
	set_cutflow_labels(cutFlow); // keep the labels close to the cuts definitions (below)
	histVector.push_back(cutFlow);
	
	cutFlowWeight = new TH1D("cut_flow_weight","cut flow with PU weight",14,-0.5,13.5);
	cutFlowWeight->SetDirectory(0);
	set_cutflow_labels(cutFlowWeight);
	histVector.push_back(cutFlowWeight);

	genPhoRegionWeight = new TH1D("genPhoRegionWeight","GenPhoton passing fiducial cuts: barrel 0 or endcap 1",2,-0.5,1.5);
	genPhoRegionWeight->SetDirectory(0);
	histVector.push_back(genPhoRegionWeight);

	genPhoRegionWeight_1l_2l = new TH1D("genPhoRegionWeight_1l_2l","GenPhoton passing fiducial cuts with 1 or 2 gen leptons: barrel 0 or endcap 1",2,-0.5,1.5);
	genPhoRegionWeight_1l_2l->SetDirectory(0);
	histVector.push_back(genPhoRegionWeight_1l_2l);

	genPhoRegionWeight_1fiducial = new TH1D("genPhoRegionWeight_1lfid","GenPhoton passing fiducial cuts with 1 or 2 gen leptons passing fiducial cuts: barrel 0 or endcap 1",8,-0.5,7.5);
	genPhoRegionWeight_1fiducial->GetXaxis()->SetBinLabel(1,"1 e, barrel pho");
	genPhoRegionWeight_1fiducial->GetXaxis()->SetBinLabel(2,"1 e, endcap pho");
	genPhoRegionWeight_1fiducial->GetXaxis()->SetBinLabel(3,"1 e, topSel, barrel pho");
	genPhoRegionWeight_1fiducial->GetXaxis()->SetBinLabel(4,"1 e, topSel, endcap pho");
	genPhoRegionWeight_1fiducial->GetXaxis()->SetBinLabel(5,"1 #mu, barrel pho");
	genPhoRegionWeight_1fiducial->GetXaxis()->SetBinLabel(6,"1 #mu, endcap pho");
	genPhoRegionWeight_1fiducial->GetXaxis()->SetBinLabel(7,"1 #mu, topSel, barrel pho");
	genPhoRegionWeight_1fiducial->GetXaxis()->SetBinLabel(8,"1 #mu, topSel, endcap pho");

	genPhoRegionWeight_1fiducial->SetDirectory(0);
	histVector.push_back(genPhoRegionWeight_1fiducial);

	genPhoMinDR = new TH1D("genPhoMinDR", "Min DR between gen photon and other gen particles", 100, 0., 1.);
	genPhoMinDR->SetDirectory(0);
	histVector.push_back(genPhoMinDR);


	// assign cut values
	veto_jet_dR = 0.1;
	veto_jet_lep_dR = 0.4;
	veto_jet_pho_dR = 0.1;
	veto_pho_jet_dR = 0.4;
	veto_pho_lep_dR = 0.4;
	MET_cut = 20.0;
	no_trigger = false;
	Jet_Pt_cut_1 = 30; // change this to 25.
	Njet_ge = 1;
	NBjet_ge = 1;
	NBjet_tag=1;
	Nele_eq = 1;
	Nmu_eq = 1;
	NEleVeto_le = 0;
	Npho_ge = 1;
	NlooseMuVeto_le = 0;
	NlooseEleVeto_le = 0;
	//NlooseMuVeto_le = 0;
	NmediumEleVeto_le = 0;
}

EventPick::~EventPick(){
}

void EventPick::process_event(const EventTree* inp_tree, const Selector* inp_selector, double weight){
	tree = inp_tree;
	selector = inp_selector;
	clear_vectors();
	passSkim = false;
	passPreSel = false;
	passAll = false;
	for(std::vector<int>::const_iterator phoVi = selector->PhotonsPresel.begin();phoVi != selector->PhotonsPresel.end(); phoVi++){	
	bool goodPhoton = true;
	for(std::vector<int>::const_iterator muInd = selector->Muons.begin(); muInd != selector->Muons.end(); muInd++){
               if(dR_mu_pho(*muInd, *phoVi)< veto_pho_lep_dR) goodPhoton = false; 
        //          std::cout<<"muon pt "<< tree->muPt_->at(*muInd)<< std::endl;
          //        std::cout<< "delR(pho,mu)"<< dR_mu_pho(*muInd,*phoVi) << std::endl;}
         	}
	       if(goodPhoton){
                        Photons.push_back(*phoVi);
                }
	 }
		
//	std::cout <<"number of :"<< selector->Jets.size()<<std::endl;
	for(std::vector<int>::const_iterator jetInd = selector->Jets.begin(); jetInd != selector->Jets.end(); jetInd++){
		bool goodJet = true;
		//remove jets too close to electrons	
		for(std::vector<int>::const_iterator eleInd = selector->Electrons.begin(); eleInd != selector->Electrons.end(); eleInd++)
			if(dR_jet_ele(*jetInd, *eleInd) < veto_jet_lep_dR) goodJet = false;
//		std::cout <<"jet pt is"<< tree->jetPt_->at(*jetInd) << "passes jet ele dr?"<< goodJet <<std::endl;
		//remove jets too close to muons
		for(std::vector<int>::const_iterator muInd = selector->Muons.begin(); muInd != selector->Muons.end(); muInd++)
			if(dR_jet_mu(*jetInd, *muInd) < veto_jet_lep_dR) goodJet = false;
//		std::cout <<"jet pt is"<< tree->jetPt_->at(*jetInd) << "passes jet mu dr?"<< goodJet <<std::endl;
//		// remove jets too close to photons
//		std::cout << "are there photons?how many?" << Photons.size()<<std::endl;
		for(std::vector<int>::const_iterator phoVi = Photons.begin(); phoVi != Photons.end(); phoVi++)
                        if(dR_jet_pho(*jetInd, *phoVi) < veto_jet_pho_dR) goodJet = false;
//		std::cout <<"jet pt is"<< tree->jetPt_->at(*jetInd) << "passes jet pho dr?"<< goodJet <<std::endl;
			
		
				
		if(goodJet) Jets.push_back(*jetInd);
		// take care of bJet collection
		for(std::vector<int>::const_iterator bjetInd = selector->bJets.begin(); bjetInd != selector->bJets.end(); bjetInd++)
			if(*bjetInd == *jetInd && goodJet) bJets.push_back(*bjetInd);
	}
        // std::cout<<"Jetsize"<<Jets.size()<<std::endl;	
	// keep electrons that are not close to jets (veto_lep_jet_dR)
	for(std::vector<int>::const_iterator eleInd = selector->Electrons.begin(); eleInd != selector->Electrons.end(); eleInd++){
		bool goodEle = true;
		if(goodEle) Electrons.push_back(*eleInd);
	}
	
	//loose electrons
	for(std::vector<int>::const_iterator eleInd = selector->ElectronsLoose.begin(); eleInd != selector->ElectronsLoose.end(); eleInd++){
		bool goodEle = true;
		if(goodEle) ElectronsLoose.push_back(*eleInd);
	}
	//medium electrons
	 for(std::vector<int>::const_iterator eleInd = selector->ElectronsMedium.begin(); eleInd != selector->ElectronsMedium.end(); eleInd++){
                bool goodEle = true;
                if(goodEle) ElectronsMedium.push_back(*eleInd);
        }
	 //do cleaning of muons that are close to jets
	for(std::vector<int>::const_iterator muInd = selector->Muons.begin(); muInd != selector->Muons.end(); muInd++){
		bool goodMu = true;
		if(goodMu) Muons.push_back(*muInd);
	}
	
	 //photon cleaning:
//	std::cout<<"how many photons are there?"<< Photons.size() << std::endl;
	for(std::vector<int>::const_iterator phoVi = Photons.begin();phoVi != Photons.end(); phoVi++){
		bool goodPhoton = true;
//		std::cout<<"how many jets are there?"<< Jets.size()<<std::endl;
		for(std::vector<int>::const_iterator jetInd = Jets.begin();jetInd != Jets.end(); jetInd++){
		 //remove photons close to jets
			if (dR_jet_pho(*jetInd, *phoVi) < veto_pho_jet_dR) goodPhoton = false; 
                          // std::cout <<"jet with  Pt "<< tree->jetPt_->at(*jetInd)<<std::endl;
                          // std::cout << "delR(pho,jet)"<<dR_jet_pho(*jetInd, *phoVi) << std::endl;
			} 
			if(goodPhoton) PhotonSel.push_back(*phoVi);
                             
                       
		         
		
	             // }
	         }
	bool Pass_trigger=true;
	Pass_trigger = ( tree->HLTEleMuX >> 20 & 1) || (tree->HLTEleMuX >> 19 & 1) ;
	cutFlow->Fill(0.0); // Input events
	cutFlowWeight->Fill(0.0,weight);
	passPreSel = false;
	passSkim = true;
	//std::cout <<tree->HLTEleMuX<<std::endl;
	if( passSkim && Pass_trigger) {cutFlow->Fill(1); cutFlowWeight->Fill(1,weight);passSkim = true;}
	else passSkim = false;
	if( passSkim && tree->isPVGood_) {cutFlow->Fill(2); cutFlowWeight->Fill(2,weight);}
	else passSkim = false;
	if(passSkim && Muons.size() == Nmu_eq) {cutFlow->Fill(3); cutFlowWeight->Fill(3,weight);}
	else passSkim = false;
	if(passSkim && selector->MuonsLoose.size() <= 0 ) {cutFlow->Fill(4); cutFlowWeight->Fill(4,weight);}
	else passSkim = false;
	if(passSkim && selector->Electrons.size() <= 0. && selector->ElectronsLoose.size() <= 0.) {cutFlow->Fill(5); cutFlowWeight->Fill(5,weight);}
	else passSkim = false;
	if(passSkim && Jets.size() >= Njet_ge ) {cutFlow->Fill(6); cutFlowWeight->Fill(6,weight);}
	else passSkim = false;
	if(passSkim && Jets.size() >= Njet_ge+1 ) {cutFlow->Fill(7); cutFlowWeight->Fill(7,weight);}
        else passSkim = false;
	if(passSkim && Jets.size() >= Njet_ge+2 ) {cutFlow->Fill(8); cutFlowWeight->Fill(8,weight);}//std::cout << tree->event_ <<std::endl;}
        else passSkim = false;
	if ( passSkim &&bJets.size() >= NBjet_ge) {cutFlow->Fill(9); cutFlowWeight->Fill(9,weight);passPreSel=true;}
        else passSkim = false;
//	if(passPreSel && Jets.size() >= Njet_ge+3) {cutFlow->Fill(10); cutFlowWeight->Fill(10,weight);}
//	else passPreSel = false;	
        if ( passPreSel && bJets.size() == NBjet_ge  ) {cutFlow->Fill(10); cutFlowWeight->Fill(10,weight);}
        else passPreSel = false;
	if(passPreSel && PhotonSel.size() >= Npho_ge) { cutFlow->Fill(11); cutFlowWeight->Fill(11,weight);passAll = true;}// std::cout <<"Event which passes photon sel:"<<tree->event_ <<std::endl;}
	
	else passAll = false ; 
//	std::cout <<"done with cutflow"<<std::endl;	
	// require >=1 photon
	//if(passPreSel && Photons.size() >= Npho_ge){
	//	cutFlow->Fill(8);
	//	cutFlowWeight->Fill(8,weight);
	//	passAll = true;
//	}

	// saving information about Gen Level photons, if any
	// Save it only if PreSelection passed
	// Separate count for barrel and endcap (will be used separately anyway)
	
	bool foundGenPhotonBarrel = false;
	bool foundGenPhotonEndcap = false;
	if(passPreSel && !(tree->isData_)){
		for(int mcInd=0; mcInd<tree->nMC_; ++mcInd){
			if(tree->mcPID->at(mcInd) == 22 &&
			   (tree->mcParentage->at(mcInd)==2 || tree->mcParentage->at(mcInd)==10 || tree->mcParentage->at(mcInd)==26) &&
			   tree->mcPt->at(mcInd) > selector->pho_Et_cut){		
				if(secondMinDr(mcInd, tree) > 0.2){
					double fabsEta = TMath::Abs(tree->mcEta->at(mcInd));
					if(fabsEta < 1.4442){ 
						foundGenPhotonBarrel = true;
					}
					if( 1.566 < fabsEta && fabsEta < 2.5){
						 foundGenPhotonEndcap = true;
					}
				}
			}
		}
	}
	if(foundGenPhotonBarrel) genPhoRegionWeight->Fill(0.0, weight);
	if(foundGenPhotonEndcap) genPhoRegionWeight->Fill(1.0, weight);

	int EleP = 0;
	int EleM = 0;
	int MuP = 0;
	int MuM = 0;
	int TauP = 0;
	int TauM = 0;

	for( int mcI = 0; mcI < tree->nMC_; ++mcI){
	  if(abs(tree->mcMomPID->at(mcI))==24 && tree->mcParentage->at(mcI)==10){
	    if( tree->mcPID->at(mcI) == 11 ) EleP = 1;
	    if( tree->mcPID->at(mcI) == -11 ) EleM = 1;
	    if( tree->mcPID->at(mcI) == 13 ) MuP = 1;
	    if( tree->mcPID->at(mcI) == -13 ) MuM = 1;
	    if( tree->mcPID->at(mcI) == 15) TauP = 1;
	    if( tree->mcPID->at(mcI) == -15) TauM = 1;
	  }
	}
	int nEle = EleP + EleM;
	int nMu = MuP + MuM;
	int nTau = TauP + TauM;
	int nLep = nEle + nMu + nTau;

	if (nLep == 1 || nLep == 2){
	  if(foundGenPhotonBarrel) genPhoRegionWeight_1l_2l->Fill(0.0, weight);
	  if(foundGenPhotonEndcap) genPhoRegionWeight_1l_2l->Fill(1.0, weight);
	}	  

	
	int ElePfid = 0;
	int EleMfid = 0;
	int MuPfid = 0;
	int MuMfid = 0;
	int nNufid = 0;
	for( int mcI = 0; mcI < tree->nMC_; ++mcI){
	  if((abs(tree->mcMomPID->at(mcI))==24 && tree->mcParentage->at(mcI)==10) || (abs(tree->mcMomPID->at(mcI))==15 && tree->mcParentage->at(mcI)==26)){		  
	    if( tree->mcPID->at(mcI) == 11 ) {
	      if (tree->mcPt->at(mcI) > 35 && (fabs(tree->mcEta->at(mcI)) < 2.5 && !(fabs(tree->mcEta->at(mcI)) > 1.4442 && fabs(tree->mcEta->at(mcI))<1.566))) ElePfid += 1;
	    }
	    if( tree->mcPID->at(mcI) == -11 ) {
	      if (tree->mcPt->at(mcI) > 35 && (fabs(tree->mcEta->at(mcI)) < 2.5 && !(fabs(tree->mcEta->at(mcI)) > 1.4442 && fabs(tree->mcEta->at(mcI))<1.566))) EleMfid += 1;
	    }
	    if( tree->mcPID->at(mcI) == 13 ) {
	      if (tree->mcPt->at(mcI) > 26 && fabs(tree->mcEta->at(mcI)) < 2.1) MuPfid += 1;
	    }
	    if( tree->mcPID->at(mcI) == -13 ) {
	      if (tree->mcPt->at(mcI) > 26 && fabs(tree->mcEta->at(mcI)) < 2.1) MuMfid += 1;
	    }
	  }
	  if( fabs(tree->mcPID->at(mcI)) == 12 || fabs(tree->mcPID->at(mcI)) == 14 || fabs(tree->mcPID->at(mcI)) == 16 ) {
	    if (tree->mcPt->at(mcI) > 20) nNufid += 1;
	  }
	}
	int nElefid = ElePfid + EleMfid;
	int nMufid = MuPfid + MuMfid;
	int nJetsfid = 0;
	if ((nElefid + nMufid)==1 && nNufid == 1){
	  for ( int jetI = 0; jetI < tree->nJet_; jetI++){
	    if (tree->jetGenPt_->at(jetI) >= 30) {
	      if ( fabs(tree->jetGenEta_->at(jetI)) < 2.4) nJetsfid += 1;
	    }
	  }
	}
	

	if (nElefid==1 && nMufid==0){
	  if(foundGenPhotonBarrel) genPhoRegionWeight_1fiducial->Fill(0.0, weight);
	  if(foundGenPhotonEndcap) genPhoRegionWeight_1fiducial->Fill(1.0, weight);
	  if(nJetsfid >=3 && nNufid == 1){
	    if (foundGenPhotonBarrel) genPhoRegionWeight_1fiducial->Fill(2.0, weight);
	    if (foundGenPhotonEndcap) genPhoRegionWeight_1fiducial->Fill(3.0, weight);
	  }
	}	  
	if (nElefid==0 && nMufid==1){
	  if(foundGenPhotonBarrel) genPhoRegionWeight_1fiducial->Fill(4.0, weight);
	  if(foundGenPhotonEndcap) genPhoRegionWeight_1fiducial->Fill(5.0, weight);
	  if(nJetsfid >=3 && nNufid == 1){
	    if (foundGenPhotonBarrel) genPhoRegionWeight_1fiducial->Fill(6.0, weight);
	    if (foundGenPhotonEndcap) genPhoRegionWeight_1fiducial->Fill(7.0, weight);
	  }
	}	  

	if(passPreSel && !(tree->isData_)){
	  double minDR = 999.;
	  for(int mcInd=0; mcInd<tree->nMC_; ++mcInd){
	    if(tree->mcPID->at(mcInd) == 22 &&
	       (tree->mcParentage->at(mcInd)==2 || tree->mcParentage->at(mcInd)==10 || tree->mcParentage->at(mcInd)==26)){
	      double dr = secondMinDr(mcInd, tree);
	      if (dr < minDR) minDR = dr;
	    }
	  }
	  genPhoMinDR->Fill(minDR, weight);
	}

}

void EventPick::print_cutflow(){
	std::cout << "Cut-Flow for the event selector: " << title << std::endl;
	std::cout << "Input Events :                " << cutFlow->GetBinContent(1) << std::endl;
	std::cout << "Passing Trigger               " << cutFlow->GetBinContent(2) << std::endl;
	std::cout << "Has Good Vtx               " << cutFlow->GetBinContent(3) << std::endl;
	std::cout << "Events with = " << Nmu_eq << " muon     " << cutFlow->GetBinContent(4) << std::endl;
	std::cout << "Events with <= " << NlooseMuVeto_le << " loose muons " << cutFlow->GetBinContent(5) << std::endl;
	std::cout << "Events with <= " << NEleVeto_le << " electrons " << cutFlow->GetBinContent(6) << std::endl;
	std::cout << "Events with >= " << Njet_ge << " jets   " << cutFlow->GetBinContent(7) << std::endl;
	std::cout << "Events with >= " << Njet_ge+1 << " jets " << cutFlow->GetBinContent(8) << std::endl;
 	std::cout << "Events with >= " << Njet_ge+2 << " jets " << cutFlow->GetBinContent(9) << std::endl;
	std::cout << "Events with >= " << NBjet_ge << " bjets " << cutFlow->GetBinContent(10) << std::endl;
//	std::cout << "Events with >= " << Njet_ge+3 << " jets "<< cutFlow->GetBinContent(11) << std::endl;
        std::cout << "Events with >= " << NBjet_ge << " bjets "<< cutFlow->GetBinContent(11) << std::endl;
	std::cout << "Events with >= 1 photon " << cutFlow->GetBinContent(12) << std::endl; //change this later
	std::cout << std::endl;
}

void EventPick::set_cutflow_labels(TH1D* hist){
	hist->GetXaxis()->SetBinLabel(1,"Input");
	hist->GetXaxis()->SetBinLabel(2,"passTrigger");
	hist->GetXaxis()->SetBinLabel(3,"hasGoodVtx");
	hist->GetXaxis()->SetBinLabel(4,"==1Mu");
	hist->GetXaxis()->SetBinLabel(5,"0 LooseMu");
 	hist->GetXaxis()->SetBinLabel(6,"0 Electrons");
	hist->GetXaxis()->SetBinLabel(7,">=1jet");
	hist->GetXaxis()->SetBinLabel(8,">=2jets");
	hist->GetXaxis()->SetBinLabel(9,">=3jets");
	hist->GetXaxis()->SetBinLabel(10,">=1 btags");
//	hist->GetXaxis()->SetBinLabel(11,">=4jets");
	hist->GetXaxis()->SetBinLabel(11,">=1 btags");
	hist->GetXaxis()->SetBinLabel(12,">=1 Photon"); //change this later
	hist->GetXaxis()->SetBinLabel(1,"");
}

void EventPick::clear_vectors(){
	Electrons.clear();
	ElectronsLoose.clear();
	ElectronsMedium.clear();
	Muons.clear();
	MuonsLoose.clear();
	Jets.clear();
	bJets.clear();
	Photons.clear();
	PhotonSel.clear();
	PhoPassChHadIso.clear();
	PhoPassPhoIso.clear();
	PhoPassSih.clear();
}

double EventPick::dR_jet_ele(int jetInd, int eleInd){
	return dR(tree->jetEta_->at(jetInd), tree->jetPhi_->at(jetInd), tree->eleEta_->at(eleInd), tree->elePhi_->at(eleInd));
}
double EventPick::dR_jet_mu(int jetInd, int muInd){
	return dR(tree->jetEta_->at(jetInd), tree->jetPhi_->at(jetInd), tree->muEta_->at(muInd), tree->muPhi_->at(muInd));
}
double EventPick::dR_jet_pho(int jetInd, int phoInd){
	return dR(tree->jetEta_->at(jetInd), tree->jetPhi_->at(jetInd), tree->phoEta_->at(phoInd), tree->phoPhi_->at(phoInd));
}
double EventPick::dR_ele_pho(int eleInd, int phoInd){
	return dR(tree->eleEta_->at(eleInd), tree->elePhi_->at(eleInd), tree->phoEta_->at(phoInd), tree->phoPhi_->at(phoInd));
}
double EventPick::dR_mu_pho(int muInd, int phoInd){
	return dR(tree->muEta_->at(muInd), tree->muPhi_->at(muInd), tree->phoEta_->at(phoInd), tree->phoPhi_->at(phoInd));
}

