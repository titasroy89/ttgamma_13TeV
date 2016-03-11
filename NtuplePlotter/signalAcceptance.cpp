#include<iostream>
#include<cmath>

#include"EventTree.h"
#include"Selector.h"
#include"EventPick.h"
#include"TVector2.h"

void saveHist(TH1F* hist, TFile* file);
void doJER(EventTree* tree);
double JERcorrection(double eta);
int minDrIndex(double myEta, double myPhi, std::vector<float> *etas, std::vector<float> *phis);
int secondMinDrIndex(int myInd, EventTree* tree);

bool overlapMadGraph(EventTree* tree);

double dPhi(double dphi){	
	return fabs( fabs( fabs(dphi) - M_PI ) - M_PI);
}

void fillCategory(EventTree* tree, TH1F* hist, double weight){
	int EleP = 0;
	int EleM = 0;
	int MuP = 0;
	int MuM = 0;
	int TauP = 0;
	int TauM = 0;
	int ElePfid = 0;
	int EleMfid = 0;
	int MuPfid = 0;
	int MuMfid = 0;
	int nNufid = 0;

	int nPhofid = 0;
	int nJetfid = 0;
	int nBJetfid = 0;

	TVector2 MET = TVector2(0,0);
	TVector2 tempNu = TVector2(0,0);

	for( int mcI = 0; mcI < tree->nMC_; ++mcI){
	  if(abs(tree->mcMomPID->at(mcI))==24 && tree->mcParentage->at(mcI)==10){
	    if( tree->mcPID->at(mcI) == 11 ) EleP = 1;
	    if( tree->mcPID->at(mcI) == -11 ) EleM = 1;
	    if( tree->mcPID->at(mcI) == 13 ) MuP = 1;
	    if( tree->mcPID->at(mcI) == -13 ) MuM = 1;
	    if( tree->mcPID->at(mcI) == 15) TauP = 1;
	    if( tree->mcPID->at(mcI) == -15) TauM = 1;
	  }
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
	    if (tree->mcPt->at(mcI) > 20){
	      nNufid += 1;
	      //cout << nNufid << "\t" << tree->run_ << tree->lumis_ << tree->event_ << "\t" << mcI << "\t" <<  tree->mcPt->at(mcI) << "\t" << tree->mcEta->at(mcI) << "\t" << tree->mcPhi->at(mcI) << endl;
	    }
	    tempNu.SetMagPhi(tree->mcPt->at(mcI),tree->mcPhi->at(mcI));
	    MET += tempNu;
	    //	    cout << nNufid << "\t" << tree->run_ << tree->lumis_ << tree->event_ << "\t" << mcI << "\t" << tree->mcPID->at(mcI) << "\t" <<  tree->mcPt->at(mcI) << "\t" << tree->mcEta->at(mcI) << "\t" << tree->mcPhi->at(mcI) << "\t" << MET.Mod() << endl;
	  }

	  if(tree->mcPID->at(mcI) == 22 && 
	     (tree->mcParentage->at(mcI)==2 || tree->mcParentage->at(mcI)==10 || tree->mcParentage->at(mcI)==26) && 
	     tree->mcPt->at(mcI) > 25 && 
	     fabs(tree->mcEta->at(mcI)) < 1.4442){
	    nPhofid += 1;
	  }

	  if( abs(tree->mcPID->at(mcI)) < 6 && (abs(tree->mcMomPID->at(mcI))==24 || abs(tree->mcMomPID->at(mcI))==6 ) && tree->mcPt->at(mcI) > 30 && abs(tree->mcEta->at(mcI))<2.4 ) {
	    nJetfid += 1;
	    if (abs(tree->mcMomPID->at(mcI))==6 && abs(tree->mcPID->at(mcI))==5){
	      nBJetfid += 1;
	    }
	  }
	}
	
	hist->Fill(1.0, weight); // Total
	int nEle = EleP + EleM;
	int nMu = MuP + MuM;
	int nTau = TauP + TauM;
	if( nEle + nMu + nTau == 0) hist->Fill(2.0, weight); // All Had
	if( nEle + nMu + nTau == 1) hist->Fill(3.0, weight); // Single Lepton
	if( nEle + nMu + nTau == 2) hist->Fill(4.0, weight); // Di Lepton

	int nElefid = ElePfid + EleMfid;
	int nMufid = MuPfid + MuMfid;

	if ( nEle==1 && nMu==0 && nTau==0 ) hist->Fill(6.0, weight);
	if ( nEle==0 && nMu==1 && nTau==0 ) hist->Fill(7.0, weight);
	if ( nEle==0 && nMu==0 && nTau==1 ) hist->Fill(8.0, weight);
	if ( nElefid==1 && nMufid==0) hist->Fill(9.0, weight); //Ejets final state (generated)
	if ( nElefid==0 && nMufid==1) hist->Fill(10.0, weight); //Mujets final state (generated)

	// if ( nElefid==1 && nMufid==0) hist->Fill(6.0, weight); //Ejets final state (generated)
	// if ( nElefid==2 && nMufid==0) hist->Fill(7.0, weight); //Ejets final state (generated)
	// if ( nElefid==0 && nMufid==1) hist->Fill(8.0, weight); //Mujets final state (generated)
	// if ( nElefid==0 && nMufid==2) hist->Fill(9.0, weight); //Ejets final state (generated)
	// if ( nElefid==1 && nMufid==1) hist->Fill(10.0, weight); //Ejets final state (generated)

	int nJetsfid = 0;
	int nBJetsfid = 0;
	
	if ((nElefid + nMufid)==1){
	  // for ( int jetI = 0; jetI < tree->genJetPt_->size(); jetI++){
	  //   if (tree->genJetPt_->at(jetI) >= 30 && fabs(tree->genJetEta_->at(jetI)) < 2.4) nJetsfid += 1;	    
	  // }
	  for ( int jetI = 0; jetI < tree->nJet_; jetI++){
	    //	    if (tree->jetGenPt_->at(jetI) >= 30 && fabs(tree->jetGenEta_->at(jetI)) < 2.4) nJetsfid += 1;
	    if (tree->jetGenJetPt_->at(jetI) >= 30 && fabs(tree->jetGenEta_->at(jetI)) < 2.4){
	      nJetsfid += 1;
	      if (abs(tree->jetGenPartonID_->at(jetI))==5) nBJetsfid += 1;
	    }
	  }
	}


	if(nElefid==1 && nMufid==0 && nJetsfid >=3){
	  hist->Fill(11.0, weight);
	  if (nBJetsfid >= 1){
	    hist->Fill(12.0, weight);
	    if (MET.Mod() > 20){
	      //	    if (tree->genMET_ > 20){
	      hist->Fill(13.0, weight);
	      if (nPhofid > 0) hist->Fill(14.0, weight);
	    }
	  }
	}


	// //TESTING
	// if(nElefid==0 && nMufid==1 && nJetsfid >=3){
	//   hist->Fill(11.0, weight);
	//   if (nNufid==1){
	//     hist->Fill(12.0,weight);
	//   }
	//   if (nNufid > 0){
	//     hist->Fill(13.0, weight);
	//   }
	//   if (MET.Mod() > 20) hist->Fill(14.0,weight);
	//   if (nBJetsfid >= 1){
	//     //	    hist->Fill(12.0, weight);
	//     if (tree->genMET_ > 20){
	//       //	      hist->Fill(13.0, weight);
	//       // if (nPhofid > 0) hist->Fill(14.0, weight);
	//     }
	//   }
	// }

	if(nElefid==0 && nMufid==1 && nJetsfid >=3){
	  hist->Fill(16.0, weight);
	  if (nBJetsfid >= 1){
	    hist->Fill(17.0, weight);
	    //	    if (tree->genMET_ > 20){
	    if (MET.Mod() > 20){
	      hist->Fill(18.0, weight);
	      if (nPhofid > 0) hist->Fill(19.0, weight);
	    }
	  }
	}
	// if(nElefid==0 && nMufid==1 && nJetfid >=3 && nBJetfid >= 1 && nNufid == 1) {
	//   hist->Fill(12.0, weight);
	//   if (nPhofid > 0) hist->Fill(15.0, weight);
	// }

	// if(nElefid==1 && nMufid==0 && nJetsfid >=3 && nNufid == 1){
	//   hist->Fill(11.0, weight);
	//   if (nPhofid > 0) hist->Fill(14.0, weight);
	// }
	// if(nElefid==0 && nMufid==1 && nJetsfid >=3 && nNufid == 1) {
	//   hist->Fill(12.0, weight);
	//   if (nPhofid > 0) hist->Fill(15.0, weight);
	// }

	return;
}

int main(int ac, char** av){
	if(ac < 2){
		std::cout << "usage: ./signalAcceptance inputFile[s]" << std::endl;
		return -1;
	}
	
	TH1F* allCategory = new TH1F("allCategory","all Category",19,0.5,19.5);
	allCategory->GetXaxis()->SetBinLabel(1,"Total");
	allCategory->GetXaxis()->SetBinLabel(2,"AllHad");
	allCategory->GetXaxis()->SetBinLabel(3,"1 lepton");
	allCategory->GetXaxis()->SetBinLabel(4,"2 leptons");
	allCategory->GetXaxis()->SetBinLabel(5,"");
	allCategory->GetXaxis()->SetBinLabel(6,"1 e");
	allCategory->GetXaxis()->SetBinLabel(7,"2 e");
	allCategory->GetXaxis()->SetBinLabel(8,"1 #mu");
	allCategory->GetXaxis()->SetBinLabel(9,"2 #mu");
	allCategory->GetXaxis()->SetBinLabel(10,"1 e 1 #mu");
	allCategory->GetXaxis()->SetBinLabel(11,"1 e, 3 jets");
	allCategory->GetXaxis()->SetBinLabel(12,"1 e, 3 jets 1 b");
	allCategory->GetXaxis()->SetBinLabel(13,"1 e, 3 jets 1 b, MET");
	allCategory->GetXaxis()->SetBinLabel(14,"1 e, 3 jets 1 b, MET, 1 pho");
	allCategory->GetXaxis()->SetBinLabel(16,"1 #mu, 3 jets");
	allCategory->GetXaxis()->SetBinLabel(17,"1 #mu, 3 jets 1 b");
	allCategory->GetXaxis()->SetBinLabel(18,"1 #mu, 3 jets 1 b, MET");
	allCategory->GetXaxis()->SetBinLabel(19,"1 #mu, 3 jets 1 b, MET, 1 pho");

	TH1F* preselCategory = new TH1F("preselCategory","presel Category",19,0.5,19.5);
	preselCategory->GetXaxis()->SetBinLabel(1,"Total");
	preselCategory->GetXaxis()->SetBinLabel(2,"AllHad");
	preselCategory->GetXaxis()->SetBinLabel(3,"1 lepton");
	preselCategory->GetXaxis()->SetBinLabel(4,"2 leptons");
	preselCategory->GetXaxis()->SetBinLabel(5,"");
	preselCategory->GetXaxis()->SetBinLabel(6,"1 e");
	preselCategory->GetXaxis()->SetBinLabel(7,"2 e");
	preselCategory->GetXaxis()->SetBinLabel(8,"1 #mu");
	preselCategory->GetXaxis()->SetBinLabel(9,"2 #mu");
	preselCategory->GetXaxis()->SetBinLabel(10,"1 e 1 #mu");
	preselCategory->GetXaxis()->SetBinLabel(11,"1 e, 3 jets");
	preselCategory->GetXaxis()->SetBinLabel(12,"1 e, 3 jets 1 b");
	preselCategory->GetXaxis()->SetBinLabel(13,"1 e, 3 jets 1 b, MET");
	preselCategory->GetXaxis()->SetBinLabel(14,"1 e, 3 jets 1 b, MET, 1 pho");
	preselCategory->GetXaxis()->SetBinLabel(16,"1 #mu, 3 jets");
	preselCategory->GetXaxis()->SetBinLabel(17,"1 #mu, 3 jets 1 b");
	preselCategory->GetXaxis()->SetBinLabel(18,"1 #mu, 3 jets 1 b, MET");
	preselCategory->GetXaxis()->SetBinLabel(19,"1 #mu, 3 jets 1 b, MET, 1 pho");
	

	TH1F* photonCategory = new TH1F("photonCategory","reco photon Category",19,0.5,19.5);
	photonCategory->GetXaxis()->SetBinLabel(1,"Total");
	photonCategory->GetXaxis()->SetBinLabel(2,"AllHad");
	photonCategory->GetXaxis()->SetBinLabel(3,"1 lepton");
	photonCategory->GetXaxis()->SetBinLabel(4,"2 leptons");
	photonCategory->GetXaxis()->SetBinLabel(5,"");
	photonCategory->GetXaxis()->SetBinLabel(6,"1 e");
	photonCategory->GetXaxis()->SetBinLabel(7,"2 e");
	photonCategory->GetXaxis()->SetBinLabel(8,"1 #mu");
	photonCategory->GetXaxis()->SetBinLabel(9,"2 #mu");
	photonCategory->GetXaxis()->SetBinLabel(10,"1 e 1 #mu");
	photonCategory->GetXaxis()->SetBinLabel(11,"1 e, 3 jets");
	photonCategory->GetXaxis()->SetBinLabel(12,"1 e, 3 jets 1 b");
	photonCategory->GetXaxis()->SetBinLabel(13,"1 e, 3 jets 1 b, MET");
	photonCategory->GetXaxis()->SetBinLabel(14,"1 e, 3 jets 1 b, MET, 1 pho");
	photonCategory->GetXaxis()->SetBinLabel(16,"1 #mu, 3 jets");
	photonCategory->GetXaxis()->SetBinLabel(17,"1 #mu, 3 jets 1 b");
	photonCategory->GetXaxis()->SetBinLabel(18,"1 #mu, 3 jets 1 b, MET");
	photonCategory->GetXaxis()->SetBinLabel(19,"1 #mu, 3 jets 1 b, MET, 1 pho");
	
	TH1F* VisAllCategory = new TH1F("VisAllCategory","all Category, Vis",15,0.5,15.5);
	TH1F* VisPreselCategory = new TH1F("VisPreselCategory","presel Category, Vis",15,0.5,15.5);
	TH1F* VisPhotonCategory = new TH1F("VisPhotonCategory","reco photon Category, Vis",15,0.5,15.5);
	
	
	TH1F* dROtherGen = new TH1F("dROtherGen", "dROtherGen", 800, 0.0, 4.0);
	TH1F* parentage = new TH1F("parentage","parentage",30, 0, 30);	
	TH1F* dptOverpt = new TH1F("dptOverpt","dptOverpt", 400, -2.0, 2.0);
	TH1F* dRrecoGen = new TH1F("dRrecoGen","dRrecoGen", 200, 0.0, 0.2);
	TH1F* dPhiRecoGen = new TH1F("dPhiRecoGen","dPhiRecoGen", 400, 0.0, 0.2);
	TH1F* dEtaRecoGen = new TH1F("dEtaRecoGen","dEtaRecoGen", 800, -0.2, 0.2);
	
	TH1F* dRGenNearJet = new TH1F("dRGenNearJet","dRGenNearJet", 200, 0.0, 1.0);
	TH1F* dPhiGenNearJet = new TH1F("dPhiGenNearJet","dPhiGenNearJet", 100, 0.0, 0.5);
	TH1F* dEtaGenNearJet = new TH1F("dEtaGenNearJet","dEtaGenNearJet", 200, -0.5, 0.5);

	//TH1F* dRGenNextNearJet = new TH1F("dRGenNextNearJet","dRGenNextNearJet", 600, 0.0, 6.0);
	//TH1F* dPhiGenNextNearJet = new TH1F("dPhiGenNextNearJet","dPhiGenNextNearJet", 300, 0.0, 3.0);
	//TH1F* dEtaGenNextNearJet = new TH1F("dEtaGenNextNearJet","dEtaGenNextNearJet", 600, -3.0, 3.0);

	// object selector
	Selector* selectorLoose = new Selector();
	// create event selectors here
	EventPick* evtPickLoose = new EventPick("LoosePhotonID");
	// do not do jet to photon dR cleaning
	evtPickLoose->veto_pho_jet_dR = 0.0;
	
	
	EventTree* tree = new EventTree(ac-1, av+1);
	double PUweight = 1.0;
	
	Long64_t nEntr = tree->GetEntries();
	for(Long64_t entry=0; entry<nEntr; entry++){
		if(entry%10000 == 0) std::cout << "processing entry " << entry << " out of " << nEntr << std::endl;
		tree->GetEntry(entry);
		
		doJER(tree);

		selectorLoose->process_objects(tree);
		evtPickLoose->process_event(tree, selectorLoose, PUweight);
		// fill the histograms
		fillCategory(tree, allCategory, PUweight);
		if(evtPickLoose->passPreSel) fillCategory(tree, preselCategory, PUweight);
		if(evtPickLoose->passAll) fillCategory(tree, photonCategory, PUweight);

		// fill histograms for gen photon passing the acceptance cuts defined in analysis
		bool inAcc = false;
		for(int mcInd=0; mcInd<tree->nMC_; ++mcInd){
			if(tree->mcPID->at(mcInd) == 22 && 
			(tree->mcParentage->at(mcInd)==2 || tree->mcParentage->at(mcInd)==10 || tree->mcParentage->at(mcInd)==26) && 
			tree->mcPt->at(mcInd) > 25 && 
			fabs(tree->mcEta->at(mcInd)) < 1.4442){
				inAcc = true;
			}
		}
		if(inAcc){
			fillCategory(tree, VisAllCategory, PUweight);
			if(evtPickLoose->passPreSel) fillCategory(tree, VisPreselCategory, PUweight);
			if(evtPickLoose->passAll) fillCategory(tree, VisPhotonCategory, PUweight);
		}
		
		// have at least one good photon
		if(!evtPickLoose->passAll) continue;

		// test
		//if(overlapMadGraph(tree)) continue;

		int phoInd = evtPickLoose->Photons.at(0);
		// experiment with delta R cuts for photons
		for(int mcInd=0; mcInd<tree->nMC_; ++mcInd){
			bool etetamatch = dR(tree->mcEta->at(mcInd),tree->mcPhi->at(mcInd),tree->phoEta_->at(phoInd),tree->phoPhi_->at(phoInd)) < 0.2 && 
			(fabs(tree->phoEt_->at(phoInd) - tree->mcPt->at(mcInd)) / tree->mcPt->at(mcInd)) < 1.0;
			if( etetamatch && tree->mcPID->at(mcInd) == 22){
				// test
				if(!(tree->mcParentage->at(mcInd)==2 || tree->mcParentage->at(mcInd)==10 || tree->mcParentage->at(mcInd)==26)) continue;
				
				// fill histograms for mathced photon candidate
				parentage->Fill( tree->mcParentage->at(mcInd) );
				dptOverpt->Fill( (tree->phoEt_->at(phoInd) - tree->mcPt->at(mcInd)) / tree->mcPt->at(mcInd));
				dRrecoGen->Fill( dR(tree->mcEta->at(mcInd),tree->mcPhi->at(mcInd),tree->phoEta_->at(phoInd),tree->phoPhi_->at(phoInd)) );
				dPhiRecoGen->Fill( dPhi( tree->phoPhi_->at(phoInd) - tree->mcPhi->at(mcInd) ) );
				dEtaRecoGen->Fill( tree->phoEta_->at(phoInd) - tree->mcEta->at(mcInd) );
				
				int closestGenInd = secondMinDrIndex( mcInd, tree );
				if(dR(tree->mcEta->at(mcInd), tree->mcPhi->at(mcInd), tree->mcEta->at(closestGenInd), tree->mcPhi->at(closestGenInd)) < 0.01){
					std::cout << "closest PID " << tree->mcPID->at(closestGenInd) << "  MomPID " << tree->mcMomPID->at(closestGenInd) << std::endl;
					std::cout << "photon mother PID " << tree->mcMomPID->at(mcInd) << std::endl;
				}
				dROtherGen->Fill( dR(tree->mcEta->at(mcInd), tree->mcPhi->at(mcInd), tree->mcEta->at(closestGenInd), tree->mcPhi->at(closestGenInd)) );
	
				int closestJetInd = minDrIndex( tree->mcEta->at(mcInd), tree->mcPhi->at(mcInd), tree->jetEta_, tree->jetPhi_ );
				dRGenNearJet->Fill( dR(tree->mcEta->at(mcInd), tree->mcPhi->at(mcInd), tree->jetEta_->at(closestJetInd), tree->jetPhi_->at(closestJetInd) ) );
				dPhiGenNearJet->Fill( dPhi( tree->jetPhi_->at(closestJetInd) - tree->mcPhi->at(mcInd) ) );
				dEtaGenNearJet->Fill( tree->jetEta_->at(closestJetInd) - tree->mcEta->at(mcInd) );
				
				//closestJetInd = secondMinDrIndex( tree->mcEta->at(mcInd), tree->mcPhi->at(mcInd), tree->jetEta_, tree->jetPhi_ );	
				//dRGenNextNearJet->Fill( dR(tree->mcEta->at(mcInd), tree->mcPhi->at(mcInd), tree->jetEta_->at(closestJetInd), tree->jetPhi_->at(closestJetInd) ) );
				//dPhiGenNextNearJet->Fill( dPhi( tree->jetPhi_->at(closestJetInd) - tree->mcPhi->at(mcInd) ) );
				//dEtaGenNextNearJet->Fill( tree->jetEta_->at(closestJetInd) - tree->mcEta->at(mcInd) );
			}
		}
	}

	evtPickLoose->print_cutflow();

	// write histograms
	TFile outFile("signalAcc.root","RECREATE");
	
	saveHist(allCategory, &outFile);
	saveHist(preselCategory, &outFile);
	saveHist(photonCategory, &outFile);
	
	saveHist(VisAllCategory, &outFile);
	saveHist(VisPreselCategory, &outFile);
	saveHist(VisPhotonCategory, &outFile);
	
	saveHist(dROtherGen, &outFile);
	saveHist(parentage, &outFile);
	saveHist(dptOverpt, &outFile);
	saveHist(dRrecoGen, &outFile);
	saveHist(dPhiRecoGen, &outFile);
	saveHist(dEtaRecoGen, &outFile);
	saveHist(dRGenNearJet, &outFile);
	saveHist(dPhiGenNearJet, &outFile);
	saveHist(dEtaGenNearJet, &outFile);
	
	//saveHist(dRGenNextNearJet, &outFile);
	//saveHist(dPhiGenNextNearJet, &outFile);
	//saveHist(dEtaGenNextNearJet, &outFile);

	outFile.Close();

	delete tree;
	return 0;
}

void saveHist(TH1F* hist, TFile* file){
	hist->SetDirectory(file->GetDirectory(""));
	hist->Write();
	hist->SetDirectory(0);
}

int minDrIndex(double myEta, double myPhi, std::vector<float> *etas, std::vector<float> *phis){
	double mindr = 999.0;
	double dr;
	int bestInd = -1;
	for( int oind = 0; oind < etas->size(); oind++){
		dr = dR(myEta, myPhi, etas->at(oind), phis->at(oind));
		if( mindr > dr ) {
			mindr = dr;
			bestInd = oind;
		}
	}
	return bestInd;
}

int secondMinDrIndex(int myInd, EventTree* tree){
	double myEta = tree->mcEta->at(myInd);
	double myPhi = tree->mcPhi->at(myInd);
	int myPID = tree->mcPID->at(myInd);
	
	double mindr = 999.0;
	double dr;
	int bestInd = -1;
	for( int oind = 0; oind < tree->nMC_; oind++){
		if(oind == myInd) continue;
		if(tree->mcMass->at(oind) > 10.0) continue;
		int opid = abs(tree->mcPID->at(oind));
		if(opid == 12 || opid == 14 || opid == 16) continue;
		dr = dR(myEta, myPhi, tree->mcEta->at(oind), tree->mcPhi->at(oind));
		if( mindr > dr ) {
			mindr = dr;
			bestInd = oind;
		}
	}
	return bestInd;
}




// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
void doJER(EventTree* tree){
	// correct MET when jets are smeared
	TLorentzVector tMET;
	tMET.SetPtEtaPhiM(tree->pfMET_,0.0,tree->pfMETPhi_,0.0);
	//std::cout << "before correction MET " << tree->pfMET_ << " " << tree->pfMETPhi_ << "    ";
	// scale jets
	for(int jetInd = 0; jetInd < tree->nJet_ ; ++jetInd){
		if(tree->jetPt_->at(jetInd) < 20) continue;
		if(tree->jetGenJetIndex_->at(jetInd)>0){
			TLorentzVector tjet;
			tjet.SetPtEtaPhiM(tree->jetPt_->at(jetInd), tree->jetEta_->at(jetInd), tree->jetPhi_->at(jetInd), 0.0);
			tMET+=tjet;
			double oldPt = tree->jetPt_->at(jetInd);
			double genPt = tree->jetGenJetPt_->at(jetInd);
			double eta = tree->jetEta_->at(jetInd);
			tree->jetPt_->at(jetInd) = std::max(0.0, genPt + JERcorrection(eta)*(oldPt-genPt));
			tjet.SetPtEtaPhiM(tree->jetPt_->at(jetInd), tree->jetEta_->at(jetInd), tree->jetPhi_->at(jetInd), 0.0);			
			tMET-=tjet;
			//std::cout << "old " << oldPt << "  new " << tree->jetPt_->at(jetInd) << std::endl;
		}
	}
	// save updated MET values
	tree->pfMET_ = tMET.Pt();
	tree->pfMETPhi_ = tMET.Phi();
	//std::cout << "after corrections " << tree->pfMET_ << " " << tree->pfMETPhi_ << std::endl;
}


double JERcorrection(double JetEta){
	double eta = TMath::Abs(JetEta);
	static const double corr[5] = {1.052, 1.057, 1.096, 1.134, 1.288};
	static const double corrDown[5] = {0.990, 1.001, 1.032, 1.042, 1.089};
	static const double corrUp[5] = {1.115, 1.114, 1.161, 1.228, 1.488};
	int region = 0;
	if( eta >= 0.5 ) region++;
	if( eta >= 1.1 ) region++;
	if( eta >= 1.7 ) region++;
	if( eta >= 2.3 ) region++;
	return corr[region];
	
	// should not get here
	return 1.0;
}

