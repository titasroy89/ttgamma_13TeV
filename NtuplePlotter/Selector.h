#ifndef SELECTOR_H
#define SELECTOR_H

#include<vector>
#include<iostream>
#include<algorithm>
#include<TH1F.h>
#include<TMath.h>
#include<TLorentzVector.h>
#include"EventTree.h"

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonID2012
// photon ID is not going to be changed every time this code runs
// barrel/endcap, Loose/Medium/Tight
const int    photonID_IsConv[2][3]                = { {0, 0, 0} ,             {0, 0, 0}             };
const double photonID_HoverE[2][3]                = { {0.05, 0.05, 0.05} ,    {0.05, 0.05, 0.05}    };
const double photonID_SigmaIEtaIEta[2][3]         = { {0.012, 0.011, 0.011} , {0.034, 0.033, 0.031} };
const double photonID_RhoCorrR03ChHadIso[2][3]    = { {2.6, 1.5, 0.7} ,       {2.3, 1.2, 0.5}       };
const double photonID_RhoCorrR03NeuHadIso_0[2][3] = { {3.5, 1.0, 0.4} ,       {2.9, 1.5, 1.5}       };
const double photonID_RhoCorrR03NeuHadIso_1[2][3] = { {0.04, 0.04, 0.04} ,    {0.04, 0.04, 0.04}    };
const double photonID_RhoCorrR03PhoIso_0[2][3]    = { {1.3, 0.7, 0.5} ,       {999, 1.0, 1.0}       };
const double photonID_RhoCorrR03PhoIso_1[2][3]    = { {0.005, 0.005, 0.005} , {0.005, 0.005, 0.005} };

double dR(double eta1, double phi1, double eta2, double phi2);

class Selector{
public:
	Selector();
	~Selector();
	
	void process_objects(const EventTree* inp_tree);
	
	// selected object indices
	std::vector<int> PhotonsPresel;
	std::vector<bool> PhoPassChHadIso;
	std::vector<bool> PhoPassPhoIso;
	std::vector<bool> PhoPassSih;
	std::vector<int> Electrons;
	std::vector<int> ElectronsLoose;
	std::vector<int> ElectronsMedium;
	std::vector<int> Muons;
	std::vector<int> MuonsLoose;
	std::vector<int> Jets;
	std::vector<int> bJets;
	
	// calculated rho corrected PF isolations
	std::vector<double> Ele03RelIso;
	std::vector<double> Mu04RelIso;
	std::vector<double> Pho03ChHadIso;
	std::vector<double> Pho03ChHadSCRIso;
	std::vector<double> Pho03NeuHadIso;
	std::vector<double> Pho03PhoIso;
	std::vector<double> Pho03PhoSCRIso;
	std::vector<double> Pho03RandPhoIso;
	std::vector<double> Pho03RandChHadIso;
	
	// jets
	double jet_Pt_cut;
	double btag_cut;

	// electrons
	double ele_Pt_cut;
	double ele_PtLoose_cut;
	double ele_Eta_tight;
	double ele_Eta_loose;
	double mu_Eta_loose;
	double mu_Eta_tight;
	double mu_Pt_cut;
	double ele_Ptmedium_cut;
	double ele_RelIso_range[2];
	double ele_RelIsoLoose_cut;
	double ele_MVA_range[2];
	double ele_cutbased_range[2];
	double ele_MVALoose_cut;
	double ele_Dxy_cut;
	int    ele_MissInnHit_cut;
	bool   ele_Iso_MVA_invert;
	
	// photons
	double pho_Et_cut;
	int    pho_ID_ind; // 0 - Loose, 1 - Medium, 2 - Tight
	bool   pho_noPixelSeed_cut;
	bool   pho_noEleVeto_cut;

	// muons
	double mu_PtLoose_cut;
	double mu_RelIsoLoose_cut;
	double mu_RelIso_range[2];
 	double mu_MVA_range[2];
	bool   mu_Iso_invert;

private:
	const EventTree* tree;
	void clear_vectors();
	void filter_photons();
	void filter_electrons();
	void filter_muons();
	void filter_jets();
	
	bool fidEtaPass(double Eta);
	
	// effective areas, see Selector.cpp for more information
	double eleEffArea03(double SCEta);
	double muEffArea04(double muEta);
	double phoEffArea03ChHad(double phoEta);
	double phoEffArea03NeuHad(double phoEta);
	double phoEffArea03Pho(double phoEta);
	int phoRegion(double absEta);
	bool passPhoMediumID(int phoInd);
};

#endif
