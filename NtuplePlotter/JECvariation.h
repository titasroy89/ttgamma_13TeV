#ifndef JECVARIATION_H
#define JECVARIATION_H
        
#include<vector>
#include<iostream>
#include<algorithm>
#include<TH1F.h>
#include<TMath.h>
#include<TLorentzVector.h>
#include"/uscms_data/d3/troy2012/cmssw/CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include"/uscms_data/d3/troy2012/cmssw/CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include"/uscms_data/d3/troy2012/cmssw/CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include"EventTree.h"

class JECvariation{
public: 
	//JetCorrectionUncertainty();
    	//JetCorrectionUncertainty(const std::string& fDataFile);
    	//JetCorrectionUncertainty(const JetCorrectorParameters& fParameters);
    	//~JetCorrectionUncertainty();

	//JetCorrectorParameters() { valid_ = false;}
    	//JetCorrectorParameters(const std::string& fFile, const std::string& fSection = "");
    	//JetCorrectorParameters(const JetCorrectorParameters::Definitions& fDefinitions,
          //               const std::vector<JetCorrectorParameters::Record>& fRecords)
      	//: mDefinitions(fDefinitions),mRecords(fRecords) { valid_ = true;}


	JECvariation(std::string inputPrefix, bool isMC);
        ~JECvariation();

        void applyJEC(EventTree* tree, int scaleDownNormUp012);

private:
        JetCorrectionUncertainty *jecUnc;

        JetCorrectorParameters *ResJetPar;
        JetCorrectorParameters *L3JetPar;
        JetCorrectorParameters *L2JetPar;
	JetCorrectorParameters *L1JetPar;

        FactorizedJetCorrector *JetCorrector;
};


#endif
