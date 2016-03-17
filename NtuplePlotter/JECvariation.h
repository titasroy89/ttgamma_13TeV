#ifndef JECVARIATION_H
#define JECVARIATION_H
        
#include<vector>
#include<iostream>
#include<algorithm>
#include<TH1F.h>
#include<TMath.h>
#include<TLorentzVector.h>
#include"JetMETObjects/JetCorrectorParameters.h"
#include"JetMETObjects/FactorizedJetCorrector.h"
#include"JetMETObjects/JetCorrectionUncertainty.h"
#include"EventTree.h"

class JECvariation{
public:
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
