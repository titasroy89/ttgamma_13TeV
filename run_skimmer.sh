#!/bin/bash

job=$1
cd ${_CONDOR_SCRATCH_DIR}
echo ${_CONDOR_SCRATCH_DIR}
echo "tar -xvf CMSSW_8_0_26_patch1.tar.gz"
tar -xzf CMSSW_8_0_26_patch1.tar.gz
cd CMSSW_8_0_26_patch1/src/
source /cvmfs/cms.cern.ch/cmsset_default.sh

eval `scramv1 runtime -sh`
cd  TTGammaSemiLep/






outputdir="root://cmseos.fnal.gov//store/user/troy2012/skims_13TeV_V08_26/"

outputfiles=("QCD_120to170EM_skim.root" \
"QCD_120to170Mu_skim.root" \
"QCD_170to300EM_skim.root" \
"QCD_170to300Mu_skim.root" \
"QCD_20to30EM_skim.root" \
"QCD_20to30Mu_skim.root" \
"QCD_300to470Mu_skim.root" \
"QCD_300toInfEM_skim.root" \
"QCD_30to50EM_skim.root" \
"QCD_30to50Mu_skim.root" \
"QCD_470to600Mu_skim.root" \
"QCD_50to80EM_skim.root" \
"QCD_50to80Mu_skim.root" \
"QCD_600to800Mu_skim.root" \
"QCD_80to120EM_skim.root" \
"QCD_80to120Mu_skim.root" \
"QCD_800to1000Mu_skim.root" \
"QCD_1000toInfMu_skim.root" \
"Zgamma_skim.root" \
"Wgamma_skim.root" \
"Data_SingleMu_b_skim.root" \
"Data_SingleMu_c_skim.root" \
"Data_SingleMu_d_skim.root" \
"Data_SingleMu_e_skim.root" \
"Data_SingleMu_f1_skim.root" \
"Data_SingleMu_f2_skim.root" \
"Data_SingleMu_g_skim.root" \
"Data_SingleMu_hv2_skim.root" \
"Data_SingleMu_hv3_skim.root" \
"TTbar_skim.root" \
"ST_s-channel_skim.root" \
"ST_tW-channel_skim.root" \
"ST_tbarW-channel_skim.root" \
"ST_t-channel_skim.root" \
"ST_tbar-channel_skim.root" \
"W1jets_skim.root" \
"W2jets_skim.root" \
"W3jets_skim.root" \
"W4jets_skim.root" \
"TTGamma_Dilepton_skim.root" \
"TTGamma_SingleLeptFromTbar_skim.root" \
"TTGamma_Hadronic_skim.root" \
"TTGamma_SingleLeptFromT_skim.root" \
"TTW_skim.root" \
"TTZ_skim.root" \
"DYJets_skim.root")


inputfiles=("QCD_120to170EM_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_120to170EM.root" \
"QCD_120to170Mu_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_120to170Mu.root" \
"QCD_170to300EM_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_170to300EM.root" \
"QCD_170to300Mu_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_170to300Mu.root" \
"QCD_20to30EM_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_20to30EM.root" \
"QCD_20to30Mu_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_20to30Mu.root" \
"QCD_300to470Mu_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_300to470Mu.root" \
"QCD_300toInfEM_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_300toInfEM.root" \
"QCD_30to50EM_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_30to50EM.root" \
"QCD_30to50Mu_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_30to50Mu.root" \
"QCD_470to600Mu_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_470to600Mu.root" \
"QCD_50to80EM_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_50to80EM.root" \
"QCD_50to80Mu_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_50to80Mu.root" \
"QCD_600to800Mu_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_600to800Mu.root" \
"QCD_80to120EM_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_80to120EM.root" \
"QCD_80to120Mu_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_80to120Mu.root" \
"QCD_800to1000Mu_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_800to1000Mu.root" \
"QCD_1000toInfMu_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_1000toInfMu.root" \
"Zgamma_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ZGTo2LG_aMCatNLO_summer16.root" \
"Wgamma_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/WGTo2LG_aMCatNLO_summer16.root" \
"Data_SingleMu_b_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016B_FebReminiAOD.root" \
"Data_SingleMu_c_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016C_FebReminiAOD.root" \
"Data_SingleMu_d_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016D_FebReminiAOD.root" \
"Data_SingleMu_e_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016E_FebReminiAOD.root" \
"Data_SingleMu_f1_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016F_FebReminiAOD1.root" \
"Data_SingleMu_f2_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016F_FebReminiAOD2.root" \
"Data_SingleMu_g_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016G_FebReminiAOD.root" \
"Data_SingleMu_hv2_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016H_FebReminiAODv2.root" \
"Data_SingleMu_hv3_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016H_FebReminiAODv3.root" \
"TTbar_skim.root root://cmseos.fnal.gov//store/user/yumiceva/ntuples_2016/TTbar.root" \
"ST_s-channel_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ST_s.root" \
"ST_tW-channel_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ST_tW_top.root" \
"ST_tbarW-channel_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ST_t_tWbar.root" \
"ST_t-channel_skim.root root://cmseos.fnal.gov//store/user/yumiceva/ntuples_2016/ST_t_top.root" \
"ST_tbar-channel_skim.root root://cmseos.fnal.gov//store/user/yumiceva/ntuples_2016/ST_t_bar.root" \
"W1jets_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/W1jets.root" \
"W2jets_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/W2jets.root" \
"W3jets_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/W3jets.root" \
"W4jets_skim.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/W4jets.root" \
"TTGamma_Dilepton_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ttgamma_dilept.root" \
"TTGamma_SingleLeptFromTbar_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ttgamma_SingleLeptFromTbar.root" \
"TTGamma_Hadronic_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ttgamma_hadronic.root" \
"TTGamma_SingleLeptFromT_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ttgamma_semileptfromT.root" \
"TTW_skim.root root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/TTWJets.root" \
"TTZ_skim.root root://cmseos.fnal.gov//store/user/yumiceva/ntuples_2016/TTZ.root" \
"DYJets_skim.root root://cmseos.fnal.gov///store/user/yumiceva/ntuples_2016//ntuple_DYjets.root")

echo "NtuplePlotter/makeSkim ${inputfiles[job]}"
NtuplePlotter/makeSkim ${inputfiles[job]}

echo "xrdcp -f ${outputfiles[job]} ${outputdir}"
xrdcp -f ${outputfiles[job]} ${outputdir}

