#!/bin/bash

job=$1
cd ${_CONDOR_SCRATCH_DIR}
echo ${_CONDOR_SCRATCH_DIR}
echo "tar -xvf CMSSW_8_0_11.tar.gz"
tar -xzf CMSSW_8_0_11.tar.gz
cd CMSSW_8_0_11/src/
source /cvmfs/cms.cern.ch/cmsset_default.sh

eval `scramv1 runtime -sh`
cd  TTGammaSemiLep/






outputdir="root://cmseos.fnal.gov//store/user/troy2012/QCD/"

outputfiles=("data_mu_b.root" \
"data_mu_c.root" \
"data_mu_d.root" \
"DYJets.root" \
"ST_s.root" \
"ST_t_bar.root" \
"ST_t_top.root" \
"ST_tW_antitop.root" \
"ttgamma.root" \
"Wjets.root" \
"TTbar_1.root")


inputfiles=("data_mu_b.root root://cmsxrootd.fnal.gov//store/user/troy2012/ntuples_2016/job_SingleMu_Run2016B_PRv2.root" \
"data_mu_c.root root://cmsxrootd.fnal.gov//store/user/troy2012/ntuples_2016/job_SingleMu_Run2016C_PRv2.root" \
"data_mu_d.root root://cmsxrootd.fnal.gov//store/user/troy2012/ntuples_2016/job_SingleMu_Run2016D_PRv2.root" \
"DYJets.root root://cmsxrootd.fnal.gov//store/user/troy2012/ntuples_2016/job_spring16_DYJetsToLL_m50_aMCatNLO.root" \
"ST_s.root root://cmsxrootd.fnal.gov//store/user/troy2012/ntuples_2016/ST_s.root" \
"ST_t_bar.root root://cmsxrootd.fnal.gov//store/user/troy2012/ntuples_2016/ST_t_antitop.root" \
"ST_t_top.root root://cmsxrootd.fnal.gov//store/user/troy2012/ntuples_2016/ST_t_top.root" \
"ST_tW_antitop.root root://cmsxrootd.fnal.gov//store/user/troy2012/ntuples_2016/ST_tW_antitop.root" \
"ttgamma.root root://cmsxrootd.fnal.gov//store/user/troy2012/ntuples_2016/ttgamma.root" \
"Wjets.root root://cmsxrootd.fnal.gov//store/user/troy2012/ntuples_2016/job_spring16_WJetsToLNu_aMCatNLO.root" \
"TTbar_1.root root://cmsxrootd.fnal.gov//store/user/troy2012/ntuples_2016/job_spring16_TT_powheg_ext3.root")

echo "NtuplePlotter/makeSkim ${inputfiles[job]}"
NtuplePlotter/makeSkim ${inputfiles[job]}

echo "xrdcp -f ${outputfiles[job]} ${outputdir}"
xrdcp -f ${outputfiles[job]} ${outputdir}

