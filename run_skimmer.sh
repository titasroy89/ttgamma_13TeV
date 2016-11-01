#!/bin/bash

job=$1

workingnode=$(pwd)

echo $workingnode

tar -xzf Summer15_25nsV6_MC.tar.gz

cd /uscms_data/d3/troy2012/ttgamma_13TeV/CMSSW_7_4_14/src/TTGammaSemiLep/
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

cd $workingnode
outputdir="root://cmseos.fnal.gov//store/user/troy2012/skims_13TeV/"
outputfiles=("data_d_PR.root" \
"st_t.root")
inputfiles=("root://cmseos.fnal.gov//store/user/troy2012/ntuples_13TeV/job_SingleMu_Run2015D_PR_v4_miniAOD.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_13TeV/job_spring15_ST_t_top_4f_miniAOD.root")
echo "./makeSkim ${outputfiles[job]} ${inputfiles[job]}"
./makeSkim ${outputfiles[job]} ${inputfiles[job]}

echo " xrdcp -f ${outputfiles[job]} ${outputdir} "
xrdcp -f ${outputfiles[job]} ${outputdir}

