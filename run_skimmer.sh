#!/bin/bash

job=$1
#job=5

cd /uscms/home/troy2012/ttgamma_13TeV/CMSSW_7_4_14/src/
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
cd /uscms/home/troy2012/ttgamma_13TeV/CMSSW_7_4_14/src/TTGammaSemiLep/

outputDir="/uscms/home/troy2012/ttgamma_13TeV/CMSSW_7_4_14/src/TTGammaSemiLep/"
inputfiles=("sync.root /uscms/home/troy2012/ttgamma_13TeV/CMSSW_7_4_14/src/ggAnalysis/ggNtuplizer/test/ggtree_sync.root")

echo "NtuplePlotter/makeSkim ${outputDir}${inputfiles[job]}"

NtuplePlotter/makeSkim ${outputDir}${inputfiles[job]}
