#!/bin/csh

# set hist="QCD_bin"
# set dir="/eos/uscms/store/user/makouski/PreSelSkimQCD/"

set hist="/eos/uscms/store/user/dnoonan/MuHists_looseVeto/QCD_bins/"
set dir="/eos/uscms/store/user/dnoonan/QCDskims_mu/"


mkdir ${hist}

# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates Data_a $hist ${dir}skim_data_mu_a.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates Data_b $hist ${dir}skim_data_mu_b.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates Data_c $hist ${dir}skim_data_mu_c.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates Data_d $hist ${dir}skim_data_mu_d.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates W4Jets $hist ${dir}skim_W4jets.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates W3Jets $hist ${dir}skim_W3jets.root
# #/uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/#NtuplePlotter/makeTemplates W2Jets $hist ${dir}skim_W2jets.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates TTJets2l $hist ${dir}skim_ttjets_2l.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates TTJets1l $hist ${dir}skim_ttjets_1l.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates TTJetsHad $hist ${dir}skim_ttjets_0l.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates TTGamma $hist ${dir}skim_ttg.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates ZJets $hist ${dir}skim_DYJetsToLL.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates Zgamma $hist ${dir}skim_Zg.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates Wgamma $hist ${dir}skim_Wg.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates SingleT_t $hist ${dir}skim_t_t.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates SingleT_s $hist ${dir}skim_t_s.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates SingleT_tw $hist ${dir}skim_t_tW.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates SingleTbar_t $hist ${dir}skim_tbar_t.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates SingleTbar_s $hist ${dir}skim_tbar_s.root
# /uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/makeTemplates SingleTbar_tw $hist ${dir}skim_tbar_tW.root


NtuplePlotter/makeTemplates Data_a $hist ${dir}skim_data_mu_a.root
NtuplePlotter/makeTemplates Data_b $hist ${dir}skim_data_mu_b.root
NtuplePlotter/makeTemplates Data_c $hist ${dir}skim_data_mu_c.root
NtuplePlotter/makeTemplates Data_d $hist ${dir}skim_data_mu_d.root
NtuplePlotter/makeTemplates W4Jets $hist ${dir}skim_W4jets.root
NtuplePlotter/makeTemplates W3Jets $hist ${dir}skim_W3jets.root
NtuplePlotter/makeTemplates TTJets2l $hist ${dir}skim_ttjets_2l.root
NtuplePlotter/makeTemplates TTJets1l $hist ${dir}skim_ttjets_1l.root
NtuplePlotter/makeTemplates TTJetsHad $hist ${dir}skim_ttjets_0l.root
NtuplePlotter/makeTemplates TTGamma $hist ${dir}skim_ttg.root
NtuplePlotter/makeTemplates ZJets $hist ${dir}skim_DYJetsToLL.root
NtuplePlotter/makeTemplates Zgamma $hist ${dir}skim_Zg.root
NtuplePlotter/makeTemplates Wgamma $hist ${dir}skim_Wg.root
NtuplePlotter/makeTemplates SingleT_t $hist ${dir}skim_t_t.root
NtuplePlotter/makeTemplates SingleT_s $hist ${dir}skim_t_s.root
NtuplePlotter/makeTemplates SingleT_tw $hist ${dir}skim_t_tW.root
NtuplePlotter/makeTemplates SingleTbar_t $hist ${dir}skim_tbar_t.root
NtuplePlotter/makeTemplates SingleTbar_s $hist ${dir}skim_tbar_s.root
NtuplePlotter/makeTemplates SingleTbar_tw $hist ${dir}skim_tbar_tW.root
