#!/bin/csh

set hist="QCD_bin"
set dir="/eos/uscms/store/user/makouski/PreSelSkimQCD/"

# set hist="histTest"
# set dir="/eos/uscms/store/user/makouski/PreSelSkim/"

NtuplePlotter/makeTemplates Data_a        $hist ${dir}1electron_a_presel.root
NtuplePlotter/makeTemplates Data_b        $hist ${dir}1electron_b_presel.root
NtuplePlotter/makeTemplates Data_c        $hist ${dir}1electron_c_presel.root
NtuplePlotter/makeTemplates Data_d        $hist ${dir}1electron_d_presel.root

NtuplePlotter/makeTemplates W4Jets        $hist ${dir}job_summer12_W4jets_presel.root
NtuplePlotter/makeTemplates W3Jets        $hist ${dir}job_summer12_W3jets_presel.root
NtuplePlotter/makeTemplates TTJets2l      $hist ${dir}job_summer12_ttjets_2l_presel.root 
NtuplePlotter/makeTemplates TTJets1l      $hist ${dir}job_summer12_ttjets_1l_presel.root 
NtuplePlotter/makeTemplates TTJetsHad     $hist ${dir}job_summer12_ttjets_had_presel.root
NtuplePlotter/makeTemplates TTGamma       $hist ${dir}job_summer12_ttg_presel.root
NtuplePlotter/makeTemplates ZJets         $hist ${dir}job_summer12_Zjets_presel.root          
NtuplePlotter/makeTemplates Zgamma        $hist ${dir}job_summer12_Zg_presel.root	      
NtuplePlotter/makeTemplates Wgamma        $hist ${dir}job_summer12_Wg_presel.root	      
NtuplePlotter/makeTemplates SingleT_t     $hist ${dir}job_summer12_t_t_presel.root
NtuplePlotter/makeTemplates SingleT_s     $hist ${dir}job_summer12_t_s_presel.root
NtuplePlotter/makeTemplates SingleT_tw    $hist ${dir}job_summer12_t_tW_presel.root
NtuplePlotter/makeTemplates SingleTbar_t  $hist ${dir}job_summer12_tbar_t_presel.root
NtuplePlotter/makeTemplates SingleTbar_s  $hist ${dir}job_summer12_tbar_s_presel.root
NtuplePlotter/makeTemplates SingleTbar_tw $hist ${dir}job_summer12_tbar_tW_presel.root


NtuplePlotter/makeTemplates WW_2l2nu      $hist ${dir}job_summer12_WW_2l2nu_presel.root
NtuplePlotter/makeTemplates WZ_2l2q       $hist ${dir}job_summer12_WZ_2l2q_presel.root
NtuplePlotter/makeTemplates WZ_3lnu       $hist ${dir}job_summer12_WZ_3lnu_presel.root
NtuplePlotter/makeTemplates ZZ_2e2mu      $hist ${dir}job_summer12_ZZ_2e2mu_presel.root
NtuplePlotter/makeTemplates ZZ_2e2tau     $hist ${dir}job_summer12_ZZ_2e2tau_presel.root
NtuplePlotter/makeTemplates ZZ_2mu2tau    $hist ${dir}job_summer12_ZZ_2mu2tau_presel.root
NtuplePlotter/makeTemplates ZZ_4e         $hist ${dir}job_summer12_ZZ_4e_presel.root
NtuplePlotter/makeTemplates ZZ_4mu        $hist ${dir}job_summer12_ZZ_4mu_presel.root
NtuplePlotter/makeTemplates ZZ_4tau       $hist ${dir}job_summer12_ZZ_4tau_presel.root
