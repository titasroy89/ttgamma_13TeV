#!/bin/bash

set hist="hist"
set dir="root://cmseos.fnal.gov//store/user/troy2012/skims_13TeV/"

NtuplePlotter/makeTemplates data_mu_b $hist ${dir}data_mu_b.root
NtuplePlotter/makeTemplates data_mu_c $hist ${dir}data_mu_c.root
NtuplePlotter/makeTemplates data_mu_d $hist ${dir}data_mu_d.root
NtuplePlotter/makeTemplates DYJets $hist ${dir}DYJets.root
NtuplePlotter/makeTemplates ST_s $hist ${dir}ST_s.root
NtuplePlotter/makeTemplates ST_t_bar $hist ${dir}ST_t_bar.root
NtuplePlotter/makeTemplates ST_t_top $hist ${dir}ST_t_top.root
NtuplePlotter/makeTemplates ST_tW_antitop $hist ${dir}ST_tW_antitop.root
NtuplePlotter/makeTemplates ttgamma $hist ${dir}ttgamma.root
NtuplePlotter/makeTemplates Wjets $hist ${dir}Wjets.root
NtuplePlotter/makeTemplates TTbar $hist ${dir}TTbar_1.root

