# parameters used in chi square calculation
e_photnPurity	       =  0.574992272157
e_photnPurityErr       =  0.0611339548487
e_M3_photon_topFrac    =  0.729435920088
e_M3_photon_topFracErr =  0.0610220962652
e_Ndata		       =  935.0
e_NdataErr             =  30.5777697028

# parameters used in chi square calculation
mu_photnPurity	        =  0.530947306374
mu_photnPurityErr       =  0.0559716612967
mu_M3_photon_topFrac    =  0.705992344212
mu_M3_photon_topFracErr =  0.0583684417455
mu_Ndata	        =  1136.0
mu_NdataErr             =  33.7045990927


import likelihoodCombination

likelihoodCombination.e_data  = {'photnPurity'          : e_photnPurity	        ,     
                                 'photnPurityErr'       : e_photnPurityErr      ,
                                 'M3_photon_topFrac'    : e_M3_photon_topFrac   ,
                                 'M3_photon_topFracErr' : e_M3_photon_topFracErr,
                                 'Ndata'                : e_Ndata	        ,
                                 'NdataErr'             : e_NdataErr            ,
                                 }

likelihoodCombination.mu_data  = {'photnPurity'          : mu_photnPurity         ,     
                                  'photnPurityErr'       : mu_photnPurityErr      ,
                                  'M3_photon_topFrac'    : mu_M3_photon_topFrac   ,
                                  'M3_photon_topFracErr' : mu_M3_photon_topFracErr,
                                  'Ndata'                : mu_Ndata               ,
                                  'NdataErr'             : mu_NdataErr            ,
                                  }

eFile = 'ratioFiles_ele/templates_barrel_scaled_afterPhotonM3_nominal.root'
muFile = 'ratioFiles_mu/templates_barrel_scaled_afterPhotonM3_nominal.root'

efficiencies = {'TTgammaPhoEffAccErr': 0.0004750033079833124, 'TTgammaPhoEffAccVisErr': 0.0010776668819454565, 'TTJets_topEffAcc': 0.078461458262999995, 'TTgammaPhoEffAccVis': 0.046971944490986239, 'topPreselInt': 378648.02083299996, 'TTgammaPhoEffAcc': 0.018382076117990417, 'topPreselErr': 2003.5308289752568, 'TTJets_topEffAccErr': 3.6002086032723661e-05}

print '-'*30
print 'ele_value'
print '-'*30

likelihoodCombination.mu_data['photnPurityErr'] = 9999999999.
likelihoodCombination.mu_data['M3_photon_topFracErr'] = 9999999999.
likelihoodCombination.mu_data['NdataErr'] = 9999999999.


result_ele = likelihoodCombination.calculateTTGamma(eFile, muFile, efficiencies,saveFitPlots = False, verbose = False, progressBar = True)

print
print result_ele


likelihoodCombination.e_data['photnPurityErr'] = 0.0000000001
result_photn_fixed_ele = likelihoodCombination.calculateTTGamma(eFile, muFile, efficiencies,saveFitPlots = False, verbose = False, progressBar = True)
likelihoodCombination.e_data['photnPurityErr'] = e_photnPurityErr
print
print 'Fix Photon Purity'
print
print result_photn_fixed_ele
print

likelihoodCombination.e_data['M3_photon_topFracErr'] = 0.0000000001
result_topFrac_fixed_ele = likelihoodCombination.calculateTTGamma(eFile, muFile, efficiencies,saveFitPlots = False, verbose = False, progressBar = True)
likelihoodCombination.e_data['M3_photon_topFracErr'] = e_M3_photon_topFracErr
print
print 'Fix Photon Purity'
print
print result_topFrac_fixed_ele
print

likelihoodCombination.e_data['NdataErr'] = 0.0000000001
result_Ndata_fixed_ele = likelihoodCombination.calculateTTGamma(eFile, muFile, efficiencies,saveFitPlots = False, verbose = False, progressBar = True)
likelihoodCombination.e_data['NdataErr'] = e_NdataErr
print
print 'Fix Photon Purity'
print
print result_Ndata_fixed_ele
print

print '-'*30
print 'mu_value'
print '-'*30

likelihoodCombination.e_data['photnPurityErr'] = 9999999999.
likelihoodCombination.e_data['M3_photon_topFracErr'] = 9999999999.
likelihoodCombination.e_data['NdataErr'] = 9999999999.

likelihoodCombination.mu_data['photnPurityErr'] = mu_photnPurityErr
likelihoodCombination.mu_data['M3_photon_topFracErr'] = mu_M3_photon_topFracErr
likelihoodCombination.mu_data['NdataErr'] = mu_NdataErr


result_mu = likelihoodCombination.calculateTTGamma(eFile, muFile, efficiencies,saveFitPlots = False, verbose = False, progressBar = True)

print
print result_mu


likelihoodCombination.mu_data['photnPurityErr'] = 0.0000000001
result_photn_fixed_mu = likelihoodCombination.calculateTTGamma(eFile, muFile, efficiencies,saveFitPlots = False, verbose = False, progressBar = True)
likelihoodCombination.mu_data['photnPurityErr'] = mu_photnPurityErr
print
print 'Fix Photon Purity'
print
print result_photn_fixed_mu
print

likelihoodCombination.mu_data['M3_photon_topFracErr'] = 0.0000000001
result_topFrac_fixed_mu = likelihoodCombination.calculateTTGamma(eFile, muFile, efficiencies,saveFitPlots = False, verbose = False, progressBar = True)
likelihoodCombination.mu_data['M3_photon_topFracErr'] = mu_M3_photon_topFracErr
print
print 'Fix Photon Purity'
print
print result_topFrac_fixed_mu
print

likelihoodCombination.mu_data['NdataErr'] = 0.0000000001
result_Ndata_fixed_mu = likelihoodCombination.calculateTTGamma(eFile, muFile, efficiencies,saveFitPlots = False, verbose = False, progressBar = True)
likelihoodCombination.mu_data['NdataErr'] = mu_NdataErr
print
print 'Fix Photon Purity'
print
print result_Ndata_fixed_mu
print

print
print
print
total_e = result_ele[2][0][1]
photn_e = result_photn_fixed_ele[2][0][1]
topfrac_e = result_topFrac_fixed_ele[2][0][1]
nData_e = result_Ndata_fixed_ele[2][0][1]
print 'e photon:', (total_e**2 - photn_e**2)**0.5
print 'e topFrac:', (total_e**2 - topfrac_e**2)**0.5
print 'e nData:', (total_e**2 - nData_e**2)**0.5

print
print
print
total_mu = result_mu[2][0][1]
photn_mu = result_photn_fixed_mu[2][0][1]
topfrac_mu = result_topFrac_fixed_mu[2][0][1]
nData_mu = result_Ndata_fixed_mu[2][0][1]
print 'e photon:', (total_mu**2 - photn_mu**2)**0.5
print 'e topFrac:', (total_mu**2 - topfrac_mu**2)**0.5
print 'e nData:', (total_mu**2 - nData_mu**2)**0.5
