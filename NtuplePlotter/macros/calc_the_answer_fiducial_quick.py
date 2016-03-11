import ROOT
import signal_finder_quick

def integral(file,histname):
	err = ROOT.Double(0.0)
	hist = file.Get(histname)
	integr = hist.IntegralAndError(0,hist.GetNbinsX() + 1, err)
#	print integr,err
	return integr,err

def integral_bins(file,histname,bin1,bin2):
	err = ROOT.Double(0.0)
	hist = file.Get(histname)
	integr = hist.IntegralAndError(bin1, bin2, err)
#	print integr,err
	return integr,err

def integral_one_bin(file,histname):
	err = ROOT.Double(0.0)
	hist = file.Get(histname)
	integr = hist.IntegralAndError(1, 1, err)
#	print integr,err
	return integr,err

# for ttbar acceptance
TTJets1l_num        = 24849110
TTJets2l_num        = 12086717
TTJetsHad_num       = 31178278

## >=3j
photnPurity = 0.657
photnPurityErr = 0.065

M3TopSF = 0.994
M3TopSFErr = 0.0075

M3WJetsSF = 2.527
M3WJetsSFErr = 0.044

eleFakeSF = 1.5 
eleFakeSFErr = 0.2

# experiment
M3_photon_topFrac = 0.619
M3_photon_topFracErr = 0.089

fidSelectionEff = -1
fidSelectionEffErr = -1

channel = "ele"
# preselFileName = 'templates_presel_scaled.root'
# barrelFileName = 'templates_barrel_scaled.root'
preselFileName = 'templates_presel_scaled_MCCat.root'
barrelFileName = 'templates_barrel_scaled_MCCat.root'
barrelFileName_M3fitscaled = 'templates_barrel_scaled_MCCat.root'

ttjetsSyst = ''

ttjetsNum ={}

def doTheCalculation(verbose = True, progressBar = False, saveDataFile = False):
	if verbose:
		print '#$'*50
		print 'input data'
		print 'photon Purity: ',photnPurity,' +-',photnPurityErr
		print 'TopSF: ',M3TopSF,' +-',M3TopSFErr
		print 'WJetsSF: ',M3WJetsSF,' +-',M3WJetsSFErr
		print 'eleFakeSF: ',eleFakeSF,' +-',eleFakeSFErr
		
		print
		print 'Opening file with pre-selection histograms', preselFileName
	preselfile = ROOT.TFile(preselFileName,'READ')

	if verbose:
		print 'File was last modified: ',preselfile.GetModificationDate().AsString()
		print

	# integraph of TTJets passing the preselection (ttbar plus ttgamma)

	#Get number of ttbar events passing the preselection
	topPreselInt,topPreselErr = integral(preselfile,'TTJets_MET')
	if verbose:
		print 'TTJets presel events expected ',topPreselInt,' +-',topPreselErr,'(MC stat)'
		print

	#Get number of ttgamma events passing the preselection
	whizPreselInt,whizPreselErr = integral(preselfile,'TTGamma_MET')

	if verbose:
		print 'TTGamma presel events expected ',whizPreselInt,' +-',whizPreselErr,'(MC stat)'
		print
		print
	# total number of events with top quark pair is the sum of them
	topPreselInt += whizPreselInt
	topPreselErr = (topPreselErr**2 + whizPreselErr**2 + (topPreselInt*M3TopSFErr/M3TopSF)**2 )**0.5

	if verbose:
		print 'taking into account fit error and ttgamma contribution: ',topPreselInt,' +-',topPreselErr,'(MC stat+fit)'
		print

	# number of ttgamma events in preselection, generated in the 1l or 2l final state
	ttgammaSigPreselInt, ttgammaSigPreselErr = integral_bins(preselfile,'TTGamma_MCcategory',3,4)

	if verbose:
		print 'number of signal events with one or two leptons in preselection ', ttgammaSigPreselInt, ' +-', ttgammaSigPreselErr,'(MC stat)'
	
		print 'Opening file with photon in barrel', barrelFileName

	barrelfile = ROOT.TFile(barrelFileName,'READ')
	barrelfile2 = ROOT.TFile(barrelFileName_M3fitscaled,'READ')
	if verbose:
		print 'File was last modified: ',barrelfile.GetModificationDate().AsString()
		print
	# Get number of ttgamma events, passing the photon selection, with a "signal" photon, generated in the 1l or 2l final state
	ttgammaBarrInt,ttgammaBarrErr = integral_bins(barrelfile,'TTGamma_signal_MCcategory',3,4)
	
	if verbose:
		print 'ttgamma signal events with one or two leptons and matched photon in barrel expected ',ttgammaBarrInt, ' +-',ttgammaBarrErr,'(MC stat)'
		print
	
	#photon acceptance is the number of ttbar events passing photon selection divided by number passing the preselection (all generated 1l or 2l final states)
	phoAcc = ttgammaBarrInt/ttgammaSigPreselInt
	# These errors are not independent
	phoAccErr = phoAcc * ((ttgammaSigPreselErr/ttgammaSigPreselInt) + (ttgammaBarrErr/ttgammaBarrInt))
	if verbose:
		print 'signal photon acceptance',phoAcc,' +-',phoAccErr
		print
	
		print 'reconstruction efficiency: matched signal events in barrel'

	ttgammaSigBarrInt,ttgammaSigBarrErr = integral(barrelfile,'TTGamma_signal_MET')
	
	if verbose: print 'fiducial photons in barrel after pre-selection'

	genPhoInAcc,genPhoInAccErr = integral_one_bin(preselfile,'TTGamma_genPhoRegionWeight')
	phoRecoEff = ttgammaSigBarrInt/genPhoInAcc
	phoRecoEffErr = phoRecoEff * ((genPhoInAccErr/genPhoInAcc)**2 + (ttgammaSigBarrErr/ttgammaSigBarrInt)**2)**0.5

	if verbose:
		print 'signal photon reconstruction efficiency (with truth matching to photon)', phoRecoEff,' +-',phoRecoEffErr
		print

	if verbose: print 'TTJets acceptance calculation'
	accfile = ROOT.TFile('ttbar_acceptance.root','READ')
	if verbose: print 'File was last modified: ',accfile.GetModificationDate().AsString()

	ttjets1l_top, ttjets1l_topErr = integral_one_bin(accfile,'TTJets1l_presel_MCcategory')
	ttjets2l_top, ttjets2l_topErr = integral_one_bin(accfile,'TTJets2l_presel_MCcategory')
	ttjetsHad_top, ttjetsHad_topErr = integral_one_bin(accfile,'TTJetsHad_presel_MCcategory')
	Whad = 0.676
	TTJets_topEffAcc = Whad*(1.0-Whad)*2 * ttjets1l_top / TTJets1l_num + (1.0-Whad)*(1.0-Whad) * ttjets2l_top / TTJets2l_num + Whad*Whad * ttjetsHad_top / TTJetsHad_num
	TTJets_topEffAccErr = ( (Whad*(1.0-Whad)*2 * ttjets1l_topErr / TTJets1l_num)**2 + ((1.0-Whad)*(1.0-Whad) * ttjets2l_topErr / TTJets2l_num)**2 + (Whad*Whad * ttjetsHad_topErr / TTJetsHad_num)**2 )**0.5

	if not ttjetsSyst=='':
		ttjets_top, ttjets_topErr = integral_one_bin(accfile,'TTJets_presel_MCcategory')
		TTJets_topEffAcc = ttjets_top/ttjetsNum[ttjetsSyst]
		TTJets_topEffAccErr = ttjets_topErr/ttjetsNum[ttjetsSyst]


	if verbose:
		print 'Inclusive ttbar acceptance ',TTJets_topEffAcc,' +-',TTJets_topEffAccErr
		print 
		print 'TTgamma acceptance calculation'

	sigaccfile = ROOT.TFile('signalAcc.root','READ')
	if verbose: print 'File was last modified: ',sigaccfile.GetModificationDate().AsString()

	ttgamma_1l_2l_num, ttgamma_1l_2l_numErr = integral_bins(sigaccfile,'allCategory',3,4)
	
	ttgamma_1l_2l_presel, ttgamma_1l_2l_preselErr = integral_bins(accfile, 'TTGamma_presel_MCcategory',3,4)
	ttgamma_1l_2l_sig, ttgamma_1l_2l_sigErr = integral_bins(accfile, 'TTGamma_signal_MCcategory',3,4)
	
	TTGamma_topEffAcc = ttgamma_1l_2l_presel / ttgamma_1l_2l_num
	TTGamma_topEffAccErr = TTGamma_topEffAcc * ((ttgamma_1l_2l_preselErr/ttgamma_1l_2l_presel)**2 + (ttgamma_1l_2l_numErr/ttgamma_1l_2l_num)**2)**0.5

	if verbose: print 'Acceptance for ttgamma after top selection ',TTGamma_topEffAcc, ' +/-',TTGamma_topEffAccErr

	ttgammaAccFile = ROOT.TFile('ttgamma_acceptance.root','READ')
	ttgammaAccFile2 = ROOT.TFile('ttgamma_acceptance2.root','READ')

	ttgamma_1l_2l_Vispresel, ttgamma_1l_2l_VispreselErr = integral_one_bin(ttgammaAccFile2,'TTGamma_presel_genPhoRegionWeight_1l_2l')
	ttgamma_1l_2l_Visnum, ttgamma_1l_2l_VisnumErr = integral_bins(sigaccfile,'VisAllCategory',3,4)
	
#	TTGammaVis_topAcc = ttgamma_1l_2l_presel / ttgamma_1l_2l_Visnum
	TTGammaVis_topAcc = ttgamma_1l_2l_Vispresel / ttgamma_1l_2l_Visnum

	TTGammaVis_topAccErr = TTGammaVis_topAcc * ((ttgamma_1l_2l_VispreselErr/ttgamma_1l_2l_Vispresel)**2 + (ttgamma_1l_2l_VisnumErr/ttgamma_1l_2l_Visnum)**2)**0.5

	if verbose:
		print 'Acceptance for visible ttgamma after top selection ',TTGammaVis_topAcc, ' +/-',TTGammaVis_topAccErr
		print
		print 'Product of ttgamma acceptance for top and photon calculated at once',ttgamma_1l_2l_sig/ttgamma_1l_2l_num
		print 
	
	
	DataBarrInt,DataBarrErr = integral(barrelfile2,'Data_MET')
	if verbose: print 'Data events',DataBarrInt,' +-',DataBarrErr
	





	
	###################
	### DEFINE FIDUCIAL CROSS SECTION VALUES
	###################

	#FidCX = Nttgamma/FidTopEff/FidPhoEff
	

	#FidTopEff is the percentage of events generated in the e+jets (or mu+jets) final state (including tau-> mu or tau-> ele) which pass the top selection

	#FidPhoEff is the percentage of events that pass the preselection, which have a generated photon in the photon phase space, which pass the cuts
	# This shouldn't really depend on the generated state (in principle), but we still consider only events in the e+jets or mu+jets final state for this


	if channel=='ele':
		categoryBin = 14
	if channel=='mu':
		categoryBin = 19


	#Get the number of events in the ttgamma MC with e or mu+jets final state in fiducial area (generated)
	#Get number of ttgamma events generated with ttgamma acceptance cuts applied (1 lepton , topsel, phosel)
	nGeneratedSemiLepChannel, nGeneratedSemiLepChannelErr = integral_bins(sigaccfile,'allCategory',categoryBin, categoryBin)

	#Get number of ttgamma events after preselection with ttgamma acceptance cuts applied, SF applied, but not scaled to lumi
	nPreselSemiLepChannel, nPreselSemiLepChannelErr = integral_bins(ttgammaAccFile,'TTGamma_presel_MCcategoryfid',categoryBin, categoryBin)
	
	nPhotonSelSemiLepChannel, nPhotonSelSemiLepChannelErr = integral_bins(ttgammaAccFile,'TTGamma_photonSel_MCcategoryfid',categoryBin, categoryBin)


	if verbose: print nPreselSemiLepChannelErr,nPreselSemiLepChannel
	if verbose: print nGeneratedSemiLepChannelErr,nGeneratedSemiLepChannel

	FidTopEff = nPreselSemiLepChannel/nGeneratedSemiLepChannel
	FidTopEffErr = ((nPreselSemiLepChannelErr/nPreselSemiLepChannel)**2 + (nGeneratedSemiLepChannelErr/nGeneratedSemiLepChannel)**2)**0.5

	FidEff = nPhotonSelSemiLepChannel/nGeneratedSemiLepChannel
	FidEffErr = FidEff*((nPhotonSelSemiLepChannelErr/nPhotonSelSemiLepChannel)**2 + (nGeneratedSemiLepChannelErr/nGeneratedSemiLepChannel)**2)**0.5

	print 'initial eff of ttgamma fid:', FidEff, FidEffErr
	print nPhotonSelSemiLepChannel
	print nGeneratedSemiLepChannel

	if not fidSelectionEff == -1:
		FidEff = fidSelectionEff
		FidEffErr = fidSelectionEffErr


	#number of events passing the preselection which pass the top and photon cuts at the generator level, scaled to lumi
	nPreselSemiLepGenPhoton, nPreselSemiLepGenPhotonErr = integral_bins(preselfile,"TTGamma_MCcategoryfid", categoryBin, categoryBin)
	#number of events passing the photon selection which pass the top and photon cuts at the generator level, scaled to lumi
	nPhotonSelSemiLepGenPhoton, nPhotonSelSemiLepGenPhotonErr = integral_bins(barrelfile,"TTGamma_MCcategoryfid", categoryBin, categoryBin)

	FidPhoEff = nPhotonSelSemiLepGenPhoton/nPreselSemiLepGenPhoton
	FidPhoEffErr = ((nPhotonSelSemiLepGenPhotonErr/nPhotonSelSemiLepGenPhoton)**2 + (nPreselSemiLepGenPhotonErr/nPreselSemiLepGenPhoton)**2)**0.5

	nGenerated, nGeneratedErr = integral_bins(sigaccfile,'allCategory',1,1)

	FidAcc = nGeneratedSemiLepChannel/nGenerated
	FidAccErr = ((nGeneratedSemiLepChannelErr/nGeneratedSemiLepChannel)**2 + (nGeneratedErr/nGenerated)**2)**0.5
	
	# here we need a procedure to find NttgammaSignal
	signal_finder_quick.barrelFileName = barrelFileName_M3fitscaled



	signal_finder_quick.eleFakeSF = eleFakeSF
	signal_finder_quick.eleFakeSFErr = eleFakeSFErr
	signal_finder_quick.photnPurity = photnPurity
	signal_finder_quick.photnPurityErr = photnPurityErr
	signal_finder_quick.M3_photon_topFrac = M3_photon_topFrac
	signal_finder_quick.M3_photon_topFracErr = M3_photon_topFracErr
	signal_finder_quick.Ndata = DataBarrInt
	signal_finder_quick.NdataErr = DataBarrErr


	if verbose:
		print '*'*80
		print 'eleFakeSF	    = ', signal_finder_quick.eleFakeSF		   
		print 'eleFakeSFErr	    = ', signal_finder_quick.eleFakeSFErr	   
		print 'photnPurity	    = ', signal_finder_quick.photnPurity	   
		print 'photnPurityErr	    = ', signal_finder_quick.photnPurityErr	   
		print 'M3_photon_topFrac    = ', signal_finder_quick.M3_photon_topFrac	   
		print 'M3_photon_topFracErr = ', signal_finder_quick.M3_photon_topFracErr 
		print 'Ndata		    = ', signal_finder_quick.Ndata		   
		print 'NdataErr             = ', signal_finder_quick.NdataErr             
		
	ttgammaSig,ttgammaSigErr, ttgammaFull,ttgammaFullErr, bestttgSF, bestttgSFErr, bestvgammaSF, bestvgammaSFErr, bestjgSF, bestjgSFErr = signal_finder_quick.calculateTTGamma(progressBar,saveDataFile)





	#####################
	### Calculate cross section directly (not ratio)
	#####################

	from SF import luminosity
	xsDirect = ttgammaSig / phoAcc / TTGamma_topEffAcc / luminosity
	xsDirectErr =  xsDirect * ( (ttgammaSigErr/ttgammaSig)**2 + 
				    (phoAccErr/phoAcc)**2 + 
				    (TTGamma_topEffAccErr/TTGamma_topEffAcc)**2 +
				    (0.022)**2
				    )**0.5
	
	if verbose:
		print 
		print '*'*80
		print 'Direct Cross Section Value: %0.5f  +- %0.5f' %( xsDirect, xsDirectErr)
		print '*'*80
		print

#	xsDirect = ttgammaFull / FidPhoEff / FidTopEff / luminosity
	# xsDirectErr =  xsDirect * ( (ttgammaSigErr/ttgammaSig)**2 + 
	# 			    (phoAccErr/phoAcc)**2 + 
	# 			    (TTGamma_topEffAccErr/TTGamma_topEffAcc)**2 +
	# 			    (0.022)**2
	# 			    )**0.5

	xsDirect = ttgammaFull / FidEff / luminosity
	xsDirectErr =  xsDirect * ( (ttgammaFullErr/ttgammaFull)**2 + 
				    (FidEffErr/FidEff)**2 + 
				    (0.022)**2
				    )**0.5


	if verbose:
		print 
		print '*'*80
		print 'Direct Fiducial Cross Section Value: %0.5f  +- %0.5f' %( xsDirect, xsDirectErr)
		print '*'*80
		print

	#################################################
	### Fiducial cross section ratio
	### Ratio of fiducial ttgamma cross section to full ttbar cross section
	#################################################

	# fidxsRatio = ttgammaFull / FidPhoEff / FidTopEff / topPreselInt * TTJets_topEffAcc
	# fidxsRatioRelErr = ( (ttgammaFullErr/ttgammaFull)**2 + 
	# 		     (FidPhoEffErr)**2 + 
	# 		     (FidTopEffErr)**2 + 
	# 		     (topPreselErr/topPreselInt)**2 +
	# 		     (TTJets_topEffAccErr/TTJets_topEffAcc)**2
	# 		     )**0.5

	fidxsRatio = ttgammaFull / FidEff / topPreselInt * TTJets_topEffAcc
	fidxsRatioErr = fidxsRatio * ( (ttgammaFullErr/ttgammaFull)**2 + 
				       (FidEffErr/FidEff)**2 + 
				       (topPreselErr/topPreselInt)**2 +
				       (TTJets_topEffAccErr/TTJets_topEffAcc)**2
				       )**0.5
	
	if verbose:
		print 
		print '*'*80
		print 'Fiducial Cross Section Ratio Value: %0.7f  +- %0.7f' %( fidxsRatio, fidxsRatioErr)
		print '*'*80
		print

	
	#################################################
	
	xsRatio = ttgammaSig / phoAcc / TTGamma_topEffAcc / topPreselInt * TTJets_topEffAcc
	xsRatioRelErr = (   (ttgammaSigErr/ttgammaSig)**2 + 
						(phoAccErr/phoAcc)**2 + 
						(TTGamma_topEffAccErr/TTGamma_topEffAcc)**2 +
						(topPreselErr/topPreselInt)**2 +
						(TTJets_topEffAccErr/TTJets_topEffAcc)**2
					)**0.5

	if verbose:
		print '*'*80
		print 'final answer: cross section ratio:'
		print xsRatio,' +-',xsRatio*xsRatioRelErr
		print '*'*80
		
	
	vis_xsRatio = ttgammaSig / phoRecoEff / TTGammaVis_topAcc / topPreselInt * TTJets_topEffAcc
	vis_xsRatioErr = (  (ttgammaSigErr/ttgammaSig)**2 + 
						(phoRecoEffErr/phoRecoEff)**2 +
						(TTGammaVis_topAccErr/TTGammaVis_topAcc)**2 +
						(topPreselErr/topPreselInt)**2 +
						(TTJets_topEffAccErr/TTJets_topEffAcc)**2
					)**0.5
	
	if verbose:
		print 'visible cross section ratio:'
		print vis_xsRatio,' +-',vis_xsRatio*vis_xsRatioErr
		print '*'*80

		print '*'*80
		print "ttgammaSig        = ", ttgammaSig  
		print "ttgammaSigErr        = ", ttgammaSigErr  
		print
		print "ttgammaFull        = ", ttgammaFull  
		print "ttgammaFullErr        = ", ttgammaFullErr  
		print
		print "FidEff            = ", FidEff
		print "FidEffErr         = ", FidEffErr
		print "FidTopEff         = ", FidTopEff
		print "FidTopEffErr      = ", FidTopEffErr
		print "FidPhoEff         = ", FidPhoEff
		print "FidPhoEffErr      = ", FidPhoEffErr
		print "FidAcc            = ", FidAcc
		print "FidAccErr         = ", FidAccErr
		print
		print "nPreselSemiLepChannel      = ", nPreselSemiLepChannel
		print "nGeneratedSemiLepChannel   = ", nGeneratedSemiLepChannel
		print "nPhotonSelSemiLepGenPhoton = ", nPhotonSelSemiLepGenPhoton
		print "nPreselSemiLepGenPhoton    = ", nPreselSemiLepGenPhoton
		print "nGenerated                 = ", nGenerated
		print "nPreselSemiLepChannelErr      = ", nPreselSemiLepChannelErr
		print "nGeneratedSemiLepChannelErr   = ", nGeneratedSemiLepChannelErr
		print "nPhotonSelSemiLepGenPhotonErr = ", nPhotonSelSemiLepGenPhotonErr
		print "nPreselSemiLepGenPhotonErr    = ", nPreselSemiLepGenPhotonErr
		print "nGeneratedErr                 = ", nGeneratedErr
		print
		

		print "phoAcc            = ", phoAcc 	       
		print "TTGamma_topEffAcc = ", TTGamma_topEffAcc 
		print "topPreselInt      = ", topPreselInt     
		print "TTJets_topEffAcc  = ", TTJets_topEffAcc 
		print 'phoRecoEff        = ', phoRecoEff
		print 'TTGammaVis_topAcc = ', TTGammaVis_topAcc
		print
		print "phoAccErr            = ", phoAccErr  
		print "TTGamma_topEffAccErr = ", TTGamma_topEffAccErr
		print "topPreselErr         = ", topPreselErr    
		print "TTJets_topEffAccErr  = ", TTJets_topEffAccErr
		print 'phoRecoEffErr        = ', phoRecoEffErr
		print 'TTGammaVis_topAccErr = ', TTGammaVis_topAccErr

	return xsRatio, xsRatio*xsRatioRelErr,bestttgSF, bestttgSFErr, bestvgammaSF, bestvgammaSFErr, bestjgSF, bestjgSFErr
