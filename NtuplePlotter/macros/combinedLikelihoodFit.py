import ROOT
from math import exp

import sys
e_directory = "ElePlots/"
mu_directory = "MuPlots/"

preselFileName = 'templates_presel_scaled.root'
barrelFileName = 'templates_barrel_scaled.root'
barrelFileName_M3fitscaled = 'templates_barrel_scaled_afterPhotonM3.root'


myfile = None

def integral(histname, fileName):
	global myfile
	# if myfile == None:
	# 	myfile = ROOT.TFile(fileName,'READ')		
	myfile = ROOT.TFile(fileName,'READ')
	err = ROOT.Double(0.0)
	hist = myfile.Get(histname)
	integr = hist.IntegralAndError(0,hist.GetNbinsX() + 1, err)
	print integr,err
	return integr,err

eleFakeSF = 1.2 
eleFakeSFErr = 0.15

# parameters used in chi square calculation
e_photnPurity	       =  0.564158170272
e_photnPurityErr       =  0.0646043755
e_M3_photon_topFrac    =  0.683315978674
e_M3_photon_topFracErr =  0.0701062046615
e_Ndata		       =  977.0
e_NdataErr             =  31.2569992162

# parameters used in chi square calculation
mu_photnPurity	        =  0.531874490818
mu_photnPurityErr       =  0.0575003678412
mu_M3_photon_topFrac    =  0.69654529757
mu_M3_photon_topFracErr =  0.0610381089575
mu_Ndata	        =  1173.0
mu_NdataErr             =  34.2490875791

def readSamples(suffix, fileName):
	var = 'MET'
	######## Added QCD to the list of samples ########
	samplnames = ['TTGamma', 'TTJets', 'Vgamma', 'SingleTop', 'WJets', 'ZJets', 'QCD']
	samples = {}
	for n in samplnames:
		if n=='QCD' and (suffix=='signal' or suffix=='electron'):
			samples[n] = (0.0,0.0)
			continue
		print 'getting counts for ',n,suffix###, 'in', fileName
		if n=='QCD':
			int,err = integral(n+'_'+var, fileName)
		else:
			int,err = integral(n+'_'+suffix+'_'+var, fileName)
		samples[n] = (int,err)
	return samples

def Negamma(pho,ele,fake, ttgSF, VgSF, jgSF):
	sum = 0.0
	err2 = 0.0
	for sampln in pho:
		npho,nphoErr = pho[sampln]
		if sampln == 'TTGamma':
			npho *= ttgSF
			nphoErr *= ttgSF
		if sampln == 'Vgamma':
			npho *= VgSF
			nphoErr *= VgSF
		sum += npho
		err2 += nphoErr**2
	
	for sampln in ele:
		nele,neleErr = ele[sampln]
		nele *= eleFakeSF
		neleErr *= eleFakeSF
		if sampln == 'TTGamma':
			nele *= ttgSF
			neleErr *= ttgSF
		if sampln == 'Vgamma':
			nele *= VgSF
			neleErr *= VgSF
		sum += nele
		err2 += neleErr**2
	
	return sum,err2**0.5

def NallMC(pho,ele,fake, ttgSF, VgSF, jgSF):
	# all MC is egamma plus fakes
	sum,err = Negamma(pho,ele,fake, ttgSF, VgSF, jgSF)
	err2 = err**2
	for sampln in fake:
		nfake,nfakeErr = fake[sampln]
		nfake *= jgSF
		nfakeErr *= jgSF
		if sampln == 'TTGamma':
			nfake *= ttgSF
			nfakeErr *= ttgSF
		if sampln == 'Vgamma':
			nfake *= VgSF
			nfakeErr *= VgSF
		sum += nfake
		err2 += nfakeErr**2
	
	return sum,err2**0.5


def Ntop(pho,ele,fake, ttgSF, VgSF, jgSF):
	toppho = {}
	topele = {}
	topfake = {}
	for t in ['TTGamma', 'TTJets']:
		toppho[t] = pho[t]
		topele[t] = ele[t]
		topfake[t] = fake[t]
	
	return NallMC(toppho,topele,topfake, ttgSF, VgSF, jgSF)



def getChiSq_e(pho,ele,fake, ttgSF, VgSF, jgSF):
	allmc,allmcErr = NallMC(pho,ele,fake, ttgSF, VgSF, jgSF)
	top,topErr = Ntop(pho,ele,fake, ttgSF, VgSF, jgSF)
	egamma,egammaErr = Negamma(pho,ele,fake, ttgSF, VgSF, jgSF)
	cs = 0.0
	cs += (allmc - e_Ndata)**2 / (e_NdataErr**2 + allmcErr**2)
	
	mcPurErr = (egamma/allmc) * ((egammaErr/egamma)**2 + (allmcErr/allmc)**2)**0.5 
	cs += (egamma/allmc - e_photnPurity)**2 / (e_photnPurityErr**2 + mcPurErr**2)
	
	mcTopErr = (top/allmc) * ((topErr/top)**2 + (allmcErr/allmc)**2)**0.5
	cs += (top/allmc - e_M3_photon_topFrac)**2 / (e_M3_photon_topFracErr**2 + mcTopErr**2)
	return cs

def getChiSq_mu(pho,ele,fake, ttgSF, VgSF, jgSF):
	allmc,allmcErr = NallMC(pho,ele,fake, ttgSF, VgSF, jgSF)
	top,topErr = Ntop(pho,ele,fake, ttgSF, VgSF, jgSF)
	egamma,egammaErr = Negamma(pho,ele,fake, ttgSF, VgSF, jgSF)
	cs = 0.0
	cs += (allmc - mu_Ndata)**2 / (mu_NdataErr**2 + allmcErr**2)
	
	mcPurErr = (egamma/allmc) * ((egammaErr/egamma)**2 + (allmcErr/allmc)**2)**0.5 
	cs += (egamma/allmc - mu_photnPurity)**2 / (mu_photnPurityErr**2 + mcPurErr**2)
	
	mcTopErr = (top/allmc) * ((topErr/top)**2 + (allmcErr/allmc)**2)**0.5
	cs += (top/allmc - mu_M3_photon_topFrac)**2 / (mu_M3_photon_topFracErr**2 + mcTopErr**2)
	return cs


def getLikelihood(e_pho,e_ele,e_fake, mu_pho, mu_ele,mu_fake, ttgSF, VgSF, jgSF):
	return exp(-0.5*(getChiSq_e(e_pho,e_ele,e_fake, ttgSF, VgSF, jgSF)+getChiSq_mu(mu_pho,mu_ele,mu_fake, ttgSF, VgSF, jgSF)))

def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])
    else:
        return([])


def calculateTTGamma():
	# read all the numbers from the file
	e_pho = readSamples('signal', e_directory+barrelFileName_M3fitscaled)
	e_ele = readSamples('electron', e_directory+barrelFileName_M3fitscaled)
	e_fake = readSamples('fake', e_directory+barrelFileName_M3fitscaled)

	mu_pho = readSamples('signal', mu_directory+barrelFileName_M3fitscaled)
	mu_ele = readSamples('electron', mu_directory+barrelFileName_M3fitscaled)
	mu_fake = readSamples('fake', mu_directory+barrelFileName_M3fitscaled)
	
	ttghist = ROOT.TH1F('ttghist','Marginalized Likelihood', 141, 0.2-0.005, 1.6+0.005)
	vghist = ROOT.TH1F('vghist','Marginalized Likelihood', 61, 0-0.025, 3.0+0.025)
	jghist = ROOT.TH1F('jghist','Marginalized Likelihood', 131, 0.5-0.005, 1.8+0.005)

	maxlk = -1.0
	bestttgSF = -1.0
	bestVgSF = -1.0
	bestjgSF = -1.0
	
	percent = 0
	i = 0
	for ttgSF in seq(0.2, 1.6, 0.01):		
		if (i%7==0):
			sys.stdout.write("\r[" + "=" * (i / 7) + " " * (20- (i/7)) + "]" + str(percent) + "%")
			sys.stdout.flush()
			percent += 5
		i += 1
		for VgSF in seq(0.0, 3.0, 0.05):
			for jgSF in seq(0.5, 1.8, 0.01):
				lk = getLikelihood(e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake,ttgSF, VgSF, jgSF)
				ttghist.Fill(ttgSF,lk*0.0005)
				vghist.Fill(VgSF,lk*0.0001)
				jghist.Fill(jgSF,lk*0.0005)
				if lk > maxlk:
					maxlk = lk
					bestttgSF = ttgSF
					bestVgSF = VgSF
					bestjgSF = jgSF
	
	print 'max likelihood ',maxlk
	print 'best ttgSF',bestttgSF
	print 'best VgSF',bestVgSF
	print 'best jgSF',bestjgSF
	b_e_allmc,e_allmcErr = NallMC(e_pho,e_ele,e_fake, bestttgSF, bestVgSF, bestjgSF)
	b_e_top,e_topErr = Ntop(e_pho,e_ele,e_fake, bestttgSF, bestVgSF, bestjgSF)
	b_e_egamma,e_egammaErr = Negamma(e_pho,e_ele,e_fake, bestttgSF, bestVgSF, bestjgSF)

	b_mu_allmc,mu_allmcErr = NallMC(mu_pho,mu_ele,mu_fake, bestttgSF, bestVgSF, bestjgSF)
	b_mu_top,mu_topErr = Ntop(mu_pho,mu_ele,mu_fake, bestttgSF, bestVgSF, bestjgSF)
	b_mu_egamma,mu_egammaErr = Negamma(mu_pho,mu_ele,mu_fake, bestttgSF, bestVgSF, bestjgSF)

	print 'total number of electron MC events ',b_e_allmc, '  Data:',e_Ndata
	print 'photon purity in electron MC ',b_e_egamma/b_e_allmc, '  template fit:',e_photnPurity
	print 'top fraction in electron MC ',b_e_top/b_e_allmc, '  M3 fit:',e_M3_photon_topFrac

	print 'total number of muon MC events ',b_mu_allmc, '  Data:',mu_Ndata
	print 'photon purity in muon MC ',b_mu_egamma/b_mu_allmc, '  template fit:',mu_photnPurity
	print 'top fraction in muon MC ',b_mu_top/b_mu_allmc, '  M3 fit:',mu_M3_photon_topFrac
	
	ROOT.gStyle.SetOptFit(111)
	ccc = ROOT.TCanvas('ccc','ccc',800,800)
	ttghist.Draw()
	ttghist.GetXaxis().SetTitle('TTGamma Scale Factor')
	ttghist.Fit('gaus')
	ccc.SaveAs('TTGamma_SF_Lkhood.png')
	fit = ttghist.GetFunction('gaus')
	bestttgSFErr = fit.GetParameter(2)
	
	vghist.Draw()
	vghist.GetXaxis().SetTitle('Vgamma Scale Factor')
	vghist.SetMinimum(0.0)
	vghist.Fit('gaus')
	ccc.SaveAs('Vgamma_SF_Lkhood.png')
	
	jghist.Draw()
	jghist.GetXaxis().SetTitle('Jet to Photon Scale Factor')
	jghist.Fit('gaus')
	ccc.SaveAs('jet_gamma_SF_Lkhood.png')
	
	ttgammaSig = bestttgSF*(e_pho['TTGamma'][0]+mu_pho['TTGamma'][0])
	ttgammaSigErr = bestttgSFErr*(e_pho['TTGamma'][0]+mu_pho['TTGamma'][0])
	print 'number of signal events', ttgammaSig, ' +/-',ttgammaSigErr
	return ttgammaSig,ttgammaSigErr

combined_ttgammaSig, combined_ttgammaSigErr = calculateTTGamma()


#combined_ttgammaSig, combined_ttgammaSigErr = 762.213858252, 118.341451818
	
mu_phoAcc 	     =  0.13477462661
mu_TTGamma_topEffAcc =  0.079146136607
mu_topPreselInt      =  222697.219876
mu_TTJets_topEffAcc  =  0.0450047771279
mu_phoRecoEff        =  0.271603162212
mu_TTGammaVis_topAcc =  0.202572797378

mu_phoAccErr            =  0.00218140809392
mu_TTGamma_topEffAccErr =  0.000357352741188
mu_topPreselErr         =  1604.07716655
mu_TTJets_topEffAccErr  =  2.71182285127e-05
mu_phoRecoEffErr        =  0.00361732138042
mu_TTGammaVis_topAccErr =  0.000965526239153


e_phoAcc 	    =  0.131769190217
e_TTGamma_topEffAcc =  0.0635214446893
e_topPreselInt      =  160093.700896
e_TTJets_topEffAcc  =  0.0340351862929
e_phoRecoEff        =  0.278322649164
e_TTGammaVis_topAcc =  0.162581741774

e_phoAccErr            =  0.00240040666203
e_TTGamma_topEffAccErr =  0.000317815765302
e_topPreselErr         =  1365.70588597
e_TTJets_topEffAccErr  =  2.36146892948e-05
e_phoRecoEffErr        =  0.00419737860878
e_TTGammaVis_topAccErr =  0.00085048266326

###Combine the efficiencies between the two channels
combined_topPreselInt = mu_topPreselInt + e_topPreselInt
combined_TTgammaPhoEffAcc = mu_TTGamma_topEffAcc*mu_phoAcc + e_TTGamma_topEffAcc*e_phoAcc
combined_TTJets_topEffAcc = mu_TTJets_topEffAcc + e_TTJets_topEffAcc

combined_TTgammaPhoEffAccVis = mu_TTGamma_topEffAcc*mu_phoRecoEff + e_TTGammaVis_topAcc*e_phoRecoEff

###Calculate the combined uncertainties
combined_topPreselErr = (mu_topPreselErr**2+e_topPreselErr**2)**0.5
combined_TTgammaPhoEffAccErr = ((mu_TTGamma_topEffAccErr/mu_TTGamma_topEffAcc)**2 + 
				(mu_phoAccErr/mu_phoAcc)**2 + 
				(e_TTGamma_topEffAccErr/e_TTGamma_topEffAcc)**2 + 
				(e_phoAccErr/e_phoAcc)**2)**0.5
combined_TTJets_topEffAccErr = (mu_TTJets_topEffAccErr**2 + e_TTJets_topEffAccErr**2)**0.5
combined_TTgammaPhoEffAccVisErr =  ((mu_TTGammaVis_topAccErr/mu_TTGammaVis_topAcc)**2 + 
				    (mu_phoRecoEffErr/mu_phoRecoEff)**2 + 
				    (e_TTGammaVis_topAccErr/e_TTGammaVis_topAcc)**2 + 
				    (e_phoRecoEffErr/e_phoRecoEff)**2)**0.5



xsRatio = combined_ttgammaSig / combined_TTgammaPhoEffAcc / combined_topPreselInt * combined_TTJets_topEffAcc
xsRatioRelErr = (   (combined_ttgammaSigErr/combined_ttgammaSig)**2 + 
		    (combined_TTgammaPhoEffAccErr/combined_TTgammaPhoEffAcc)**2 +
		    (combined_topPreselErr/combined_topPreselInt)**2 +
		    (combined_TTJets_topEffAccErr/combined_TTJets_topEffAcc)**2
		    )**0.5

vis_xsRatio = combined_ttgammaSig / combined_TTgammaPhoEffAccVis / combined_topPreselInt * combined_TTJets_topEffAcc
vis_xsRatioRelErr = (   (combined_ttgammaSigErr/combined_ttgammaSig)**2 + 
		       (combined_TTgammaPhoEffAccVisErr/combined_TTgammaPhoEffAccVis)**2 +
		       (combined_topPreselErr/combined_topPreselInt)**2 +
		       (combined_TTJets_topEffAccErr/combined_TTJets_topEffAcc)**2
		       )**0.5


print '*'*80
print 'final answer: cross section ratio:'
print xsRatio,' +-',xsRatio*xsRatioRelErr
print '*'*80
print
print '*'*80
print 'visible cross section ratio:'
print vis_xsRatio,' +-',vis_xsRatio*vis_xsRatioRelErr
print '*'*80
