import ROOT
import sys
from math import exp

ROOT.gROOT.SetBatch(0)

myfile = None

eleFakeSF = 1.5
eleFakeSFErr = 0.15

# parameters used in chi square calculation
e_photnPurity	       =  0.564158170272
e_photnPurityErr       =  0.0646043755
e_M3_photon_topFrac    =  0.720712773763
e_M3_photon_topFracErr =  0.0609031509759
e_Ndata		       =  977.0
e_NdataErr             =  31.2569992162

# parameters used in chi square calculation
mu_photnPurity	        =  0.531874490818
mu_photnPurityErr       =  0.0575003678412
mu_M3_photon_topFrac    =  0.696751649028
mu_M3_photon_topFracErr =  0.0585215226403
mu_Ndata	        =  1173.0
mu_NdataErr             =  34.2490875791

e_data  = {'photnPurity'          : e_photnPurity	  ,     
	   'photnPurityErr'       : e_photnPurityErr      ,
	   'M3_photon_topFrac'    : e_M3_photon_topFrac   ,
	   'M3_photon_topFracErr' : e_M3_photon_topFracErr,
	   'Ndata'                : e_Ndata	          ,
	   'NdataErr'             : e_NdataErr            ,
	   }

mu_data  = {'photnPurity'          : mu_photnPurity         ,     
	    'photnPurityErr'       : mu_photnPurityErr      ,
	    'M3_photon_topFrac'    : mu_M3_photon_topFrac   ,
	    'M3_photon_topFracErr' : mu_M3_photon_topFracErr,
	    'Ndata'                : mu_Ndata               ,
	    'NdataErr'             : mu_NdataErr            ,
	    }

def integral(histname, fileName, verbose = True):
	global myfile
	# if myfile == None:
	# 	myfile = ROOT.TFile(fileName,'READ')		
	myfile = ROOT.TFile(fileName,'READ')
	err = ROOT.Double(0.0)
	hist = myfile.Get(histname)
	integr = hist.IntegralAndError(0,hist.GetNbinsX() + 1, err)
	if verbose: print integr,err
	return integr,err


def readSamples(suffix, fileName, verbose = True):
	var = 'MET'
	######## Added QCD to the list of samples ########
	samplnames = ['TTGamma', 'TTJets', 'Vgamma', 'SingleTop', 'WJets', 'ZJets', 'QCD']
	samples = {}
	for n in samplnames:
		if n=='QCD' and (suffix=='signal' or suffix=='electron'):
			samples[n] = (0.0,0.0)
			continue
		if verbose: print 'getting counts for ',n,suffix###, 'in', fileName
		if n=='QCD':
			int,err = integral(n+'_'+var, fileName, verbose)
		else:
			int,err = integral(n+'_'+suffix+'_'+var, fileName, verbose)
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
	cs += (allmc - e_data['Ndata'])**2 / (e_data['NdataErr']**2 + allmcErr**2)
	
	mcPurErr = (egamma/allmc) * ((egammaErr/egamma)**2 + (allmcErr/allmc)**2)**0.5 
	cs += (egamma/allmc - e_data['photnPurity'])**2 / (e_data['photnPurityErr']**2 + mcPurErr**2)
	
	mcTopErr = (top/allmc) * ((topErr/top)**2 + (allmcErr/allmc)**2)**0.5
	cs += (top/allmc - e_data['M3_photon_topFrac'])**2 / (e_data['M3_photon_topFracErr']**2 + mcTopErr**2)
	return cs

def getChiSq_mu(pho,ele,fake, ttgSF, VgSF, jgSF):
	allmc,allmcErr = NallMC(pho,ele,fake, ttgSF, VgSF, jgSF)
	top,topErr = Ntop(pho,ele,fake, ttgSF, VgSF, jgSF)
	egamma,egammaErr = Negamma(pho,ele,fake, ttgSF, VgSF, jgSF)
	cs = 0.0
	cs += (allmc - mu_data['Ndata'])**2 / (mu_data['NdataErr']**2 + allmcErr**2)
	
	mcPurErr = (egamma/allmc) * ((egammaErr/egamma)**2 + (allmcErr/allmc)**2)**0.5 
	cs += (egamma/allmc - mu_data['photnPurity'])**2 / (mu_data['photnPurityErr']**2 + mcPurErr**2)
	
	mcTopErr = (top/allmc) * ((topErr/top)**2 + (allmcErr/allmc)**2)**0.5
	cs += (top/allmc - mu_data['M3_photon_topFrac'])**2 / (mu_data['M3_photon_topFracErr']**2 + mcTopErr**2)
	return cs


def getLikelihood(e_pho,e_ele,e_fake, mu_pho, mu_ele,mu_fake, ttgSF, VgSF, jgSF):
	return exp(-0.5*(getChiSq_e(e_pho,e_ele,e_fake, ttgSF, VgSF, jgSF)+getChiSq_mu(mu_pho,mu_ele,mu_fake, ttgSF, VgSF, jgSF)))

def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])
    else:
        return([])


def calculateTTGamma(e_templateFile, mu_templateFile, combined_eff, saveFitPlots = False, verbose = True, progressBar = False):
	# read all the numbers from the file
	e_pho = readSamples('signal', e_templateFile, verbose)
	e_ele = readSamples('electron', e_templateFile, verbose)
	e_fake = readSamples('fake', e_templateFile, verbose)

	mu_pho = readSamples('signal', mu_templateFile, verbose)
	mu_ele = readSamples('electron', mu_templateFile, verbose)
	mu_fake = readSamples('fake', mu_templateFile, verbose)
	
	ttghist = ROOT.TH1F('ttghist','Marginalized Likelihood', 141, 0.2-0.005, 1.6+0.005)
	vghist = ROOT.TH1F('vghist','Marginalized Likelihood', 61, 0-0.025, 3.0+0.025)
	jghist = ROOT.TH1F('jghist','Marginalized Likelihood', 131, 0.5-0.005, 1.8+0.005)

	if saveFitPlots:
		ttg_vg_hist = ROOT.TH2F('ttg_vg','Marginalized Likelihood', 141, 0.2-0.005, 1.6+0.005, 61, 0-0.025, 3.0+0.025)
		vg_jg_hist  = ROOT.TH2F('vg_jg','Marginalized Likelihood', 61, 0-0.025, 3.0+0.025, 131, 0.5-0.005, 1.8+0.005)
		jg_ttg_hist = ROOT.TH2F('jg_ttg','Marginalized Likelihood', 131, 0.5-0.005, 1.8+0.005, 141, 0.2-0.005, 1.6+0.005)
		
		ttg_vg_jg_hist = ROOT.TH3F('ttg_vg_jg','Marginalized Likelihood', 141, 0.2-0.005, 1.6+0.005, 61, 0-0.025, 3.0+0.025, 131, 0.5-0.005, 1.8+0.005)



	maxlk = -1.0
	bestttgSF = -1.0
	bestvgSF = -1.0
	bestjgSF = -1.0
	
	percent = 0
	i = 0
	for ttgSF in seq(0.2, 1.6, 0.01):		
		if progressBar:
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

				if saveFitPlots:
					ttg_vg_hist.Fill(ttgSF, VgSF, lk/131.)
					vg_jg_hist.Fill(VgSF,   jgSF, lk/141.)
					jg_ttg_hist.Fill(jgSF, ttgSF, lk/61.)
					ttg_vg_jg_hist.Fill(ttgSF, VgSF, jgSF, lk)


				if lk > maxlk:
					maxlk = lk
					bestttgSF = ttgSF
					bestvgSF = VgSF
					bestjgSF = jgSF
	
	if verbose:
		print 'max likelihood ',maxlk
		print "best SF's:", bestttgSF, bestvgSF, bestjgSF

	b_e_allmc,e_allmcErr = NallMC(e_pho,e_ele,e_fake, bestttgSF, bestvgSF, bestjgSF)
	b_e_top,e_topErr = Ntop(e_pho,e_ele,e_fake, bestttgSF, bestvgSF, bestjgSF)
	b_e_egamma,e_egammaErr = Negamma(e_pho,e_ele,e_fake, bestttgSF, bestvgSF, bestjgSF)

	b_mu_allmc,mu_allmcErr = NallMC(mu_pho,mu_ele,mu_fake, bestttgSF, bestvgSF, bestjgSF)
	b_mu_top,mu_topErr = Ntop(mu_pho,mu_ele,mu_fake, bestttgSF, bestvgSF, bestjgSF)
	b_mu_egamma,mu_egammaErr = Negamma(mu_pho,mu_ele,mu_fake, bestttgSF, bestvgSF, bestjgSF)

	if saveFitPlots:
		out = ROOT.TFile('Likelihood_Histograms.root','recreate')
		ttghist.Write()
		vghist.Write()
		jghist.Write()
		ttg_vg_hist.Write()
		vg_jg_hist.Write()
		jg_ttg_hist.Write()
		ttg_vg_jg_hist.Write()
		out.Close()


	ROOT.gStyle.SetOptFit(111)
	ccc = ROOT.TCanvas()

	ttghist.Draw()
	ttghist.GetXaxis().SetTitle('TTGamma Scale Factor')
	ttghist.Fit('gaus',"Q")
	fit = ttghist.GetFunction('gaus')
	bestttgSFErr = fit.GetParameter(2)
	if saveFitPlots: ccc.SaveAs('TTGamma_SF_Lkhood.png')

	vghist.Draw()
	vghist.GetXaxis().SetTitle('Vgamma Scale Factor')
	vghist.SetMinimum(0.0)
	vghist.Fit('gaus',"Q")
	fit = vghist.GetFunction('gaus')
	bestvgSFErr = fit.GetParameter(2)
	if saveFitPlots: ccc.SaveAs('Vgamma_SF_Lkhood.png')
	
	jghist.Draw()
	jghist.GetXaxis().SetTitle('Jet to Photon Scale Factor')
	jghist.Fit('gaus',"Q")
	fit = jghist.GetFunction('gaus')
	bestjgSFErr = fit.GetParameter(2)
	if saveFitPlots: ccc.SaveAs('jet_gamma_SF_Lkhood.png')
	
	ttgammaSig = bestttgSF*(e_pho['TTGamma'][0]+mu_pho['TTGamma'][0])
	ttgammaSigErr = bestttgSFErr*(e_pho['TTGamma'][0]+mu_pho['TTGamma'][0])

	

	xsRatio = ttgammaSig / combined_eff['TTgammaPhoEffAcc'] * combined_eff['TTJets_topEffAcc'] / combined_eff['topPreselInt'] 
	xsRatioRelErr = (   (ttgammaSigErr/ttgammaSig)**2 + 
			    (combined_eff['TTgammaPhoEffAccErr']/combined_eff['TTgammaPhoEffAcc'])**2 + 
			    (combined_eff['topPreselErr']/combined_eff['topPreselInt'])**2 + 
			    (combined_eff['TTJets_topEffAccErr']/combined_eff['TTJets_topEffAcc'])**2
			    )**0.5

	vis_xsRatio = ttgammaSig / combined_eff['TTgammaPhoEffAccVis'] * combined_eff['TTJets_topEffAcc'] / combined_eff['topPreselInt'] 
	vis_xsRatioRelErr = ( (ttgammaSigErr/ttgammaSig)**2 + 
			      (combined_eff['TTgammaPhoEffAccVisErr']/combined_eff['TTgammaPhoEffAccVis'])**2 + 
			      (combined_eff['topPreselErr']/combined_eff['topPreselInt'])**2 + 
			      (combined_eff['TTJets_topEffAccErr']/combined_eff['TTJets_topEffAcc'])**2
			      )**0.5
	
	return (xsRatio, xsRatio*xsRatioRelErr), (vis_xsRatio, vis_xsRatio*vis_xsRatioRelErr), [(bestttgSF, bestttgSFErr), (bestvgSF, bestvgSFErr), (bestjgSF, bestjgSFErr), (ttgammaSig, ttgammaSigErr)]
