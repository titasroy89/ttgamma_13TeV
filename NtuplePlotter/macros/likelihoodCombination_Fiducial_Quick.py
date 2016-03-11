import ROOT
import sys
from math import exp
from math import log as ln

from Style import *
from array import array

import CMS_lumi
import os 
thestyle = Style()
 
HasCMSStyle = False
style = None
if os.path.isfile('tdrstyle.C'):
 	ROOT.gROOT.ProcessLine('.L tdrstyle.C')
        ROOT.setTDRStyle()
        print "Found tdrstyle.C file, using this style."
        HasCMSStyle = True
        if os.path.isfile('CMSTopStyle.cc'):
 		gROOT.ProcessLine('.L CMSTopStyle.cc+')
 		style = CMSTopStyle()
 		style.setupICHEPv1()
 		print "Found CMSTopStyle.cc file, use TOP style if requested in xml file."
if not HasCMSStyle:
 	print "Using default style defined in cuy package."
 	thestyle.SetStyle()
 
ROOT.gROOT.ForceStyle()
 #############
ROOT.gROOT.SetBatch(0)

myfile = None

eFactor = 1.
muFactor = 1.

eleFakeSF = 1.458
eleFakeSFErr = 0.1985

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

def findmax(e_pho, e_ele, e_fake, mu_pho, mu_ele, mu_fake, testttgSF = None, testvgSF = None, testjgSF = None):
    jumpVal = 0.2
    ttgSF = 1.
    vgSF  = 1.
    jgSF  = 1.
    if not testttgSF == None: ttgSF = testttgSF
    if not testvgSF == None:  vgSF = testvgSF
    if not testjgSF == None: jgSF = testjgSF
    for steps in range(25):
        localmax_lk = -1
        for i in range(3):
            if not testttgSF == None and i in [0,2]: continue
            for j in range(3):
                if not testvgSF == None and j in [0,2]: continue
                for k in range(3):
                    if not testjgSF == None and k in [0,2]: continue
                    lk = getLikelihood(e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake,ttgSF + (i-1)*jumpVal, vgSF + (j-1)*jumpVal, jgSF + (k-1)*jumpVal)
                    if lk > localmax_lk:
                        localmax_lk = lk
                        best_i = i
                        best_j = j
                        best_k = k
	if best_i==best_j==best_k==1:
            jumpVal = jumpVal/2.
        ttgSF = ttgSF + (best_i-1)*jumpVal
        vgSF = vgSF + (best_j-1)*jumpVal
        jgSF = jgSF + (best_k-1)*jumpVal
        
    return localmax_lk

def getJetGamma_OneSigma(bestjgSF, step, e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake):
	maxlk = findmax(e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake,testjgSF = bestjgSF)
	stepval = abs(step)
	testjg = bestjgSF+step
	for steps in range(10):
		best_i = -1
		minNll = 99.
		for i in range(3):
			nll = -1*ln(findmax(e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake, testjgSF = testjg+(i-1)*stepval)/maxlk)
			if abs(nll-0.5) < minNll:
				best_i = i
				minNll = abs(nll-0.5)
		if best_i==1:
			stepval = stepval/2.
		testjg = testjg+(best_i-1)*stepval
		if abs(testjg - bestjgSF)<1e-5:
			testjg = bestjgSF + stepval
			stepval = stepval/2.

	return testjg

def getVGamma_OneSigma(bestvgSF, step, e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake):
	maxlk = findmax(e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake, testvgSF =bestvgSF)
	stepval = abs(step)
	testvg = bestvgSF+step
	for steps in range(10):
		best_i = -1
		minNll = 99.
		for i in range(3):
			nll = -1*ln(findmax(e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake, testvgSF = testvg+(i-1)*stepval)/maxlk)
			if abs(nll-0.5) < minNll:
				best_i = i
				minNll = abs(nll-0.5)
		if best_i==1:
			stepval = stepval/2.
		testvg = testvg+(best_i-1)*stepval
		if abs(testvg - bestvgSF)<1e-5:
			testvg = bestvgSF + stepval
			stepval = stepval/2.

	return testvg

def getTTGamma_OneSigma(bestttgSF, step, e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake):
	maxlk = findmax(e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake, testttgSF = bestttgSF)
	stepval = abs(step)
	testttg = bestttgSF+step
	for steps in range(10):
		best_i = -1
		minNll = 99.
		for i in range(3):
			nll = -1*ln(findmax(e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake, testttgSF = testttg+(i-1)*stepval)/maxlk)
			if abs(nll-0.5) < minNll:
				best_i = i
				minNll = abs(nll-0.5)
		if best_i==1:
			stepval = stepval/2.
		testttg = testttg+(best_i-1)*stepval
		if abs(testttg - bestttgSF)<1e-5:
			testttg = bestttgSF + stepval
			stepval = stepval/2.

	return testttg



def maximizeLikelihood(e_pho, e_ele, e_fake, mu_pho, mu_ele, mu_fake):
    jumpVal = 0.2
    ttgSF = 1
    vgSF = 1
    jgSF = 1


    for steps in range(20):
        best_lk = -1
        best_vals = [-1,-1,-1]

        best_i = -1
        best_j = -1
	best_k = -1

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    lk = getLikelihood(e_pho, e_ele, e_fake, mu_pho, mu_ele, mu_fake,ttgSF+(i-1)*jumpVal, vgSF + (j-1)*jumpVal, jgSF + (k-1)*jumpVal)
                    if lk > best_lk:
			    best_lk = lk
			    best_i = i
			    best_j = j
			    best_k = k
			    
        if best_i==best_j==best_k==1:
            jumpVal = jumpVal/2.
        ttgSF = ttgSF+(best_i-1)*jumpVal
        vgSF = vgSF + (best_j-1)*jumpVal
        jgSF = jgSF + (best_k-1)*jumpVal
        maxlk = best_lk
        bestttgSF = ttgSF
        bestVgSF = vgSF
        bestjgSF = jgSF

    return bestttgSF, bestVgSF, bestjgSF, maxlk


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


ttgammaSeq = seq(0.2, 0.7, 0.01) + seq(0.705,1.15,0.005) + seq(1.16, 1.6,0.01) + seq(1.62,3,.02)

def calculateTTGamma(e_templateFile, mu_templateFile, combined_eff, saveFitPlots = False, verbose = True, progressBar = False,outputHistFileName = 'Likelihood_Histograms_Combination.root'):
	# read all the numbers from the file
	e_pho = readSamples('signal', e_templateFile, verbose)
	e_ele = readSamples('electron', e_templateFile, verbose)
	e_fake = readSamples('fake', e_templateFile, verbose)

	mu_pho = readSamples('signal', mu_templateFile, verbose)
	mu_ele = readSamples('electron', mu_templateFile, verbose)
	mu_fake = readSamples('fake', mu_templateFile, verbose)
	
#	ttgammaSeq = seq(0.2,1.6,0.2)

	# vgammaSeq = seq(0, 1., 0.05)
# 	vgammaSeq = seq(0.0,1.84,0.02)
# 	jetSeq = seq(0.5, 1.6, 0.01)

# 	ttgammaBinSize = ttgammaSeq[1]-ttgammaSeq[0]
# 	ttgammabinning = [i-0.0025 for i in ttgammaSeq]+[ttgammaSeq[-1]+0.0075]
	

# 	ttghist = ROOT.TH1F('ttghist','Marginalized Likelihood', len(ttgammaSeq), array('d',ttgammabinning))
# #	ttghist = ROOT.TH1F('ttghist','Marginalized Likelihood', 141, 0.2-0.005, 1.6+0.005)
# 	vghist = ROOT.TH1F('vghist','Marginalized Likelihood', len(vgammaSeq), vgammaSeq[0]-0.025, vgammaSeq[-1]+0.025)
# 	jghist = ROOT.TH1F('jghist','Marginalized Likelihood', len(jetSeq), jetSeq[0]-0.005, jetSeq[-1]+0.005)

# 	ttghist_nll = ROOT.TH1F('ttghist_nll','Log Likelihood', len(ttgammaSeq), array('d',ttgammabinning))
# #	ttghist = ROOT.TH1F('ttghist','Marginalized Likelihood', 141, 0.2-0.005, 1.6+0.005)
# 	vghist_nll = ROOT.TH1F('vghist_nll','Log Likelihood', len(vgammaSeq), vgammaSeq[0]-0.025, vgammaSeq[-1]+0.025)
# 	jghist_nll = ROOT.TH1F('jghist_nll','Log Likelihood', len(jetSeq), jetSeq[0]-0.005, jetSeq[-1]+0.005)

# 	if saveFitPlots:
# 		ttg_vg_hist = ROOT.TH2F('ttg_vg','Marginalized Likelihood', 141, 0.2-0.005, 1.6+0.005, 61, 0-0.025, 3.0+0.025)
# 		vg_jg_hist  = ROOT.TH2F('vg_jg','Marginalized Likelihood', 61, 0-0.025, 3.0+0.025, 131, 0.5-0.005, 1.8+0.005)
# 		jg_ttg_hist = ROOT.TH2F('jg_ttg','Marginalized Likelihood', 131, 0.5-0.005, 1.8+0.005, 141, 0.2-0.005, 1.6+0.005)
		
# #		ttg_vg_jg_hist = ROOT.TH3F('ttg_vg_jg','Marginalized Likelihood', 141, 0.2-0.005, 1.6+0.005, 61, 0-0.025, 3.0+0.025, 131, 0.5-0.005, 1.8+0.005)
# 		tempVgammaSeq = [b-.01 for b in vgammaSeq]+[2*vgammaSeq[-1]-vgammaSeq[-2]-.01]
# 		tempJetSeq = [b-.005 for b in jetSeq]+[2*jetSeq[-1]-jetSeq[-2]-.005]
# 		ttg_vg_jg_hist = ROOT.TH3F('ttg_vg_jg','Marginalized Likelihood',len(ttgammaSeq), array('d',ttgammabinning),len(vgammaSeq), array('d',tempVgammaSeq), len(jetSeq), array('d',tempJetSeq))



	maxlk = -1.0
	bestttgSF = -1.0
	bestvgSF = -1.0
	bestjgSF = -1.0
	
	percent = 0
	i = 0.
	totalTTgammaPoints = len(ttgammaSeq)
#	for ttgSF in seq(0.2, 1.6, 0.01):		
	
	

	bestttgSF, bestvgSF, bestjgSF, maxlk = maximizeLikelihood(e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake)

	up   = getTTGamma_OneSigma(bestttgSF, 0.1, e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake)
	down = getTTGamma_OneSigma(bestttgSF,-0.1, e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake)
	bestttgSFErr = (up-down)/2.

	print bestttgSF

	up   = getVGamma_OneSigma(bestvgSF, 0.1, e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake)
	down = getVGamma_OneSigma(bestvgSF,-0.1, e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake)
	bestvgSFErr = (up-down)/2.

	up   = getJetGamma_OneSigma(bestjgSF, 0.1, e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake)
	down = getJetGamma_OneSigma(bestjgSF,-0.1, e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake)
	bestjgSFErr = (up-down)/2.



	# for ttgSF in ttgammaSeq:
	# 	if progressBar:
	# 		if (i/totalTTgammaPoints)>(percent/100.):
	# 			sys.stdout.write("\r[" + "=" * int(percent/4) + " " * (25- int(percent/4)) + "]" + str(percent) + "%")
	# 			sys.stdout.flush()
	# 			percent += 1
				
	# 		# if (i%7==0):
	# 		# 	sys.stdout.write("\r[" + "=" * (i / 7) + " " * (20- (i/7)) + "]" + str(percent) + "%")
	# 		# 	sys.stdout.flush()
	# 		# 	percent += 5
	# 	i += 1.
	# 	for VgSF in vgammaSeq:
	# 		for jgSF in jetSeq:
	# 			lk = getLikelihood(e_pho,e_ele,e_fake,mu_pho,mu_ele,mu_fake,ttgSF, VgSF, jgSF)
	# 			ttghist.Fill(ttgSF,lk*0.0005)
	# 			vghist.Fill(VgSF,lk*0.0001)
	# 			jghist.Fill(jgSF,lk*0.0005)

	# 			ttghist_nll.Fill(ttgSF,-lk)
	# 			vghist_nll.Fill(VgSF,-lk)
	# 			jghist_nll.Fill(jgSF,-lk)

	# 			if saveFitPlots:
	# 				ttg_vg_hist.Fill(ttgSF, VgSF, lk/131.)
	# 				vg_jg_hist.Fill(VgSF,   jgSF, lk/141.)
	# 				jg_ttg_hist.Fill(jgSF, ttgSF, lk/61.)
	# 				ttg_vg_jg_hist.Fill(ttgSF, VgSF, jgSF, lk)


	# 			if lk > maxlk:
	# 				maxlk = lk
	# 				bestttgSF = ttgSF
	# 				bestvgSF = VgSF
	# 				bestjgSF = jgSF
	
	print
	print 'max likelihood ',maxlk
	print "best SF's:", bestttgSF, bestvgSF, bestjgSF

	b_e_allmc,e_allmcErr = NallMC(e_pho,e_ele,e_fake, bestttgSF, bestvgSF, bestjgSF)
	b_e_top,e_topErr = Ntop(e_pho,e_ele,e_fake, bestttgSF, bestvgSF, bestjgSF)
	b_e_egamma,e_egammaErr = Negamma(e_pho,e_ele,e_fake, bestttgSF, bestvgSF, bestjgSF)

	b_mu_allmc,mu_allmcErr = NallMC(mu_pho,mu_ele,mu_fake, bestttgSF, bestvgSF, bestjgSF)
	b_mu_top,mu_topErr = Ntop(mu_pho,mu_ele,mu_fake, bestttgSF, bestvgSF, bestjgSF)
	b_mu_egamma,mu_egammaErr = Negamma(mu_pho,mu_ele,mu_fake, bestttgSF, bestvgSF, bestjgSF)


# 	ROOT.gStyle.SetOptFit(111)
# 	ROOT.gStyle.SetOptStat(0)
# 	H = 600; 
# 	W = 800; 

# 	canvas = ROOT.TCanvas('c1','c1',W,H)


# 	# references for T, B, L, R
# 	T = 0.08*H
# 	B = 0.12*H 
# 	L = 0.12*W
# 	R = 0.04*W
# 	canvas.SetFillColor(0)
# 	canvas.SetBorderMode(0)
# 	canvas.SetFrameFillStyle(0)
# 	canvas.SetFrameBorderMode(0)
# 	canvas.SetLeftMargin( L/W )
# 	canvas.SetRightMargin( R/W )
# 	canvas.SetTopMargin( T/H )
# 	canvas.SetBottomMargin( B/H )
# 	canvas.SetTickx(0)
# 	canvas.SetTicky(0)

# 	ttghist.Draw()
# #	ttghist.GetXaxis().SetTitle('SF_{t#bar{t}+#gamma}')
# 	ttghist.GetXaxis().SetTitle('SF(t#bar{t}+#gamma)')
# 	ttghist.GetYaxis().SetTitle('Likelihood')
# 	ttghist.Fit('gaus',"Q")
# 	fit = ttghist.GetFunction('gaus')
# 	bestttgSFErr = fit.GetParameter(2)

# 	CMS_lumi.writeExtraText = True
# 	CMS_lumi.writeChannelText = True
# 	CMS_lumi.channelText = 'e/#mu+jets'
# 	CMS_lumi.channelTextLocation = 1

# 	CMS_lumi.CMS_lumi(canvas, 2, 11)
# 	canvas.Update()
# 	if saveFitPlots: 
# 		canvas.Print('TTGamma_SF_Lkhood.pdf','.pdf')
# 		canvas.Print('TTGamma_SF_Lkhood.png','.png')

# 	vghist.Draw()
# #	vghist.GetXaxis().SetTitle('SF_{V+#gamma}')
# 	vghist.GetXaxis().SetTitle('SF(V+#gamma)')
# #	vghist.GetXaxis().SetTitle('Vgamma Scale Factor')
# 	vghist.GetYaxis().SetTitle('Likelihood')
# 	vghist.SetMinimum(0.0)
# 	vghist.Fit('gaus',"Q")
# 	fit = vghist.GetFunction('gaus')
# 	bestvgSFErr = fit.GetParameter(2)
# 	CMS_lumi.writeExtraText = True
# 	CMS_lumi.writeChannelText = True
# 	CMS_lumi.channelText = 'e/#mu+jets'
# 	CMS_lumi.channelTextLocation = 1

# 	CMS_lumi.CMS_lumi(canvas, 2, 11)
# 	canvas.Update()
# 	if saveFitPlots: 
# 		canvas.Print('Vgamma_SF_Lkhood.pdf','.pdf')
# 		canvas.Print('Vgamma_SF_Lkhood.png','.png')

	
# 	jghist.Draw()
# #	jghist.GetXaxis().SetTitle('SF_{jet#rightarrow#gamma}')
# 	jghist.GetXaxis().SetTitle('SF(jet#rightarrow#gamma)')
# #	jghist.GetXaxis().SetTitle('Jet to Photon Scale Factor')
# 	jghist.GetYaxis().SetTitle('Likelihood')
# 	jghist.Fit('gaus',"Q")
# 	fit = jghist.GetFunction('gaus')
# 	bestjgSFErr = fit.GetParameter(2)
# 	CMS_lumi.writeExtraText = True
# 	CMS_lumi.writeChannelText = True
# 	CMS_lumi.channelText = 'e/#mu+jets'
# 	CMS_lumi.channelTextLocation = 1

# 	CMS_lumi.CMS_lumi(canvas, 2, 11)
# 	canvas.Update()
# 	if saveFitPlots: 
# 		canvas.Print('jet_gamma_SF_Lkhood.pdf','.pdf')
# 		canvas.Print('jet_gamma_SF_Lkhood.png','.png')
	
	ttgammaSig = bestttgSF*(eFactor * (e_pho['TTGamma'][0] + e_ele['TTGamma'][0] + e_fake['TTGamma'][0]) + muFactor * (mu_pho['TTGamma'][0] + mu_ele['TTGamma'][0] + mu_fake['TTGamma'][0]))
	ttgammaSigErr = bestttgSFErr*(eFactor * (e_pho['TTGamma'][0] + e_ele['TTGamma'][0] + e_fake['TTGamma'][0]) + muFactor * (mu_pho['TTGamma'][0] + mu_ele['TTGamma'][0] + mu_fake['TTGamma'][0]))

	# if saveFitPlots:
	# 	out = ROOT.TFile(outputHistFileName,'recreate')
	# 	ttghist.Write()
	# 	vghist.Write()
	# 	jghist.Write()
	# 	ttgfit = ttghist.GetFunction('gaus')
	# 	ttgfit.SetNameTitle("ttgFit","ttgFit")
	# 	vgfit = vghist.GetFunction('gaus')
	# 	vgfit.SetNameTitle("vgFit","vgFit")
	# 	jgfit = jghist.GetFunction('gaus')
	# 	jgfit.SetNameTitle("jgFit","jgFit")
	# 	ttgfit.Write()
	# 	vgfit.Write()
	# 	jgfit.Write()

	# 	ttg_vg_hist.Write()
	# 	vg_jg_hist.Write()
	# 	jg_ttg_hist.Write()
	# 	ttg_vg_jg_hist.Write()
	# 	out.Close()

	
	xsRatio = ttgammaSig / combined_eff['FidEff'] * combined_eff['TTJets_topEffAcc'] / combined_eff['topPreselInt'] 
	xsRatioRelErr = (   (ttgammaSigErr/ttgammaSig)**2 + 
			    (combined_eff['FidEffErr']/combined_eff['FidEff'])**2 + 
			    (combined_eff['topPreselErr']/combined_eff['topPreselInt'])**2 + 
			    (combined_eff['TTJets_topEffAccErr']/combined_eff['TTJets_topEffAcc'])**2
			    )**0.5

	from SF import luminosity
	xsDirect = ttgammaSig / combined_eff['FidEff'] / luminosity
	xsDirectErr =  xsDirect * ( (ttgammaSigErr/ttgammaSig)**2 + 
				    (combined_eff['FidEffErr']/combined_eff['FidEff'])**2
				    )**0.5
	
	
	return (xsRatio, xsRatio*xsRatioRelErr), (xsDirect, xsDirectErr), [(bestttgSF, bestttgSFErr), (bestvgSF, bestvgSFErr), (bestjgSF, bestjgSFErr), (ttgammaSig, ttgammaSigErr)]
