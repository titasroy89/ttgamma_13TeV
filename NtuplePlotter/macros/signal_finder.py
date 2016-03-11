import os
import ROOT
from math import exp
from math import log as ln
import sys

from Style import *

import CMS_lumi
 
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

barrelFileName = 'templates_barrel_scaled.root'
myfile = None

def integral(histname):
	global myfile
	if myfile == None:
		myfile = ROOT.TFile(barrelFileName,'READ')
	err = ROOT.Double(0.0)
	hist = myfile.Get(histname)
	integr = hist.IntegralAndError(0,hist.GetNbinsX() + 1, err)
#	print integr,err
	return integr,err

def findBounds(f, y=0.5):
        high = f.GetMinimumX()
	low = high-1
	for i in range(1,20):
                iterVal = (high-low)/2.
                midPoint = high-iterVal
                mid = f.Eval(midPoint)
                if mid > y:
                        low = midPoint
                else:
                        high = midPoint
	lowBound = midPoint

	low = f.GetMinimumX()
        high = low+1
        for i in range(1,20):
                iterVal = (high-low)/2.
                midPoint = high-iterVal
                mid = f.Eval(midPoint)
		if mid > y:
			high = midPoint
		else:
			low = midPoint
        highBound = midPoint

        return lowBound, highBound


eleFakeSF = 1.5 
eleFakeSFErr = 0.2

# parameters used in chi square calculation
photnPurity = 0.657
photnPurityErr = 0.056

M3_photon_topFrac =  0.68
M3_photon_topFracErr = 0.12

Ndata = 1203.0  
NdataErr = 34.68


def readSamples(suffix):
	var = 'MET'
	######## Added QCD to the list of samples ########
	samplnames = ['TTGamma', 'TTJets', 'Vgamma', 'SingleTop', 'WJets', 'ZJets', 'QCD']
	samples = {}
	for n in samplnames:
		if n=='QCD' and (suffix=='signal' or suffix=='electron'):
			samples[n] = (0.0,0.0)
			continue
		#print 'getting counts for ',n,suffix
		if n=='QCD':
			int,err = integral(n+'_'+var)
		else:
			int,err = integral(n+'_'+suffix+'_'+var)
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



def getChiSq(pho,ele,fake, ttgSF, VgSF, jgSF):
	allmc,allmcErr = NallMC(pho,ele,fake, ttgSF, VgSF, jgSF)
	top,topErr = Ntop(pho,ele,fake, ttgSF, VgSF, jgSF)
	egamma,egammaErr = Negamma(pho,ele,fake, ttgSF, VgSF, jgSF)
	cs = 0.0
	cs += (allmc - Ndata)**2 / (NdataErr**2 + allmcErr**2)
	
	mcPurErr = (egamma/allmc) * ((egammaErr/egamma)**2 + (allmcErr/allmc)**2)**0.5 
	cs += (egamma/allmc - photnPurity)**2 / (photnPurityErr**2 + mcPurErr**2)
	
	mcTopErr = (top/allmc) * ((topErr/top)**2 + (allmcErr/allmc)**2)**0.5
	cs += (top/allmc - M3_photon_topFrac)**2 / (M3_photon_topFracErr**2 + mcTopErr**2)
	return cs

def getLikelihood(pho,ele,fake, ttgSF, VgSF, jgSF):
	return exp(-0.5*getChiSq(pho,ele,fake, ttgSF, VgSF, jgSF))

def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])
    else:
        return([])


ttgammaSeq = seq(0.2, 0.7, 0.01) + seq(0.705,1.0,0.005) + seq(1.01, 1.6,0.01) + seq(1.62,3,.02)
vgammaSeq = seq(0.0,1.84,0.02)
jetSeq = seq(0.5, 1.6, 0.01)

def calculateTTGamma(progressBar = False, saveDataFile = False):
	# read all the numbers from the file
	pho = readSamples('signal')
	ele = readSamples('electron')
	fake = readSamples('fake')

	if saveDataFile:
		output = open('datafile.py','w')
		line  = str(pho)+'\n'
		line += str(ele)+'\n'
		line += str(fake)+'\n'
		line += '%f\n'%photnPurity
		line += '%f\n'%photnPurityErr
		line += '%f\n'%M3_photon_topFrac
		line += '%f\n'%M3_photon_topFracErr
		line += '%f\n'%Ndata
		line += '%f\n'%NdataErr
		output.write(line)
		output.close()

#	binning = [i-0.0025 for i in ttgammaSeq]+[ttgammaSeq[-1]+0.0075]
	ttgammabinning = [i-0.0025 for i in ttgammaSeq]+[ttgammaSeq[-1]+0.0075]
	tempVgammaSeq = [b-.01 for b in vgammaSeq]+[2*vgammaSeq[-1]-vgammaSeq[-2]-.01]
	tempJetSeq = [b-.005 for b in jetSeq]+[2*jetSeq[-1]-jetSeq[-2]-.005]
	

	from array import array
	ttghist = ROOT.TH1F('ttghist','Marginalized Likelihood', len(ttgammaSeq), array('d',ttgammabinning))
	vghist = ROOT.TH1F('vghist','Marginalized Likelihood', len(vgammaSeq), vgammaSeq[0]-0.025, vgammaSeq[-1]+0.025)
	jghist = ROOT.TH1F('jghist','Marginalized Likelihood', len(jetSeq), jetSeq[0]-0.005, jetSeq[-1]+0.005)

	ttg_vg_jg_hist = ROOT.TH3F('ttg_vg_jg','Marginalized Likelihood',len(ttgammaSeq), array('d',ttgammabinning),len(vgammaSeq), array('d',tempVgammaSeq), len(jetSeq), array('d',tempJetSeq))

	ttghist_NLL = ROOT.TH1F('ttghist_NLL','Negative Log Likelihood', len(ttgammaSeq), array('d',ttgammabinning))
	vghist_NLL = ROOT.TH1F('vghist_NLL','Negative Log Likelihood', len(vgammaSeq), vgammaSeq[0]-0.01, vgammaSeq[-1]+0.01)
	jghist_NLL = ROOT.TH1F('jghist_NLL','Negative Log Likelihood', len(jetSeq), jetSeq[0]-0.005, jetSeq[-1]+0.005)

	#ttghist = ROOT.TH1F('ttghist','Marginalized Likelihood', 141, 0.2-0.005, 1.6+0.005)
	# vghist = ROOT.TH1F('vghist','Marginalized Likelihood', 61, 0-0.025, 3.0+0.025)
	# jghist = ROOT.TH1F('jghist','Marginalized Likelihood', 131, 0.5-0.005, 1.8+0.005)
	
	maxlk = -1.0
	bestttgSF = -1.0
	bestVgSF = -1.0
	bestjgSF = -1.0
	

	percent = 0
	i = 0.
	totalTTgammaPoints = len(ttgammaSeq)

	#for ttgSF in seq(0.2, 1.6, 0.01):
	for ttgSF in ttgammaSeq:
		if progressBar:
			if (i/totalTTgammaPoints)>(percent/100.):
				sys.stdout.write("\r[" + "=" * int(percent/4) + " " * (25- int(percent/4)) + "]" + str(percent) + "%")
				sys.stdout.flush()
				percent += 1
		i += 1.
		for VgSF in vgammaSeq:
			for jgSF in jetSeq:
		# for VgSF in seq(0.0, 3.0, 0.05):
		# 	for jgSF in seq(0.5, 1.8, 0.01):
				lk = getLikelihood(pho,ele,fake,ttgSF, VgSF, jgSF)
				ttghist.Fill(ttgSF,lk*0.0005)
				vghist.Fill(VgSF,lk*0.0001)
				jghist.Fill(jgSF,lk*0.0005)
				ttg_vg_jg_hist.Fill(ttgSF, VgSF, jgSF, lk)
				if lk > maxlk:
					maxlk = lk
					bestttgSF = ttgSF
					bestVgSF = VgSF
					bestjgSF = jgSF
	
	print 'max likelihood ',maxlk
	print 'best ttgSF',bestttgSF
	print 'best VgSF',bestVgSF
	print 'best jgSF',bestjgSF

	b_allmc,allmcErr = NallMC(pho,ele,fake, bestttgSF, bestVgSF, bestjgSF)
	b_top,topErr = Ntop(pho,ele,fake, bestttgSF, bestVgSF, bestjgSF)
	b_egamma,egammaErr = Negamma(pho,ele,fake, bestttgSF, bestVgSF, bestjgSF)
	print 'total number of MC events ',b_allmc, '  Data:',Ndata
	print 'photon purity in MC ',b_egamma/b_allmc, '  template fit:',photnPurity
	print 'top fraction in MC ',b_top/b_allmc, '  M3 fit:',M3_photon_topFrac
	
	ROOT.gStyle.SetOptFit(111)
	ROOT.gStyle.SetOptStat(0)
	H = 600; 
	W = 800; 

	canvas = ROOT.TCanvas('c1','c1',W,H)


	# references for T, B, L, R
	T = 0.08*H
	B = 0.12*H 
	L = 0.12*W
	R = 0.04*W
	canvas.SetFillColor(0)
	canvas.SetBorderMode(0)
	canvas.SetFrameFillStyle(0)
	canvas.SetFrameBorderMode(0)
	canvas.SetLeftMargin( L/W )
	canvas.SetRightMargin( R/W )
	canvas.SetTopMargin( T/H )
	canvas.SetBottomMargin( B/H )
	canvas.SetTickx(0)
	canvas.SetTicky(0)

	ttghist.Draw()
	ttghist.GetXaxis().SetTitle('SF(t#bar{t}+#gamma)')
	ttghist.GetYaxis().SetTitle('Likelihood')
	ttghist.Fit('gaus')

	CMS_lumi.CMS_lumi(canvas, 2, 11)
	canvas.Print('plots/TTGamma_SF_Lkhood.png', '.png')
	canvas.Print('plots/TTGamma_SF_Lkhood.pdf', '.pdf')
	fit = ttghist.GetFunction('gaus')
	bestttgSFErr = fit.GetParameter(2)
	
	vghist.Draw()
	vghist.GetXaxis().SetTitle('SF(V+#gamma)')
	vghist.GetYaxis().SetTitle('Likelihood')
	vghist.SetMinimum(0.0)
	vghist.Fit('gaus')
	CMS_lumi.CMS_lumi(canvas, 2, 11)
	canvas.Print('plots/Vgamma_SF_Lkhood.png', '.png')
	canvas.Print('plots/Vgamma_SF_Lkhood.pdf', '.pdf')
	fit = vghist.GetFunction('gaus')
	bestVgSFErr = fit.GetParameter(2)
	
	jghist.Draw()
	jghist.GetXaxis().SetTitle('SF(jet#rightarrow#gamma#)')
	jghist.GetYaxis().SetTitle('Likelihood')
	jghist.Fit('gaus')
	CMS_lumi.CMS_lumi(canvas, 2, 11)
	canvas.Print('plots/jet_gamma_SF_Lkhood.png', '.png')
	canvas.Print('plots/jet_gamma_SF_Lkhood.pdf', '.pdf')
	fit = jghist.GetFunction('gaus')
	bestjgSFErr = fit.GetParameter(2)


	maxL = ttg_vg_jg_hist.GetMaximum()
	
	for ibin in range(1,ttghist_NLL.GetNbinsX()+1):
		ttg_vg_jg_hist.GetXaxis().SetRange(ibin,ibin)
		binL = ttg_vg_jg_hist.Project3D('yz').GetMaximum()
		if binL==0: binL=1e-20
		ttghist_NLL.SetBinContent(ibin,-1*ln(binL/maxL))
	ttg_vg_jg_hist.GetXaxis().SetRange(-1,-1)
	ttg_vg_jg_hist.GetYaxis().SetRange(-1,-1)
	ttg_vg_jg_hist.GetZaxis().SetRange(-1,-1)

	for ibin in range(1,vghist_NLL.GetNbinsX()+1):
		ttg_vg_jg_hist.GetYaxis().SetRange(ibin,ibin)
		binL = ttg_vg_jg_hist.Project3D('xz').GetMaximum()
		if binL==0: binL=1e-20
		vghist_NLL.SetBinContent(ibin,-1*ln(binL/maxL))
	ttg_vg_jg_hist.GetXaxis().SetRange(-1,-1)
	ttg_vg_jg_hist.GetYaxis().SetRange(-1,-1)
	ttg_vg_jg_hist.GetZaxis().SetRange(-1,-1)


	for ibin in range(1,jghist_NLL.GetNbinsX()+1):
		ttg_vg_jg_hist.GetZaxis().SetRange(ibin,ibin)
		binL = ttg_vg_jg_hist.Project3D('xy').GetMaximum()
		if binL==0: binL=1e-20
		jghist_NLL.SetBinContent(ibin,-1*ln(binL/maxL))
	ttg_vg_jg_hist.GetXaxis().SetRange(-1,-1)
	ttg_vg_jg_hist.GetYaxis().SetRange(-1,-1)
	ttg_vg_jg_hist.GetZaxis().SetRange(-1,-1)

	ttghist_NLL.SetLineColor(ROOT.kBlack)
	ttghist_NLL.GetXaxis().SetTitle('SF(t#bar{t}+#gamma)')
	ttghist_NLL.GetYaxis().SetTitle('-log(L/L_{max})')
	vghist_NLL.SetLineColor(ROOT.kBlack)
	vghist_NLL.GetXaxis().SetTitle('SF(V+#gamma)')
	vghist_NLL.GetYaxis().SetTitle('-log(L/L_{max})')
	jghist_NLL.SetLineColor(ROOT.kBlack)
	jghist_NLL.GetXaxis().SetTitle('SF(jet#rightarrow#gamma)')
	jghist_NLL.GetYaxis().SetTitle('-log(L/L_{max})')


	ttghist_NLL.Fit('pol2')
	f = ttghist_NLL.GetFunction('pol2')
#	ttghist_NLL.Draw()
	lowBound, highBound = findBounds(f)
	CMS_lumi.CMS_lumi(canvas, 2, 11)
	canvas.Print('plots/TTGamma_SF_NLL.png', '.png')
	canvas.Print('plots/TTGamma_SF_NLL.pdf', '.pdf')

	bestttgSFErr = (highBound-lowBound)/2.


	vghist_NLL.Fit('pol2')
	f = vghist_NLL.GetFunction('pol2')
#	vghist_NLL.Draw()
	lowBound, highBound = findBounds(f)
	CMS_lumi.CMS_lumi(canvas, 2, 11)
	canvas.Print('plots/Vgamma_SF_NLL.png', '.png')
	canvas.Print('plots/Vgamma_SF_NLL.pdf', '.pdf')

	bestVgSFErr = (highBound-lowBound)/2.

	jghist_NLL.Fit('pol2')
	f = jghist_NLL.GetFunction('pol2')
#	jghist_NLL.Draw()
	lowBound, highBound = findBounds(f)
	CMS_lumi.CMS_lumi(canvas, 2, 11)
	canvas.Print('plots/jet_gamma_SF_NLL.png', '.png')
	canvas.Print('plots/jet_gamma_SF_NLL.pdf', '.pdf')

	bestjgSFErr = (highBound-lowBound)/2.







	outFile = ROOT.TFile("plots/fitSF.root",'recreate')
	ttghist.Write()
	vghist.Write()
	jghist.Write()
	ttghist_NLL.Write()
	vghist_NLL.Write()
	jghist_NLL.Write()
	ttg_vg_jg_hist.Write()
	outFile.Close()

	ttgammaSig = bestttgSF*pho['TTGamma'][0]
	ttgammaSigErr = bestttgSFErr*pho['TTGamma'][0]

	ttgammaFull = bestttgSF*(pho['TTGamma'][0]+ele['TTGamma'][0]+fake['TTGamma'][0])
	ttgammaFullErr = bestttgSFErr*(pho['TTGamma'][0]+ele['TTGamma'][0]+fake['TTGamma'][0])

	print 'number of signal events', ttgammaSig, ' +/-',ttgammaSigErr
	return ttgammaSig,ttgammaSigErr, ttgammaFull, ttgammaFullErr, bestttgSF, bestttgSFErr, bestVgSF, bestVgSFErr, bestjgSF, bestjgSFErr

#calculateTTGamma()
	
