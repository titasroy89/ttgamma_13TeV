import os
from distribution_mod import distribution
import ROOT
import sys

ROOT.gROOT.SetBatch()

# initialize variables, assign values later
WJetsSF = 0.0
TopSF = 1.0
QCDSF = 0.0

import CMS_lumi

isElectron = False
isMuon = False
lep = ''

drawRatio = True
padRatio = 0.25
padOverlap = 0.15
padGap = 0.01

from Style import *

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

if len(sys.argv) > 1:
	print sys.argv
	if sys.argv[1]=='e' or 'ele' in sys.argv[1].lower():
		isElectron = True
		lep = 'ele'
	elif sys.argv[1]=='mu' or 'muon' in sys.argv[1].lower():
		isMuon = True
		lep = 'mu'
	else:
		print '#'*30
		print 'First argument must specify either electron or muon'
		print 'Allowed arguments:'
		print '   e, electron, mu, muons'
		print '#'*30
		sys.exit(1)

if isElectron and isMuon:
	print 'Error: trying to run on both electron and muon channel'
	sys.exit(1)
elif isElectron and lep =='ele':
	print 'Running on the e+jets channel'
elif isMuon and lep =='mu':
	print 'Running on the mu+jets channel'
elif lep == '':
	print 'No lepton channel specified'
	sys.exit(1)
else:
	print 'Lepton channel not properly specified'
	sys.exit(1)


systName = ''
systematicsList = ['JEC_up',  'JER_up',  'Btag_up',  'elesmear_up',  'EleEff_up',  'musmear_up',  'MuEff_up',  
		   'JEC_down','JER_down','Btag_down','elesmear_down','EleEff_down','musmear_down','MuEff_down',
		   ]

if len(sys.argv)==3:
	systName = sys.argv[2]
	if not systName in systematicsList:
		print "Unknown systematics", systName
		sys.exit(1)

#import array
#binarray = array.array('d')
#binarray.fromlist([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,300])

# load cross-sections and N gen as global variables
execfile('SF.py')


ratioPlotRanges_barrel = {'M3':0.7,
			  'photon1Et':0.7,
			  lep+'1'+lep+'2Mass':0.7,

			  }

ratioPlotRanges = {'M3':0.25,
		   'photon1Et':0.25,
		   lep+'1'+lep+'2Mass':0.7,
		       
		   }

# function definitions ####################################################

def saveTemplatesToFile(templateList, varlist, outFileName):
	outfile = ROOT.TFile(outFileName,'RECREATE')
	for template in templateList:
		for var in varlist:
			template.histList[var].SetDirectory(outfile.GetDirectory(''))
			template.histList[var].Write()
			template.histList[var].SetDirectory(0)
	outfile.Close()

def plotTemplates(dataTemplate, MCTemplateList, SignalTemplateZoomList, varlist, outDirName):
	ratioRanges = ratioPlotRanges
	if 'barrel' in outDirName:
		ratioRanges = ratioPlotRanges_barrel

	ROOT.TGaxis.SetMaxDigits(3)	

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
		
	latex = ROOT.TLatex()
	latex.SetNDC()
	latex.SetTextAlign(12)
	latex.SetTextSize(0.037)
	latex.SetLineWidth(2)
	
	canvasRatio = ROOT.TCanvas('c1Ratio','c1Ratio',W,H)
	# references for T, B, L, R
	T = 0.08*H
	B = 0.12*H 
	L = 0.12*W
	R = 0.04*W
	canvasRatio.SetFillColor(0)
	canvasRatio.SetBorderMode(0)
	canvasRatio.SetFrameFillStyle(0)
	canvasRatio.SetFrameBorderMode(0)
	canvasRatio.SetLeftMargin( L/W )
	canvasRatio.SetRightMargin( R/W )
	canvasRatio.SetTopMargin( T/H )
	canvasRatio.SetBottomMargin( B/H )
	canvasRatio.SetTickx(0)
	canvasRatio.SetTicky(0)
	pad1 = ROOT.TPad("p1","p1",0,padRatio-padOverlap,1,1)
	pad2 = ROOT.TPad("p2","p2",0,0,1,padRatio+padOverlap)
        pad1 = ROOT.TPad("p1","p1",0,max(padRatio-padOverlap+padGap/2,0),1,1)
        pad2 = ROOT.TPad("p2","p2",0,0,1,padRatio+padOverlap-padGap/2)


	pad1.SetLeftMargin( L/W )
	pad1.SetRightMargin( R/W )
	pad1.SetTopMargin( T/H/(1-padRatio) )
	pad1.SetBottomMargin( (padOverlap+padGap)/(1-padRatio+padOverlap) )
	pad2.SetLeftMargin( L/W )
	pad2.SetRightMargin( R/W )
	pad2.SetTopMargin( (padOverlap)/(padRatio+padOverlap) )
	pad2.SetBottomMargin( B/H/(padRatio+padOverlap) )

	pad1.SetFillColor(0)
	pad1.SetBorderMode(0)
	pad1.SetFrameFillStyle(0)
	pad1.SetFrameBorderMode(0)
	pad1.SetTickx(0)
	pad1.SetTicky(0)

	pad2.SetFillColor(0)
	pad2.SetFillStyle(4000)
	pad2.SetBorderMode(0)
	pad2.SetFrameFillStyle(0)
	pad2.SetFrameBorderMode(0)
	pad2.SetTickx(0)
	pad2.SetTicky(0)


	ROOT.SetOwnership(canvas, False)
	ROOT.SetOwnership(canvasRatio, False)
	ROOT.SetOwnership(pad1, False)
	ROOT.SetOwnership(pad2, False)
	
	canvasRatio.cd()
	pad1.Draw()
	pad2.Draw()

	
	for var in varlist:
		canvas.cd()
		legend = ROOT.TLegend(0.71, 0.99 - canvas.GetTopMargin() - 0.05*(1 + len(MCTemplateList) + len(SignalTemplateZoomList)), 0.99-canvas.GetRightMargin(), 0.99-canvas.GetTopMargin())
		legend.SetBorderSize(0)
		legend.SetFillColor(ROOT.kWhite)
		
		legendR = ROOT.TLegend(0.71, 1. - 0.1/(1.-padRatio) - 0.05/(1.-padRatio)*(len(MCTemplateList) + len(SignalTemplateZoomList)), 0.94, 1-0.1/(1.-padRatio))
		legendR.SetBorderSize(0)
		legendR.SetFillColor(0)

		if dataTemplate is not None:
			legend.AddEntry(dataTemplate.histList[var], dataTemplate.legName, 'pl')
			legendR.AddEntry(dataTemplate.histList[var], dataTemplate.legName, 'pl')
		
		# MC templates listed in the order they appear in legend
		for mc in MCTemplateList[::-1]:
			mcHist = mc.histList[var]
			legend.AddEntry(mcHist, mc.legName, 'f')
			legendR.AddEntry(mcHist, mc.legName, 'f')
		
		stack = ROOT.THStack('stack_'+var,var)
		# reverse order for stack to be consistent with legend

		for mc in MCTemplateList:
			mcHist = mc.histList[var]
			#if var == 'M3pho':
			#	mcHist.Rebin(2)
			stack.Add(mcHist)

		if dataTemplate is not None:
			#if var == 'M3pho':
			#	dataTemplate.histList[var].Rebin(2)
			if dataTemplate.histList[var].GetMaximum() > stack.GetMaximum():
				stack.SetMaximum(dataTemplate.histList[var].GetMaximum())

		if 'cut_flow' in var: # or 'MET' in var:
			canvas.SetLogy(1)
			stack.SetMinimum(100)
		else:
			canvas.SetLogy(0)
		
		stack.Draw('HIST')
		
		if 'barrel' in outDirName and 'photon1SigmaIEtaIEta' in var:
			stack.GetXaxis().SetRangeUser(0.0,0.025)
		
		if dataTemplate is not None:
			stack.GetXaxis().SetTitle(dataTemplate.histList[var].GetXaxis().GetTitle())
			stack.GetYaxis().SetTitle(dataTemplate.histList[var].GetYaxis().GetTitle())
			if 'mu1mu2Mass' in var:
				stack.GetXaxis().SetTitle("M(#mu,#mu) (GeV)")
				stack.GetXaxis().SetRangeUser(20,180)
				dataTemplate.histList[var].GetXaxis().SetRangeUser(20,180)
			if 'ele1ele2Mass' in var:
				stack.GetXaxis().SetTitle("M(e,e) (GeV)")
				stack.GetXaxis().SetRangeUser(20,180)
				dataTemplate.histList[var].GetXaxis().SetRangeUser(20,180)
				
#			stack.GetYaxis().SetTitleOffset(1.8)
		stack.SetTitle('')

		if dataTemplate is not None:
			dataTemplate.histList[var].SetMarkerStyle(20)
			dataTemplate.histList[var].SetMarkerSize(1.2)
			dataTemplate.histList[var].SetMarkerColor(ROOT.kBlack)
			dataTemplate.histList[var].SetLineColor(ROOT.kBlack)
			dataTemplate.histList[var].Draw('ESAME')
			dataTemplate.histList[var].Draw('ESAME')
			
		for signal,zoom in SignalTemplateZoomList:
			sigHist = signal.histList[var].Clone()
			#sigHist.SetFillStyle(3244)
			sigHist.Scale(zoom)
			sigHist.Draw('HISTSAME')
			if zoom != 1:
				legend.AddEntry(sigHist, signal.legName + ' x ' + str(zoom), 'f')
			else:
				legend.AddEntry(sigHist, signal.legName, 'f')
		if 'cut_flow' not in var:
			legend.Draw()

		channelText = ""
		if isMuon: channelText = "#mu#mu"
		if isElectron: channelText = "ee"

		CMS_lumi.channelText = channelText
		CMS_lumi.writeExtraText = True
		CMS_lumi.writeChannelText = True
		
		CMS_lumi.CMS_lumi(canvas, 2, 11)
		canvas.Update();
		canvas.RedrawAxis();
		canvas.Print(outDirName+'/'+var+".pdf",".pdf");
		canvas.Print(outDirName+'/'+var+".png",".png");
		
		

		####RATIO PLOT
		if drawRatio:
			ratio = dataTemplate.histList[var].Clone("temp")
			ratio.Divide(stack.GetStack().Last())
			# ratio = stack.GetStack().Last()
        		# ratio.Divide(dataTemplate.histList[var])
			
        		canvasRatio.cd()
        		pad1.cd()
        
        		stack.Draw('HIST')
        		if 'ele1MVA' in var:
        			stack.GetXaxis().SetRangeUser(-1.1,0.0)
        
        		if 'barrel' in outDirName and 'photon1SigmaIEtaIEta' in var:
        			stack.GetXaxis().SetRangeUser(0.0,0.025)
        		pad1.Update()
			y2 = pad1.GetY2()
			print y2
			stack.SetMinimum(-0.02*y2)
			pad1.Update()
			
			pad1.Update()
        		if dataTemplate is not None:
        			stack.GetXaxis().SetTitle('')
        			stack.GetYaxis().SetTitle(dataTemplate.histList[var].GetYaxis().GetTitle())
        		stack.SetTitle('')
			stack.GetXaxis().SetLabelSize(0)
			stack.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()/(1.-padRatio+padOverlap))
			stack.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()/(1.-padRatio+padOverlap))
			print stack.GetYaxis().GetTitleOffset()
			stack.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleYOffset()*(1.-padRatio+padOverlap))
        
        		if dataTemplate is not None:
        			dataTemplate.histList[var].Draw('ESAME')
        					
        		for signal,zoom in SignalTemplateZoomList:
        			sigHist = signal.histList[var].Clone()
        			#sigHist.SetFillStyle(3244)
        			sigHist.Scale(zoom)
        			sigHist.Draw('HISTSAME')
        			if zoom != 1:
        				legendR.AddEntry(sigHist, signal.legName + ' x ' + str(zoom), 'f')
        			else:
        				legendR.AddEntry(sigHist, signal.legName, 'f')
        		if 'cut_flow' not in var:
        			legendR.Draw()
        
        		if dataTemplate is not None:
        			ratio.GetXaxis().SetTitle(dataTemplate.histList[var].GetXaxis().GetTitle())
        			ratio.GetYaxis().SetTitle('Data/Theory')
				if 'mu1mu2Mass' in var:
					ratio.GetXaxis().SetTitle("M(#mu,#mu) (GeV)")
					print var
				if 'ele1ele2Mass' in var:
					ratio.GetXaxis().SetTitle("M(e,e) (GeV)")
				ratio.GetYaxis().CenterTitle()
        		ratio.SetTitle('')
			ratio.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()/(padRatio+padOverlap))
			ratio.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()/(padRatio+padOverlap))
			ratio.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()/(padRatio+padOverlap))
			ratio.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()/(padRatio+padOverlap))
			# ratio.GetXaxis().SetTitleOffset(ROOT.gStyle.GetTitleXOffset()*(1-padRatio))
			# ratio.GetXaxis().SetTitleOffset(ROOT.gStyle.GetTitleXOffset()*(1-padRatio))
			ratio.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleYOffset()*(padRatio+padOverlap))

			ratio.GetYaxis().SetRangeUser(0.75,1.25)
			if var in ratioRanges:
				span = ratioRanges[var]
				ratio.GetYaxis().SetRangeUser(1-span, 1+span)
			ratio.GetYaxis().SetNdivisions(503)
				
        
        		pad2.cd()
        		ratio.SetMarkerStyle(2)		
        		ratio.SetLineColor(ROOT.kBlack)
			ratio.SetLineWidth(1)
			oneLine = ROOT.TF1("oneLine","1",-999,9999)
#			oneLine = ROOT.TLine(ratio.GetXaxis().GetXmin(),1.,ratio.GetXaxis().GetXmax(),1.)
			oneLine.SetLineColor(ROOT.kBlack)
			oneLine.SetLineWidth(1)
			oneLine.SetLineStyle(2)

        		ratio.Draw()        		
			oneLine.Draw("same")

			pad2.Update()
        		CMS_lumi.CMS_lumi(canvasRatio, 2, 11)
        
			canvasRatio.Update();
			canvasRatio.RedrawAxis();
			
			canvasRatio.Print(outDirName+'/'+var+"_ratio.pdf",".pdf");
			canvasRatio.Print(outDirName+'/'+var+"_ratio.png",".png");

		
def loadDataTemplate(varlist, inputDir, prefix):
	templPrefix = inputDir+prefix
	DataTempl = distribution('Data', 'Data', [
		(templPrefix+'Data_a.root', 1),
		(templPrefix+'Data_b.root', 1),
		(templPrefix+'Data_c.root', 1),
		(templPrefix+'Data_d.root', 1),
		], varlist)
	return DataTempl


def loadMCTemplates(varList, inputDir, prefix, titleSuffix, fillStyle):
	templPrefix = inputDir+prefix
	
	MCtemplates = {}
	
	#MCtemplates['mst_510_200'] = distribution('mst_510_200'+titleSuffix,[
	#		(templPrefix+'mst_510_M3_5050_M1_200.root', gSF*0.0751004/15000)
	#	],varList, 620,3244)
	
	#MCtemplates['WHIZARD'] = distribution('TTGamma'+titleSuffix, [
	#	(templPrefix+'WHIZARD.root', TopSF*gSF*TTgamma_xs/WHIZARD_num)
	#	], varList, 98, fillStyle)
	
	MCtemplates['WHIZARD'] = distribution('TTGamma'+titleSuffix, 't#bar{t}+#gamma',[
		(templPrefix+'TTGamma.root', TopSF*gSF*newTTgamma_xs/newTTgamma_num)
		], varList, ROOT.kRed +1, fillStyle)
	
	MCtemplates['TTJets'] = distribution('TTJets'+titleSuffix, 't#bar{t}+jets',[
		(templPrefix+'TTJets1l.root', TopSF*gSF*TTJets1l_xs/TTJets1l_num),
		(templPrefix+'TTJets2l.root', TopSF*gSF*TTJets2l_xs/TTJets2l_num),
		(templPrefix+'TTJetsHad.root', TopSF*gSF*TTJetsHad_xs/TTJetsHad_num),
		], varList ,ROOT.kRed-7, fillStyle)
	
	###################################
	#return MCtemplates
	###################################
	nonWJetsSF = 1.0
	#nonWJetsSF = WJetsSF
	
	MCtemplates['Wgamma'] = distribution('Wgamma'+titleSuffix, 'W+#gamma', [
        (templPrefix+'Wgamma.root', nonWJetsSF*gSF*Wgamma_xs/Wgamma_num),
    #    (templPrefix+'WWgamma.root', gSF*WWgamma_xs/WWgamma_num),
        ], varList, ROOT.kGray, fillStyle)

	MCtemplates['Zgamma'] = distribution('Zgamma'+titleSuffix, 'Z+#gamma', [
        (templPrefix+'Zgamma.root', nonWJetsSF*gSF*Zgamma_xs/Zgamma_num),
    #    (templPrefix+'WWgamma.root', gSF*WWgamma_xs/WWgamma_num),
        ], varList, ROOT.kAzure+3, fillStyle)

	MCtemplates['SingleTop'] = distribution('SingleTop'+titleSuffix, 'Single Top', [
		(templPrefix+'SingleT_t.root',      nonWJetsSF*gSF*SingTopT_xs/SingTopT_num),
        (templPrefix+'SingleT_s.root',      nonWJetsSF*gSF*SingTopS_xs/SingTopS_num),
        (templPrefix+'SingleT_tw.root',     nonWJetsSF*gSF*SingToptW_xs/SingToptW_num),
        (templPrefix+'SingleTbar_t.root',   nonWJetsSF*gSF*SingTopbarT_xs/SingTopbarT_num),
        (templPrefix+'SingleTbar_s.root',   nonWJetsSF*gSF*SingTopbarS_xs/SingTopbarS_num),
        (templPrefix+'SingleTbar_tw.root',  nonWJetsSF*gSF*SingTopbartW_xs/SingTopbartW_num),
		], varList, ROOT.kMagenta, fillStyle)
	
	MCtemplates['WJets'] = distribution('WJets'+titleSuffix, 'W+jets', [
        #(templPrefix+'WJets.root', WJetsSF*gSF*WJets_xs/WJets_num),
		(templPrefix+'W3Jets.root', WJetsSF*gSF*W3Jets_xs/W3Jets_num),
		(templPrefix+'W4Jets.root', WJetsSF*gSF*W4Jets_xs/W4Jets_num),
		], varList, ROOT.kGreen-3, fillStyle)
		
	MCtemplates['ZJets'] = distribution('ZJets'+titleSuffix, 'Z+jets', [
		(templPrefix+'ZJets.root', nonWJetsSF*gSF*ZJets_xs/ZJets_num)], varList, ROOT.kAzure-2, fillStyle)

	
	# MCtemplates['Other'] = distribution('Diboson'+titleSuffix, [

        # (templPrefix+'WZ_3lnu.root', nonWJetsSF*gSF*WZ_3lnu_xs/WZ_3lnu_num),
        # (templPrefix+'WZ_2l2q.root', nonWJetsSF*gSF*WZ_2l2q_xs/WZ_2l2q_num),
        
        # (templPrefix+'ZZ_2e2mu.root', nonWJetsSF*gSF*ZZ_2e2mu_xs/ZZ_2e2mu_num),
        # (templPrefix+'ZZ_2e2tau.root', nonWJetsSF*gSF*ZZ_2e2tau_xs/ZZ_2e2tau_num),
        # (templPrefix+'ZZ_2mu2tau.root', nonWJetsSF*gSF*ZZ_2mu2tau_xs/ZZ_2mu2tau_num),
        # (templPrefix+'ZZ_4e.root', nonWJetsSF*gSF*ZZ_4e_xs/ZZ_4e_num),
        # (templPrefix+'ZZ_4mu.root', nonWJetsSF*gSF*ZZ_4mu_xs/ZZ_4mu_num),
        # (templPrefix+'ZZ_4tau.root', nonWJetsSF*gSF*ZZ_4tau_xs/ZZ_4tau_num),
        
        # (templPrefix+'WW_2l2nu.root', nonWJetsSF*gSF*WW_2l2nu_xs/WW_2l2nu_num),

        #   #(templPrefix+'TTW.root', gSF*TTW_xs/TTW_num),
        #   #(templPrefix+'TTZ.root', gSF*TTZ_xs/TTZ_num),
	# 	], varList, 49, fillStyle)

	return MCtemplates


def makeAllPlots(varList, inputDir, dataDir, outDirName,savePlots = True):
	# load templates PreSel	
	DataTempl = loadDataTemplate(varList, dataDir, 'hist_1pho_top_')

	MCTemplDict = loadMCTemplates(varList, inputDir, 'hist_1pho_top_','',1001)
	MCTempl = []
	MCTempl.append(MCTemplDict['WHIZARD'])
	MCTempl.append(MCTemplDict['TTJets'])
	MCTempl.append(MCTemplDict['SingleTop'])
	MCTempl.append(MCTemplDict['WJets'])
	MCTempl.append(MCTemplDict['Zgamma'])
	MCTempl.append(MCTemplDict['Wgamma'])
	MCTempl.append(MCTemplDict['ZJets'])
#	MCTempl.append(MCTemplDict['Other'])

	saveTemplatesToFile([DataTempl] + MCTempl, ['MET',lep+'1'+lep+'2Mass'], outDirName+'/templates_presel.root')
	
	if savePlots: plotTemplates( DataTempl, MCTempl, [], varList, outDirName+'/presel')
	


varList_all = ['nVtx',
			'MET','Ht','WtransMass','M3','M3first','minM3','M3pho','dRpho3j','M3phoMulti', 
			#'M3_0_30', 'M3_30_100', 'M3_100_200', 'M3_200_300', 'M3_300_up', #'M3minPt',
			lep+'1Pt',lep+'1Eta',lep+'1RelIso',
#			lep+'1D0',lep+'1MVA',lep+'1Dz',
			lep+'2Pt',lep+'2RelIso',
			lep+'1'+lep+'2Mass',
			# lep+'1sigmaIetaIeta',lep+'1EoverP',
			# lep+'1DrJet',lep+'1pho1Mass',
			# 'looseEleDrGenPho',
			# 'cut_flow',
			# 'genPhoRegionWeight',
			# 'nJets',
			# 'jet1Pt','jet2Pt','jet3Pt','jet4Pt','jet1Eta','jet2Eta','jet3Eta','jet4Eta',
			# 'photon1Et','photon1Eta','photon1HoverE','photon1SigmaIEtaIEta',
			# 'photon1DrElectron','photon1DrJet',
			# 'photon1ChHadIso','photon1NeuHadIso','photon1PhoIso',
			# 'photon1ChHadSCRIso','photon1PhoSCRIso',
			# 'photon1ChHadRandIso','photon1PhoRandIso',
			# 'photon1MotherID','photon1GMotherID','photon1DrMCbquark','GenPhotonEt',
			#'photon1_Sigma_ChSCRIso'
			]
# main part ##############################################################################################



if isElectron:
	if systName == '':
		InputHist = '/uscms_data/d2/dnoonan/Electrons/TwoEleHists/hist_bins_twoEle/'
		DataHist = '/uscms_data/d2/dnoonan/Electrons/TwoEleHists/hist_bins_twoEle/'
		# InputHist = '/eos/uscms/store/user/dnoonan/EleHists_looseVeto_Oct27/hist_bins_twoEle/'
		# DataHist = '/eos/uscms/store/user/dnoonan/EleHists_looseVeto_Oct27/hist_bins_twoEle/'
		
		makeAllPlots(varList_all, InputHist, DataHist, 'di_ele_cross_check/plots')
		
		InputHist = '/eos/uscms/store/user/dnoonan/EleHists_looseVeto_Oct27/hist_bins_twoEle_zeroB/'
		DataHist = '/eos/uscms/store/user/dnoonan/EleHists_looseVeto_Oct27/hist_bins_twoEle_zeroB/'
		
		makeAllPlots(varList_all, InputHist, DataHist, 'di_ele_cross_check_zeroB/plots')
	else:
		InputHist = '/uscms_data/d2/dnoonan/Electrons/TwoEleHists/hist_bins_twoEle_%s/'%systName
		DataHist = '/uscms_data/d2/dnoonan/Electrons/TwoEleHists/hist_bins_twoEle/'

		if not os.path.exists('di_ele_cross_check/%s'%systName):
			os.mkdir('di_ele_cross_check/%s'%systName)
			os.mkdir('di_ele_cross_check/%s/plots'%systName)
			os.mkdir('di_ele_cross_check/%s/plots/presel'%systName)
		makeAllPlots(varList_all, InputHist, DataHist, 'di_ele_cross_check/%s/plots'%systName, savePlots = False)

		# makeAllPlots(varList_all, InputHist, DataHist, 'di_ele_cross_check/plots')
		

if isMuon:
	if systName == '':
		# InputHist = '/eos/uscms/store/user/dnoonan/MuHists_looseVeto_Oct27/hist_bins_twoMu/'
		# DataHist = '/eos/uscms/store/user/dnoonan/MuHists_looseVeto_Oct27/hist_bins_twoMu/'
		InputHist = '/uscms_data/d2/dnoonan/TTGammaElectrons/TwoMuHists/hist_bins_twoMu/'
		DataHist = '/uscms_data/d2/dnoonan/TTGammaElectrons/TwoMuHists/hist_bins_twoMu/'

		makeAllPlots(varList_all, InputHist, DataHist, 'di_mu_cross_check/plots')
	else:
		InputHist = '/uscms_data/d2/dnoonan/TTGammaElectrons/TwoMuHists/hist_bins_twoMu_%s/'%systName
		DataHist = '/uscms_data/d2/dnoonan/TTGammaElectrons/TwoMuHists/hist_bins_twoMu/'
#		DataHist = '/eos/uscms/store/user/dnoonan/MuHists_looseVeto_Oct27/hist_bins_twoMu/'
		
		if not os.path.exists('di_mu_cross_check/%s'%systName):
			os.mkdir('di_mu_cross_check/%s'%systName)
			os.mkdir('di_mu_cross_check/%s/plots'%systName)
			os.mkdir('di_mu_cross_check/%s/plots/presel'%systName)

		makeAllPlots(varList_all, InputHist, DataHist, 'di_mu_cross_check/%s/plots'%systName, savePlots=False)

