from distribution_mod import distribution
import ROOT
import sys, os

import templateFits



drawRatio = True
padRatio = 0.25
padOverlap = 0.05
padGap = 0.01

ROOT.gROOT.SetBatch()
#########
#Style
from Style import *

#import CMS_lumi
 
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
isElectron = False
isMuon = False
lep = ''

#Add a couple of flags for stopping at various points in the code
skipPhoton = False #stops before the photon fitting
skipAfterMET = False #stops before the photon fitting
skipAfterM3 = False #stops before the photon fitting
skipCalc = False #stops before the calc_the_answer step
skipMET = False #stops before the MET fitting, just does the photon purity

skipQCDphoton = False

SaveOutput = False

 ######## Add an argument to determine if running on electrons or muons ######## 
isSyst = False
systematic = ''
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

else:
	print '#'*30
	print 'At least one argument is required,'
	print 'Must specify if begin run on electrons (e or electron) or muons (mu or muon)'
	print '#'*30
	sys.exit(1)

saveSystName = 'nominal'
if isSyst: saveSystName = systematic

if SaveOutput:
	if isElectron: folder = 'ratioFiles_ele/'
	if isMuon: folder = 'ratioFiles_mu/'
	sys.stdout = open(folder+'ratio_'+saveSystName+'.txt','w')
	# if isSyst:
	# 	sys.stdout = open(folder+'ratio_'+systematic+'.txt','w')
	# else:
	# 	sys.stdout = open(folder+'ratio_nominal.txt','w')


 ######## Error checking that a lepton channel was selected'
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


# load cross-sections and N gen as global variables
execfile('SF.py')


ratioPlotRanges_barrel = {'M3':0.7,
			  'photon1Et':0.7,
			  
			  }	

ratioPlotRanges = {'M3':0.25,
		   'photon1Et':0.25,
		   
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
	canvas = ROOT.TCanvas('c1','c1',800,800)
	
	latex = ROOT.TLatex()
	latex.SetNDC()
	latex.SetTextAlign(12)
	latex.SetTextSize(0.03)
	latex.SetLineWidth(2)
	
	for var in varlist:
		legend = ROOT.TLegend(0.7, 1 - 0.05*(1 + len(MCTemplateList) + len(SignalTemplateZoomList)), 0.99, 1.00)
		legend.SetBorderSize(0)
		legend.SetFillColor(10)
		
		if dataTemplate is not None:
			legend.AddEntry(dataTemplate.histList[var], dataTemplate.name, 'pl')
		
		# MC templates listed in the order they appear in legend
		for mc in MCTemplateList:
			mcHist = mc.histList[var]
			legend.AddEntry(mcHist, mc.name, 'f')
		
		stack = ROOT.THStack('stack_'+var,var)
		# reverse order for stack to be consistent with legend
		MCTemplateList.reverse()
		for mc in MCTemplateList:
			mcHist = mc.histList[var]
			#if var == 'M3pho':
			#	mcHist.Rebin(2)
			stack.Add(mcHist)
		MCTemplateList.reverse()

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
		stack.SetTitle('')

		if dataTemplate is not None:
			dataTemplate.histList[var].Draw('ESAME')
					
		for signal,zoom in SignalTemplateZoomList:
			sigHist = signal.histList[var].Clone()
			#sigHist.SetFillStyle(3244)
			sigHist.Scale(zoom)
			sigHist.Draw('HISTSAME')
			if zoom != 1:
				legend.AddEntry(sigHist, signal.name + ' x ' + str(zoom), 'f')
			else:
				legend.AddEntry(sigHist, signal.name, 'f')
		if 'cut_flow' not in var:
			legend.Draw()
		
		#latex.DrawLatex(0.1,0.94,'CMS Preliminary #sqrt{s} = 8 TeV')
		if not isSyst:
			canvas.SaveAs(outDirName+'/'+var+'.pdf')
def loadDataTemplate(varlist, inputDir, prefix):
	templPrefix = inputDir+prefix
	DataTempl = distribution('data', 'data', [
		(templPrefix+'data_c.root', 1),
		(templPrefix+'data_d.root', 1),
		(templPrefix+'data_d_PR.root',1),
		], varlist)
	return DataTempl

def loadQCDTemplate(varlist, inputDir, prefix):
	templPrefix = inputDir+prefix
	QCD_sf = QCDSF
	QCDTempl = distribution('QCD', 'QCD', [
	#	(templPrefix+'data_c.root', QCD_sf),
		(templPrefix+'data_d.root', QCD_sf),
	#	(templPrefix+'data_d_PR.root', QCD_sf),
		(templPrefix+'ttjets.root', -1 * QCD_sf * TopSF * gSF * TTJets_xs/TTJets_num),
	
		######## Added in all other channels of MC, previously just ttjets 1l and 2l removed to get QCD template ######## 
		(templPrefix+'ttjets.root', -1 * gSF * TTJets_xs/TTJets_num),
		(templPrefix+'ttgamma.root', -1 * gSF * newTTgamma_xs/newTTgamma_num),
		(templPrefix+'st_t.root',     -1 * gSF * SingTopT_xs/SingTopT_num),
		(templPrefix+'st_t_anti.root',     -1 * gSF * SingTopTbar_xs/SingTopTbar_num),
		 (templPrefix+'st_s.root',     -1 * gSF * SingTopS_xs/SingTopS_num),
		(templPrefix+'st_tW.root',    -1 * gSF * SingToptW_xs/SingToptW_num),
		(templPrefix+'st_tW_anti.root',  -1 * gSF * SingTopbartW_xs/SingTopbartW_num),
		(templPrefix+'Wjets.root', -1 * gSF * WJets_xs/WJets_num),
		(templPrefix+'DYJets.root',  -1 * gSF * ZJets_xs/ZJets_num),

	], varlist, ROOT.kYellow)
	return QCDTempl
        
def loadMCTemplates(varList, inputDir, prefix, titleSuffix, fillStyle):
	templPrefix = inputDir+prefix
	
	MCtemplates = {}
	
	
	MCtemplates['TTGamma'] = distribution('TTGamma'+titleSuffix, 't#bar{t}+#gamma', [
		(templPrefix+'ttgamma.root', gSF*newTTgamma_xs/newTTgamma_num)
		], varList, ROOT.kRed +1, fillStyle)
	
	MCtemplates['TTJets'] = distribution('TTJets'+titleSuffix, 't#bar{t}+jets', [
		(templPrefix+'ttjets.root', gSF*TTJets_xs/TTJets_num),
		], varList ,ROOT.kRed -7, fillStyle)
	SF_ttjets =  gSF*TTJets_xs/TTJets_num
	print SF_ttjets



	MCtemplates['SingleTop'] = distribution('SingleTop'+titleSuffix, 'Single Top', [
		(templPrefix+'st_t.root',  gSF*SingTopT_xs/SingTopT_num),
		(templPrefix+'st_t_anti.root',    gSF*SingTopTbar_xs/SingTopTbar_num),
		(templPrefix+'st_tW.root',     gSF*SingToptW_xs/SingToptW_num),
		(templPrefix+'st_tW_anti.root',  gSF*SingTopbartW_xs/SingTopbartW_num),
		(templPrefix+'st_s.root',    gSF*SingTopS_xs/SingTopS_num),
		], varList, ROOT.kMagenta, fillStyle)
	MCtemplates['WJets'] = distribution('WJets'+titleSuffix, 'W+jets', [
		(templPrefix+'Wjets.root', gSF*WJets_xs/WJets_num),
		], varList, ROOT.kGreen -3, fillStyle)

	MCtemplates['ZJets'] = distribution('ZJets'+titleSuffix, 'Z+jets', [
		(templPrefix+'DYJets.root', gSF*ZJets_xs/ZJets_num)], varList, ROOT.kAzure-2, fillStyle)
	return MCtemplates


def saveNoMETTemplates(inputDir, inputData, outFileName, histName):
	varList = ['MET','M3',lep+'1RelIso','ele1MVA']
	DataTempl = loadDataTemplate(varList, inputData, histName)
	MCTemplDict = loadMCTemplates(varList, inputDir, histName,'',1001)
	MCTempl = []
	MCTempl.append(MCTemplDict['TTGamma'])
	MCTempl.append(MCTemplDict['TTJets'])
	MCTempl.append(MCTemplDict['SingleTop'])
	MCTempl.append(MCTemplDict['WJets'])
	MCTempl.append(MCTemplDict['ZJets'])
	print inputDir
	print histName
	saveTemplatesToFile([DataTempl] + MCTempl, varList, outFileName)


def savePreselTemplates(inputDir, inputData, outFileName):
	if WJetsSF != 1.0 or TopSF != 1.0:
		print 'We want to save templates for M3 fit, but the SFs are not 1.0'
		print 'exiting'
		return
	
	varList = ['MET','M3',]
	DataTempl = loadDataTemplate(varList, inputData, 'hist_1pho_top_')#change
	
	MCTemplDict = loadMCTemplates(varList, inputDir, 'hist_1pho_top_','',1001) #change
	MCTempl = []
	MCTempl.append(MCTemplDict['TTGamma'])
	MCTempl.append(MCTemplDict['TTJets'])
	MCTempl.append(MCTemplDict['SingleTop'])
	MCTempl.append(MCTemplDict['WJets'])
	MCTempl.append(MCTemplDict['ZJets'])
	if QCDSF > 0.0001:
		MCTempl.append(QCDTempl)
	saveTemplatesToFile([DataTempl] + MCTempl, varList, outFileName)


def makeAllPlots(varList, inputDir, dataDir, outDirName):
	# load templates PreSel	
	DataTempl = loadDataTemplate(varList, dataDir, 'hist_1pho_top_') 
	MCTemplDict = loadMCTemplates(varList, inputDir, 'hist_1pho_top_','',1001) 
	MCTempl = []
 	MCTempl.append(MCTemplDict['TTGamma'])
	MCTempl.append(MCTemplDict['TTJets'])
	MCTempl.append(MCTemplDict['SingleTop'])
	MCTempl.append(MCTemplDict['WJets'])
	MCTempl.append(MCTemplDict['ZJets'])
	
	saveTemplatesToFile([DataTempl] + MCTempl,['MET','M3','WtransMass','genPhoRegionWeight','MCcategory'], 'templates_presel_scaled.root')
	plotTemplates( DataTempl, MCTempl, [], varList, outDirName+'/presel')
	
	del MCTempl
	del MCTemplDict
	
	shortVarList = varList[:]
	shortVarList.remove('cut_flow')
	shortVarList.remove('genPhoRegionWeight')
	shortVarList.append('ele1pho1Mass')

	print "HERE"
	
	############################
	return
	############################


	
	
	############################
	############################


varList_all = [ 'nVtx','MET','Ht','WtransMass','M3', 
			#'M3_0_30', 'M3_30_100', 'M3_100_200', 'M3_200_300', 'M3_300_up', #'M3minPt',
			lep+'1Pt',lep+'1Eta',lep+'1RelIso',
			'genPhoRegionWeight', 'MCcategory',
			'cut_flow',
			'nJets',
			'jet1Pt','jet2Pt','jet3Pt','jet4Pt','jet1Eta','jet2Eta','jet3Eta','jet4Eta',
			'photon1Et','photon1Eta','photon1HoverE','photon1SigmaIEtaIEta',
			'photon1DrElectron','photon1DrJet',
			'photon1ChHadIso','photon1NeuHadIso','photon1PhoIso',
			'photon1ChHadSCRIso','photon1PhoSCRIso',
			'photon1ChHadRandIso','photon1PhoRandIso',
			'photon1MotherID','photon1GMotherID','photon1DrMCbquark','GenPhotonEt',
			#'photon1_Sigma_ChSCRIso'
			]
# main part ##############################################################################################
if systematic in ['Btag_down','Btag_up','EleEff_down','EleEff_up','JEC_down','JEC_up','JER_down','JER_up','PU_down','PU_up','elesmear_down','elesmear_up','pho_down','pho_up','toppt_down','toppt_up','MuEff_down','MuEff_up','musmear_down','musmear_up']:
	outSuffix = '_'+systematic
else:
	outSuffix = ''

if isElectron:
	InputHist = '/eos/uscms/store/user/dnoonan/EleHists_looseVeto_Nov22/hist_bins'+outSuffix+'/'
	DataHist =  '/eos/uscms/store/user/dnoonan/EleHists_looseVeto_Nov22/hist_bins/'
if isMuon:

	InputHist = '/uscms_data/d3/troy2012/ttgamma_13TeV/CMSSW_7_4_14/src/TTGammaSemiLep/nominal_noMET_2Btag'+outSuffix+'/'
	DataHist =  '/uscms_data/d3/troy2012/ttgamma_13TeV/CMSSW_7_4_14/src/TTGammaSemiLep/nominal_noMET_2Btag/'

######## Added in a printout of histogram locations, for easier tracking later on ######## 

print 'Input Histogram location:', InputHist
print 'Data Histogram location:', DataHist



makeAllPlots(varList_all, InputHist, DataHist, 'plots')




if skipMET:
	print '*'*80
	print 'Stopping code before the MET fit'
	exit()


