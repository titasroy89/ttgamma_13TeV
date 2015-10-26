from distribution_mod import distribution
import ROOT
import sys, os

import templateFits
import qcd_fit
import calc_the_answer_test
import vgamma_fit

import mcEventsTable

mcEventsTable.egammaSF      = 1.458
mcEventsTable.egammaSFerr      = 0.1985

ROOT.gROOT.SetBatch()
#########
#Style
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

#Flags to get which lepton channel is being used
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
	if 'skipMET' in sys.argv:
		skipMET = True
		sys.argv.remove('skipMET')
	if 'skipAfterMET' in sys.argv:
		skipAfterMET = True
		sys.argv.remove('skipAfterMET')
	if 'skipAfterM3' in sys.argv:
		skipAfterM3 = True
		sys.argv.remove('skipAfterM3')
	if 'skipPhoton' in sys.argv:
		skipPhoton = True
		sys.argv.remove('skipPhoton')
	if 'skipCalc' in sys.argv:
		skipCalc = True
		sys.argv.remove('skipCalc')
	if 'skipQCDphoton' in sys.argv:
		skipQCDphoton = True
		sys.argv.remove('skipQCDphoton')
	if 'SaveOutput' in sys.argv:
		SaveOutput = True
		sys.argv.remove('SaveOutput')
	if len(sys.argv) > 2:
		systematic = sys.argv[2]
		if systematic == 'zeroB':
			print 'zeroB'
		else:
			isSyst = True
#			sys.stdout = open('ratio_'+systematic+'.txt','w')

else:
	print '#'*30
	print 'At least one argument is required,'
	print 'Must specify if begin run on electrons (e or electron) or muons (mu or muon)'
	print '#'*30
	sys.exit(1)

if SaveOutput:
	if isElectron: folder = 'ratioFiles_ele/'
	if isMuon: folder = 'ratioFiles_mu/'
	if isSyst:
		sys.stdout = open(folder+'ratio_'+systematic+'.txt','w')
	else:
		sys.stdout = open(folder+'ratio_nominal.txt','w')


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

# initialize variables, assign values later
TopSF = 1.0
WJetsSF = 1.0
QCDSF = 1.0
ZJetsSF = 1.0
otherMCSF = 1.0
WgammaSF = 1.0

TopSFErr = 0.0
WJetsSFErr = 0.0
QCDSFErr = 0.0
ZJetsSFErr = 0.0
otherMCSFErr = 0.0
WgammaSFErr = 0.0

if isElectron:	
	ZJetsSF = 1.26  
	ZJetsSFErr = 0.06
if isMuon:
	ZJetsSF = 1.14 
	ZJetsSFErr = 0.06

if systematic == 'ZJetsSF_up':
	ZJetsSF += ZJetsSFErr
if systematic == 'ZJetsSF_down':
	ZJetsSF -= ZJetsSFErr
if systematic == 'zeroB':
	ZJetsSF = 1.0


if systematic == 'otherMC_up':
	otherMCSF = 1.2
if systematic == 'otherMC_down':
	otherMCSF = 0.8



#import array
#binarray = array.array('d')
#binarray.fromlist([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,300])

# load cross-sections and N gen as global variables
execfile('SF.py')

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
#	canvas = ROOT.TCanvas('c1','c1',640,800) #old value, changed to match shape of fit plots

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
#	latex.SetTextAlign(12)
	latex.SetTextSize(0.037)
#	latex.SetLineWidth(2)
	
	for var in varlist:
		legend = ROOT.TLegend(0.71, 0.9 - 0.05*(len(MCTemplateList) + len(SignalTemplateZoomList)), 0.94, 0.9)
		legend.SetBorderSize(0)
		legend.SetFillColor(0)
		
		if dataTemplate is not None:
			legend.AddEntry(dataTemplate.histList[var], dataTemplate.legName, 'pl')
		
		# MC templates listed in the order they appear in legend
		for mc in MCTemplateList:
			mcHist = mc.histList[var]
			legend.AddEntry(mcHist, mc.legName, 'f')
		
		stack = ROOT.THStack('stack_'+var,var)
		# reverse order for stack to be consistent with legend
		MCTemplateList.reverse()
		for mc in MCTemplateList:
			mcHist = mc.histList[var]
			if var == 'M3' and 'barrel' in outDirName:
				mcHist.Rebin(4)
			stack.Add(mcHist)
		MCTemplateList.reverse()

		if dataTemplate is not None:
			if var == 'M3':
				if 'barrel' in outDirName:				
					dataTemplate.GetYaxis().SetTitle('Events / 40 GeV')
					dataTemplate.histList[var].Rebin(4)
				else:
					dataTemplate.GetYaxis().SetTitle('Events / 10 GeV')
			if dataTemplate.histList[var].GetMaximum() > stack.GetMaximum():
				stack.SetMaximum(dataTemplate.histList[var].GetMaximum())

		if 'cut_flow' in var: # or 'MET' in var:
			canvas.SetLogy(1)
			stack.SetMinimum(100)
		else:
			canvas.SetLogy(0)
		
		stack.Draw('HIST')
		if lep+'1RelIso' in var:
			canvas.SetLogy()
		if 'ele1MVA' in var:
			stack.GetXaxis().SetRangeUser(-1.1,0.0)

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
				legend.AddEntry(sigHist, signal.legName + ' x ' + str(zoom), 'f')
			else:
				legend.AddEntry(sigHist, signal.legName, 'f')
		if 'cut_flow' not in var:
			legend.Draw()
	
               	# ROOT.TGaxis.SetMaxDigits(3)	
		# labelcms = ROOT.TPaveText(0.14,1.-canvas.GetTopMargin(),0.6,1.05-canvas.GetTopMargin(),"NDCBR")
		# labelcms.SetTextAlign(12);
		# labelcms.SetTextSize(0.045);
		# labelcms.SetFillColor(ROOT.kWhite);
		# labelcms.SetFillStyle(0);
		# if isMuon: labelcms.AddText("CMS Preliminary, L=19.7 fb^{-1}, #sqrt{s} = 8 TeV, #mu+jets");
		# if isElectron: labelcms.AddText("CMS Preliminary, L=19.7 fb^{-1}, #sqrt{s} = 8 TeV, e+jets");
		# labelcms.SetBorderSize(0);
		# labelcms.Draw()
		# latex.DrawLatex(0.3,0.93,'#splitline{L=19.7 fb^{-1} #sqrt{s} = 8 TeV}{CMS Preliminary}')
		channelText = ""
		if isMuon: channelText = "#mu+jets"
		if isElectron: channelText = "e+jets"

		CMS_lumi.extraText = channelText
		CMS_lumi.writeExtraText = True

		CMS_lumi.CMS_lumi(canvas, 2, 11)



		if not isSyst:
			canvas.Update();
			canvas.RedrawAxis();
#			canvas.GetFrame().Draw();

			canvas.Print(outDirName+'/'+var+".pdf",".pdf");
			canvas.Print(outDirName+'/'+var+".png",".png");
#			canvas.SaveAs(outDirName+'/'+var+'.png')
		
def loadDataTemplate(varlist, inputDir, prefix):
	templPrefix = inputDir+prefix
	DataTempl = distribution('Data', 'Data', [
		(templPrefix+'Data_a.root', 1),
		(templPrefix+'Data_b.root', 1),
		(templPrefix+'Data_c.root', 1),
		(templPrefix+'Data_d.root', 1),
		], varlist)
	return DataTempl

def loadQCDTemplate(varlist, inputDir, prefix):
	templPrefix = inputDir+prefix
	QCD_sf = QCDSF
	QCDTempl = distribution('QCD', 'QCD', [
		(templPrefix+'Data_a.root', QCD_sf),
		(templPrefix+'Data_b.root', QCD_sf),
		(templPrefix+'Data_c.root', QCD_sf),
		(templPrefix+'Data_d.root', QCD_sf),
		(templPrefix+'TTJets1l.root', -1 * QCD_sf * TopSF * gSF * TTJets1l_xs/TTJets1l_num),
		(templPrefix+'TTJets2l.root', -1 * QCD_sf * TopSF * gSF * TTJets2l_xs/TTJets2l_num),
		######## Added in all other channels of MC, previously just ttjets 1l and 2l removed to get QCD template ######## 
		(templPrefix+'TTJetsHad.root', -1 * QCD_sf * TopSF * gSF * TTJetsHad_xs/TTJetsHad_num),
		(templPrefix+'TTGamma.root', -1 * QCD_sf * TopSF * gSF * newTTgamma_xs/newTTgamma_num),
		(templPrefix+'SingleT_t.root',     -1 * QCD_sf * otherMCSF * gSF * SingTopT_xs/SingTopT_num),
		(templPrefix+'SingleT_s.root',     -1 * QCD_sf * otherMCSF * gSF * SingTopS_xs/SingTopS_num),
		(templPrefix+'SingleT_tw.root',    -1 * QCD_sf * otherMCSF * gSF * SingToptW_xs/SingToptW_num),
		(templPrefix+'SingleTbar_t.root',  -1 * QCD_sf * otherMCSF * gSF * SingTopbarT_xs/SingTopbarT_num),
		(templPrefix+'SingleTbar_s.root',  -1 * QCD_sf * otherMCSF * gSF * SingTopbarS_xs/SingTopbarS_num),
		(templPrefix+'SingleTbar_tw.root', -1 * QCD_sf * otherMCSF * gSF * SingTopbartW_xs/SingTopbartW_num),
		(templPrefix+'W2Jets.root', -1 * QCD_sf * WJetsSF * gSF * W2Jets_xs/W2Jets_num),
		(templPrefix+'W3Jets.root', -1 * QCD_sf * WJetsSF * gSF * W3Jets_xs/W3Jets_num),
		(templPrefix+'W4Jets.root', -1 * QCD_sf * WJetsSF * gSF * W4Jets_xs/W4Jets_num),
		(templPrefix+'ZJets.root',  -1 * QCD_sf * ZJetsSF * otherMCSF * gSF * ZJets_xs/ZJets_num),
		(templPrefix+'Zgamma.root', -1 * QCD_sf * otherMCSF * gSF * Zgamma_xs/Zgamma_num),
		(templPrefix+'Wgamma.root', -1 * QCD_sf * otherMCSF * gSF * Wgamma_xs/Wgamma_num),

	], varlist, ROOT.kYellow)
	return QCDTempl
        
def loadMCTemplates(varList, inputDir, prefix, titleSuffix, fillStyle):
	templPrefix = inputDir+prefix
	
	MCtemplates = {}
	
	#MCtemplates['WHIZARD'] = distribution('TTGamma'+titleSuffix, [
	#	(templPrefix+'WHIZARD.root', TopSF*gSF*TTgamma_xs/WHIZARD_num)
	#	], varList, 98, fillStyle)
	
#	MCtemplates['WHIZARD'] = distribution('TTGamma'+titleSuffix, [ ##Changed the name
	MCtemplates['WHIZARD'] = distribution('TTGamma'+titleSuffix, 't#bar{t}+#gamma', [
		(templPrefix+'TTGamma.root', TopSF*gSF*newTTgamma_xs/newTTgamma_num)
		], varList, ROOT.kRed +1, fillStyle)
	
#	MCtemplates['TTJets'] = distribution('TTJets'+titleSuffix, [ ##Changed the name
	MCtemplates['TTJets'] = distribution('TTJets'+titleSuffix, 't#bar{t}+jets', [
		(templPrefix+'TTJets1l.root', TopSF*gSF*TTJets1l_xs/TTJets1l_num),
		(templPrefix+'TTJets2l.root', TopSF*gSF*TTJets2l_xs/TTJets2l_num),
		(templPrefix+'TTJetsHad.root', TopSF*gSF*TTJetsHad_xs/TTJetsHad_num),
		], varList ,ROOT.kRed -7, fillStyle)
	###################################
	#return MCtemplates
	###################################
	nonWJetsSF = 1.0
		
#	MCtemplates['Vgamma'] = distribution('Vgamma'+titleSuffix, [ ##changed
	MCtemplates['Vgamma'] = distribution('Vgamma'+titleSuffix, 'V+#gamma', [ 
        (templPrefix+'Zgamma.root', otherMCSF*gSF*Zgamma_xs/Zgamma_num),
        (templPrefix+'Wgamma.root', otherMCSF*WgammaSF*gSF*Wgamma_xs/Wgamma_num),
    #    (templPrefix+'WWgamma.root', gSF*WWgamma_xs/WWgamma_num),
        ], varList, ROOT.kGray, fillStyle)

#	MCtemplates['Zgamma'] = distribution('Zgamma'+titleSuffix, [
	MCtemplates['Zgamma'] = distribution('Zgamma'+titleSuffix, 'Z+#gamma', [
        (templPrefix+'Zgamma.root', otherMCSF*gSF*Zgamma_xs/Zgamma_num),
        ], varList, ROOT.kAzure+3, fillStyle)

#	MCtemplates['Wgamma'] = distribution('Wgamma'+titleSuffix, [
	MCtemplates['Wgamma'] = distribution('Wgamma'+titleSuffix, 'W+#gamma', [
        (templPrefix+'Wgamma.root', otherMCSF*WgammaSF*gSF*Wgamma_xs/Wgamma_num),
        ], varList, ROOT.kGray, fillStyle)

#	MCtemplates['SingleTop'] = distribution('SingleTop'+titleSuffix, [
	MCtemplates['SingleTop'] = distribution('SingleTop'+titleSuffix, 'Single Top', [
		(templPrefix+'SingleT_t.root',      otherMCSF*gSF*SingTopT_xs/SingTopT_num),
		(templPrefix+'SingleT_s.root',      otherMCSF*gSF*SingTopS_xs/SingTopS_num),
		(templPrefix+'SingleT_tw.root',     otherMCSF*gSF*SingToptW_xs/SingToptW_num),
		(templPrefix+'SingleTbar_t.root',   otherMCSF*gSF*SingTopbarT_xs/SingTopbarT_num),
		(templPrefix+'SingleTbar_s.root',   otherMCSF*gSF*SingTopbarS_xs/SingTopbarS_num),
		(templPrefix+'SingleTbar_tw.root',  otherMCSF*gSF*SingTopbartW_xs/SingTopbartW_num),
		], varList, ROOT.kMagenta, fillStyle)
	
#	MCtemplates['WJets'] = distribution('WJets'+titleSuffix, [
	MCtemplates['WJets'] = distribution('WJets'+titleSuffix, 'W+jets', [
        #(templPrefix+'WJets.root', WJetsSF*gSF*WJets_xs/WJets_num),
		(templPrefix+'W2Jets.root', WJetsSF*gSF*W2Jets_xs/W2Jets_num),
		(templPrefix+'W3Jets.root', WJetsSF*gSF*W3Jets_xs/W3Jets_num),
		(templPrefix+'W4Jets.root', WJetsSF*gSF*W4Jets_xs/W4Jets_num),
		], varList, ROOT.kGreen -3, fillStyle)

	######## Added back in the ZJetsSF scaling ######## 
#	MCtemplates['ZJets'] = distribution('ZJets'+titleSuffix, [
	MCtemplates['ZJets'] = distribution('ZJets'+titleSuffix, 'Z+jets', [
		(templPrefix+'ZJets.root',ZJetsSF*otherMCSF*gSF*ZJets_xs/ZJets_num)], varList, ROOT.kAzure-2, fillStyle)
	return MCtemplates

def saveAccTemplates(inputDir, outFileName):
	varList = ['MCcategory']

	AccTemplates = {}
	
	AccTemplates['TTGamma'] = distribution('TTGamma_signal', 'TTGamma', [
		(inputDir+'hist_1pho_rs_barrel_top_TTGamma.root', 1.0),
		], varList, 97)
		
	AccTemplates['TTGamma_presel'] = distribution('TTGamma_presel', 'TTGamma', [
		(inputDir+'hist_1pho_top_TTGamma.root', 1.0),
		], varList, 97)
	AccTemplates['TTJets1l'] = distribution('TTJets1l_presel', '', [
		(inputDir+'hist_1pho_top_TTJets1l.root', 1.0),
		], varList ,11)
	AccTemplates['TTJets2l'] = distribution('TTJets2l_presel', '', [
		(inputDir+'hist_1pho_top_TTJets2l.root', 1.0),
		], varList ,11)
	AccTemplates['TTJetsHad'] = distribution('TTJetsHad_presel', '', [
		(inputDir+'hist_1pho_top_TTJetsHad.root', 1.0),
		], varList ,11)
	
	saveTemplatesToFile(AccTemplates.values(), varList, outFileName)

def saveTTgammaAccTemplates(inputDir, outFileName):
	varList = ['MCcategory','genPhoRegionWeight_1l_2l']

	print inputDir

	AccTemplates = {}
	
	AccTemplates['TTGamma_presel'] = distribution('TTGamma_presel', '', [
		(inputDir+'hist_1pho_top_TTGamma.root', 1.0),
		], varList, 97)
	# AccTemplates['TTJets1l'] = distribution('TTJets1l_presel', '', [
	# 	(inputDir+'hist_1pho_top_TTJets1l.root', 1.0),
	# 	], varList ,11)
	# AccTemplates['TTJets2l'] = distribution('TTJets2l_presel', '', [
	# 	(inputDir+'hist_1pho_top_TTJets2l.root', 1.0),
	# 	], varList ,11)
	# AccTemplates['TTJetsHad'] = distribution('TTJetsHad_presel', '', [
	# 	(inputDir+'hist_1pho_top_TTJetsHad.root', 1.0),
	# 	], varList ,11)
	
	saveTemplatesToFile(AccTemplates.values(), varList, outFileName)

def saveNoMETTemplates(inputDir, inputData, outFileName, histName):
	varList = ['MET','MET_low','M3',lep+'1RelIso','ele1MVA']
	DataTempl = loadDataTemplate(varList, inputData, histName)
	MCTemplDict = loadMCTemplates(varList, inputDir, histName,'',1001)
	MCTempl = []
	MCTempl.append(MCTemplDict['WHIZARD'])
	MCTempl.append(MCTemplDict['TTJets'])
	MCTempl.append(MCTemplDict['Vgamma'])
	MCTempl.append(MCTemplDict['SingleTop'])
	MCTempl.append(MCTemplDict['WJets'])
	MCTempl.append(MCTemplDict['ZJets'])
	print inputDir
	print histName
	saveTemplatesToFile([DataTempl] + MCTempl, varList, outFileName)

def saveBarrelFitTemplates(inputDir, inputData,  outFileName):
	varList = ['MET','MET_low','M3','photon1ChHadSCRIso', 'photon1ChHadRandIso', 'photon1_Sigma_ChSCRIso']
	DataTempl_b = loadDataTemplate(varList, inputData, 'hist_1pho_barrel_top_') #change 
	
	MCTempl_b = loadMCTemplates(varList, inputDir, 'hist_1pho_barrel_top_','',1001)	#change
	MCTempl_rs_b = loadMCTemplates(varList, inputDir, 'hist_1pho_rs_barrel_top_', '_signal', 1001) #change
	MCTempl_fe_b = loadMCTemplates(varList, inputDir, 'hist_1pho_fe_barrel_top_', '_electron', 3005)#change
	MCTempl_fjrb_b = loadMCTemplates(varList, inputDir, 'hist_1pho_fjrb_barrel_top_', '_fake', 3005)#change
	
	saveTemplatesToFile([DataTempl_b] +  MCTempl_b.values() + MCTempl_rs_b.values() + MCTempl_fe_b.values() + MCTempl_fjrb_b.values(), varList, outFileName)

def savePreselTemplates(inputDir, qcdDir, inputData, outFileName):
	if WJetsSF != 1.0 or TopSF != 1.0:
		print 'We want to save templates for M3 fit, but the SFs are not 1.0'
		print 'exiting'
		return
	
	varList = ['MET','MET_low','M3',]
	DataTempl = loadDataTemplate(varList, inputData, 'hist_1pho_top_')#change
	if QCDSF > 0.0001:
		QCDTempl = loadQCDTemplate(varList, qcdDir, 'hist_1pho_top_') #change
	else:
		print 'The purpose of this function is to save templates for M3 fit, without QCD it is useless'
	
	MCTemplDict = loadMCTemplates(varList, inputDir, 'hist_1pho_top_','',1001) #change
	MCTempl = []
	MCTempl.append(MCTemplDict['WHIZARD'])
	MCTempl.append(MCTemplDict['TTJets'])
	MCTempl.append(MCTemplDict['Vgamma'])
	MCTempl.append(MCTemplDict['SingleTop'])
	MCTempl.append(MCTemplDict['WJets'])
	MCTempl.append(MCTemplDict['ZJets'])
	if QCDSF > 0.0001:
		MCTempl.append(QCDTempl)
	saveTemplatesToFile([DataTempl] + MCTempl, varList, outFileName)

def makeQCDPlots(varList,qcdDir,outDir):
	DataTempl = loadDataTemplate(varList,qcdDir,'hist_1phoNoMET_top_')
	MCTemplDict = loadMCTemplates(varList, qcdDir, 'hist_1phoNoMET_top_','',1001) #NoMET change
        MCTempl = []
        MCTempl.append(MCTemplDict['WHIZARD'])
        MCTempl.append(MCTemplDict['TTJets'])
        MCTempl.append(MCTemplDict['Vgamma'])
        MCTempl.append(MCTemplDict['SingleTop'])
        MCTempl.append(MCTemplDict['WJets'])
        MCTempl.append(MCTemplDict['ZJets'])
	if WJetsSF == 1.0 and TopSF == 1.0:
                pass
        else:
        # save final templates, exactly as they are on the plots
                saveTemplatesToFile([DataTempl] + MCTempl, ['MET','MET_low','M3','WtransMass',lep+'1RelIso','genPhoRegionWeight','MCcategory'], 'templates_presel_scaled_QCD.root')
        plotTemplates( DataTempl, MCTempl, [], varList, outDir+'/presel')
	return

def makeAllPlots(varList, inputDir, qcdDir, dataDir, outDirName):
	# load templates PreSel	
	DataTempl = loadDataTemplate(varList, dataDir, 'hist_1pho_top_') #NoMET change
	if QCDSF > 0.0001:
		QCDTempl = loadQCDTemplate(varList, qcdDir, 'hist_1pho_top_') #NoMET change
	MCTemplDict = loadMCTemplates(varList, inputDir, 'hist_1pho_top_','',1001) #NoMET change
	MCTempl = []
	MCTempl.append(MCTemplDict['WHIZARD'])
	MCTempl.append(MCTemplDict['TTJets'])
	MCTempl.append(MCTemplDict['Vgamma'])
	MCTempl.append(MCTemplDict['Wgamma'])
	MCTempl.append(MCTemplDict['Zgamma'])
	MCTempl.append(MCTemplDict['SingleTop'])
	MCTempl.append(MCTemplDict['WJets'])
	MCTempl.append(MCTemplDict['ZJets'])
	if QCDSF > 0.0001:
		MCTempl.append(QCDTempl)
	
        if WJetsSF == 1.0 and TopSF == 1.0:
		pass
	else:	
	# save final templates, exactly as they are on the plots
		saveTemplatesToFile([DataTempl] + MCTempl, ['MET','MET_low','M3','WtransMass','genPhoRegionWeight','MCcategory'], 'templates_presel_scaled.root')
 	print "SF used :", "Top=", TopSF ,"WJets=",WJetsSF, "QCD=",QCDSF, "OtherMC=",otherMCSF	
	MCTempl.remove(MCTemplDict['Vgamma'])	
	plotTemplates( DataTempl, MCTempl, [], varList, outDirName+'/presel')
	
	
	shortVarList = varList[:]
	shortVarList.remove('cut_flow')
	shortVarList.remove('genPhoRegionWeight')
	shortVarList.append('ele1pho1Mass')

	region = 'barrel'
	# load templates
	DataTempl_b = loadDataTemplate(shortVarList, dataDir, 'hist_1pho_'+region+'_top_') #change
	if QCDSF > 0.0001:
		QCDTempl_b = loadQCDTemplate(shortVarList, qcdDir, 'hist_1pho_'+region+'_top_') #change
	MCTemplDict_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_'+region+'_top_','',1001)#change
	MCTempl_b = []
	MCTempl_b.append(MCTemplDict_b['WHIZARD'])
	MCTempl_b.append(MCTemplDict_b['TTJets'])
	MCTempl_b.append(MCTemplDict_b['Vgamma'])
	MCTempl_b.append(MCTemplDict_b['Wgamma'])
	MCTempl_b.append(MCTemplDict_b['Zgamma'])
	MCTempl_b.append(MCTemplDict_b['SingleTop'])
	MCTempl_b.append(MCTemplDict_b['WJets'])
	MCTempl_b.append(MCTemplDict_b['ZJets'])
	if QCDSF > 0.0001:
		MCTempl_b.append(QCDTempl_b)
	
	MCTempl_rs_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_rs_barrel_top_', '_signal', 1001)
	MCTempl_fe_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_fe_barrel_top_', '_electron', 3005)
	MCTempl_fjrb_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_fjrb_barrel_top_', '_fake', 3005)
 	print "SF after photon selection :", "Top=", TopSF ,"WJets=",WJetsSF, "QCD=",QCDSF, "OtherMC=",otherMCSF	
	# save final templates, exactly as they are on the plots and by categories
	saveTemplatesToFile([DataTempl_b] + MCTempl_b + MCTempl_rs_b.values() + MCTempl_fe_b.values() + MCTempl_fjrb_b.values(), 
		['MET','MET_low','M3','WtransMass','MCcategory','nJets','ele1pho1Mass'], 
		'templates_barrel_scaled.root'
		)
	
	MCTempl_b.remove(MCTemplDict_b['Vgamma'])
	plotTemplates( DataTempl_b, MCTempl_b, [], shortVarList, outDirName+'/'+region+'_samples')
	
	############################
	return
	############################


def makePhotonSelectionPlots(varList, inputDir, qcdDir, dataDir, outDirName):
 	print "SF used after photon M3 fit :", "Top=", TopSF ,"WJets=",WJetsSF, "QCD=",QCDSF, "OtherMC=",otherMCSF	

	shortVarList = varList[:]
	shortVarList.remove('cut_flow')
	shortVarList.remove('genPhoRegionWeight')
	
	region = 'barrel'
	# load templates
	DataTempl_b = loadDataTemplate(shortVarList, dataDir, 'hist_1pho_'+region+'_top_') #change
	if QCDSF > 0.0001:
		QCDTempl_b = loadQCDTemplate(shortVarList, qcdDir, 'hist_1pho_'+region+'_top_') #change
	MCTemplDict_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_'+region+'_top_','',1001)#change
	MCTempl_b = []
	MCTempl_b.append(MCTemplDict_b['WHIZARD'])
	MCTempl_b.append(MCTemplDict_b['TTJets'])
	MCTempl_b.append(MCTemplDict_b['Vgamma'])
	MCTempl_b.append(MCTemplDict_b['SingleTop'])
	MCTempl_b.append(MCTemplDict_b['WJets'])
	MCTempl_b.append(MCTemplDict_b['ZJets'])
	if QCDSF > 0.0001:
		MCTempl_b.append(QCDTempl_b)
	
	MCTempl_rs_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_rs_barrel_top_', '_signal', 1001)
	MCTempl_fe_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_fe_barrel_top_', '_electron', 3005)
	MCTempl_fjrb_b = loadMCTemplates(shortVarList, inputDir, 'hist_1pho_fjrb_barrel_top_', '_fake', 3005)
 	print "SF after photon selection :", "Top=", TopSF ,"WJets=",WJetsSF, "QCD=",QCDSF, "OtherMC=",otherMCSF	
	# save final templates, exactly as they are on the plots and by categories
	saveTemplatesToFile([DataTempl_b] + MCTempl_b + MCTempl_rs_b.values() + MCTempl_fe_b.values() + MCTempl_fjrb_b.values(), 
		['MET','MET_low','M3','WtransMass','MCcategory','nJets'], 
		'templates_barrel_scaled_afterPhotonM3.root'
		)
	
	plotTemplates( DataTempl_b, MCTempl_b, [], shortVarList, outDirName+'/'+region+'_samples_scaled')
	
	############################
	return
	############################


varList_all = ['nVtx',
			'MET','MET_low','Ht','WtransMass','M3', 
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
	InputHist = '/eos/uscms/store/user/dnoonan/EleHists_looseVeto/hist_bins'+outSuffix+'/'
#	QCDHist =   '/eos/uscms/store/user/dnoonan/EleHists_looseVeto/QCD_bins/'
	QCDHist =   '/eos/uscms/store/user/dnoonan/EleHists/QCD_bins/'
	DataHist =  '/eos/uscms/store/user/dnoonan/EleHists_looseVeto/hist_bins/'
if isMuon:
	InputHist = '/eos/uscms/store/user/dnoonan/MuHists_looseVeto/hist_bins'+outSuffix+'/'
#	QCDHist =   '/eos/uscms/store/user/dnoonan/MuHists_looseVeto/QCD_bins/'
	QCDHist =   '/eos/uscms/store/user/dnoonan/MuHists/QCD_bins/'
	DataHist =  '/eos/uscms/store/user/dnoonan/MuHists_looseVeto/hist_bins/'

######## Added in a printout of histogram locations, for easier tracking later on ######## 

print 'Input Histogram location:', InputHist
print 'QCD Histogram location:', QCDHist
print 'Data Histogram location:', DataHist


# TTJets and TTGamma acceptance histograms
saveAccTemplates(InputHist, 'ttbar_acceptance.root')
saveTTgammaAccTemplates(InputHist, 'ttgamma_acceptance.root')

### templates for data driven fit or closure test. No rescaling necessary
saveBarrelFitTemplates(InputHist, DataHist, 'templates_barrel.root')
templateFits.InputFilename = 'templates_barrel.root'
templateFits.fitData = False ## to do closure test
templateFits.isElectron = isElectron
templateFits.isMuon = isMuon
templateFits.NpseudoExp = 0
templateFits.doTheFit()
templateFits.fitData = True ## to do closure test
templateFits.NpseudoExp = 1500
#templateFits.doTheFit()

phoPurity,phoPurityError,MCfrac = templateFits.doTheFit()



# if isElectron:
# 	phoPurity, phoPurityError = 0.564158170272, 0.0633698012736
# if isMuon:
#phoPurity, phoPurityError = 0.531773010465, 0.0572147752547


if skipMET:
	print '*'*80
	print 'Stopping code before the MET fit'
	exit()


# for MET fit. No rescaling
if WJetsSF == 1.0 and TopSF == 1.0:
	saveNoMETTemplates(InputHist, DataHist, 'templates_presel_nomet.root', 'hist_1phoNoMET_top_')
	saveNoMETTemplates(QCDHist, QCDHist, 'templates_presel_nomet_qcd.root', 'hist_1phoNoMET_top_')

qcd_fit.qcdMETfile = 'templates_presel_nomet_qcd.root'
qcd_fit.normMETfile = 'templates_presel_nomet.root'

qcd_fit.setQCDconstantM3 = True
qcd_fit.setOtherMCconstantM3 = True

qcd_fit.M3BinWidth=10.

qcd_fit.isElectron = isElectron
qcd_fit.isMuon = isMuon

QCDSF,QCDSFerror_met = qcd_fit.doQCDfit()

print "QCD SF from MET fit is :" , QCDSF
#QCD_low_SF,QCD_low_SFerror = qcd_fit.doQCD_lowfit()

if skipAfterMET:
	print '*'*80
	print 'Stopping code after the MET fit'
	exit()

# for systematics of QCD fit
if systematic == 'QCD_up':
	QCDSF *= 1.5
if systematic == 'QCD_down':
	QCDSF *= 0.5
# save templates for M3 fit
savePreselTemplates(InputHist, QCDHist, DataHist, 'templates_presel.root')

# do M3 fit, update SF for Top and WJets
qcd_fit.M3file = 'templates_presel.root'
TopSF_m3, TopSFerror_m3, WJetsSF_m3, WJetsSFerror_m3,otherMCSF_m3,otherMCSFerror_m3, QCDSF_m3, QCDSFerror_m3 = qcd_fit.doM3fit()

if skipAfterM3:
	print '*'*80
	print 'Stopping code after the M3 fit'
	exit()

TopSF *= TopSF_m3 
WJetsSF *= WJetsSF_m3 
QCDSF *= QCDSF_m3
otherMCSF *= otherMCSF_m3

TopSFErr = (TopSFErr**2 + TopSFerror_m3**2)**0.5 
WJetsSFErr = (WJetsSFErr**2 + WJetsSFerror_m3**2)**0.5 
QCDSFErr = (QCDSFErr**2 + QCDSFerror_m3**2)**0.5 
otherMCSFErr = (otherMCSFErr**2 + otherMCSFerror_m3**2)**0.5 

print "SF used :", "Top=", TopSF ,"WJets=",WJetsSF, "QCD=",QCDSF, "OtherMC=",otherMCSF	

#QCDSF_photon,QCDSFerror_photon = vgamma_fit.doQCDfit_photon()
#TopSF_photon, TopSFerror_photon, WJetsSF_photon, WJetsSFerror_photon = vgamma_fit.doM3fit_photon()
makeAllPlots(varList_all, InputHist, QCDHist, DataHist, 'plots')
makeQCDPlots(varList_all+['ele1MVA'], QCDHist, 'QCD_plots')


mcEventsTable.inputFileName = "templates_barrel_scaled.root"
print 
print 'Event count, photon selection, before M3 fit'
print
mcEventsTable.printMCTable(TopSFErr/TopSF, WJetsSFErr/WJetsSF, WgammaSFErr/WgammaSF)
print


if skipPhoton:
	print '*'*80
	print 'Stopping code before the photon fit'
	sys.exit(0)

######## Change the vgamma fit to return also the top fraction for use in the likelihood fit ######## 
vgamma_fit.setQCDconstantM3 = True
vgamma_fit.setOtherMCconstantM3 = True

vgamma_fit.M3BinWidth=40.

vgamma_fit.isElectron = isElectron
vgamma_fit.isMuon = isMuon

######TopSF_photon, TopSFerror_photon, WJetsSF_photon, WJetsSFerror_photon, otherMCSF_photon, otherMCSFerror_photon, m3_topFrac, m3_topFracErr = vgamma_fit.doM3fit_photon()

TopSF_presel = TopSF
TopSF_presel_error = TopSFErr
WJetsSF_presel = WJetsSF
WJetsSF_presel_error = WJetsSFErr

TopSF_phoM3, TopSFerror_phoM3, BkgSF_phoM3, BkgSFerror_phoM3, otherMCSF_phoM3, otherMCSFerror_phoM3 = 1.,1.,1.,1.,1.,1.

TopSF_phoM3, TopSFerror_phoM3, BkgSF_phoM3, BkgSFerror_phoM3, otherMCSF_phoM3, otherMCSFerror_phoM3, m3_topFrac, m3_topFracErr = vgamma_fit.doM3fit_photon_3Templates()

TopSF *= TopSF_phoM3
WgammaSF *= BkgSF_phoM3
WJetsSF *= otherMCSF_phoM3
otherMCSF *= otherMCSF_phoM3

TopSFErr = (TopSFErr**2 + TopSFerror_phoM3**2)**0.5 
WJetsSFErr = (WJetsSFErr**2 + otherMCSFerror_phoM3**2)**0.5 
otherMCSFErr = (otherMCSFErr**2 + otherMCSFerror_phoM3**2)**0.5 
WgammaSFErr = BkgSFerror_phoM3

##### TopSF_phoM3, TopSFerror_phoM3, BkgSF_phoM3, BkgSFerror_phoM3, m3_topFrac, m3_topFracErr = vgamma_fit.doM3fit_photon_2Templates()

makePhotonSelectionPlots(varList_all, InputHist, QCDHist, DataHist, 'plots')

if skipQCDphoton:
	print '*'*80
	print 'Stopping code before the MET fit after photon selection'
	exit()

saveNoMETTemplates(InputHist, DataHist, 'templates_barrel_nomet.root', 'hist_1phoNoMET_barrel_top_')
saveNoMETTemplates(QCDHist, QCDHist, 'templates_barrel_nomet_qcd.root', 'hist_1phoNoMET_barrel_top_')

vgamma_fit.qcdMETfile = 'templates_barrel_nomet_qcd.root'
vgamma_fit.normMETfile = 'templates_barrel_nomet.root'
	
QCDSF_photon,QCDSFerror_photon = vgamma_fit.doQCDfit_photon_NoMET()
#QCD_low_SF_photon,QCD_low_SFerror_photon = vgamma_fit.doQCDlowfit_photon()

mcEventsTable.inputFileName = "templates_barrel_scaled_afterPhotonM3.root"
print 
print 'Event count, photon selection, after M3 fit'
print
mcEventsTable.printMCTable(TopSFErr/TopSF, WJetsSFErr/WJetsSF, WgammaSFErr/WgammaSF)

print

if skipCalc:
	print '*'*80
	print 'Stopping code before the likelihood fit (calc_the_answer)'
	sys.exit(0)

calc_the_answer_test.TTJets1l_num = TTJets1l_num
calc_the_answer_test.TTJets2l_num = TTJets2l_num
calc_the_answer_test.TTJetsHad_num = TTJetsHad_num

calc_the_answer_test.photnPurity = phoPurity
calc_the_answer_test.photnPurityErr = phoPurityError
calc_the_answer_test.eleFakeSF = 1.458
calc_the_answer_test.eleFakeSFErr = 0.1985
if systematic == 'EleFakeSF_up':
	calc_the_answer_test.eleFakeSF = 1.458+0.1985
if systematic == 'EleFakeSF_down':
	calc_the_answer_test.eleFakeSF = 1.458-0.1985

calc_the_answer_test.M3TopSF = TopSF_presel
calc_the_answer_test.M3TopSFErr = TopSF_presel_error
calc_the_answer_test.M3WJetsSF = WJetsSF_presel
calc_the_answer_test.M3WJetsSFErr = WJetsSF_presel_error
calc_the_answer_test.M3_photon_topFrac = m3_topFrac
calc_the_answer_test.M3_photon_topFracErr = m3_topFracErr

calc_the_answer_test.barrelFileName_M3fitscaled = 'templates_barrel_scaled_afterPhotonM3.root'

xsRatio, xsRatioErr, bestttgSF, bestttgSFErr, bestvgammaSF, bestvgammaSFErr, bestjgSF, bestjgSFErr = calc_the_answer_test.doTheCalculation()

mcEventsTable.ttgammaSF     = bestttgSF
mcEventsTable.vgammaSF      = bestvgammaSF
mcEventsTable.jetToPhotonSF = bestjgSF

print
print 'BestttgSF    = ',bestttgSF    , '+-', bestttgSFErr    
print 'BestvgammaSF = ',bestvgammaSF , '+-', bestvgammaSFErr 
print 'BestjgSF     = ',bestjgSF     , '+-', bestjgSFErr     

mcEventsTable.ttgammaSFerr     = bestttgSFErr
mcEventsTable.vgammaSFerr      = bestvgammaSFErr
mcEventsTable.jetToPhotonSFerr = bestjgSFErr

print 
print 'Event count, photon selection, likelihood Scale Factors Applied'
print
mcEventsTable.printMCTable(TopSFErr/TopSF, WJetsSFErr/WJetsSF, WgammaSFErr/WgammaSF)
