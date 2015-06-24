from distribution_mod import distribution
import ROOT
import sys, os

import templateFits
import qcd_fit
import calc_the_answer
import vgamma_fit

ROOT.gROOT.SetBatch()
#########
#Style
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

isElectron = False
isMuon = False
lep = ''

 ######## Add an argument to determine if running on electrons or muons ######## 
isSyst = False
systematic = ''
if len(sys.argv) > 1:
	print sys.argv
	if sys.argv[1]=='e' or 'ele' in sys.argv[1]:
		isElectron = True
		lep = 'ele'
	elif sys.argv[1]=='mu' or 'muon' in sys.argv[1]:
		isMuon = True
		lep = 'mu'
	else:
		print '#'*30
		print 'First argument must specify either electron or muon'
		print '#'*30
		sys.exit(1)
	if len(sys.argv) > 2:
		systematic = sys.argv[2]
		if systematic != 'zeroB':
			isSyst = True
			sys.stdout = open('ratio_'+systematic+'.txt','w')
else:
	print '#'*30
	print 'At least one argument is required,'
	print 'Must specify if begin run on electrons (e or electron) or muons (mu or muon)'
	print '#'*30
	sys.exit(1)

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
WJetsSF = 1.0
TopSF = 1.0
QCDSF = 1.0
ZJetsSF = 1.0 #1.20 
ZJetsSFErr = 0.06
if systematic == 'ZJetsSF_up':
	ZJetsSF += ZJetsSFErr
if systematic == 'ZJetsSF_down':
	ZJetsSF -= ZJetsSFErr
if systematic == 'zeroB':
	ZJetsSF = 1.0

VgammaSF = 1.0 
otherMCSF = 1.0
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
	canvas = ROOT.TCanvas('c1','c1',640,800)
	
	latex = ROOT.TLatex()
	latex.SetNDC()
#	latex.SetTextAlign(12)
	latex.SetTextSize(0.037)
#	latex.SetLineWidth(2)
	
	for var in varlist:
		legend = ROOT.TLegend(0.71, 1 - 0.05*(1 + len(MCTemplateList) + len(SignalTemplateZoomList)), 0.95, 0.93)
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
		if lep+'1RelIso' in var:
			canvas.SetLogy()
		
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
	
               	ROOT.TGaxis.SetMaxDigits(3)	
		latex.DrawLatex(0.3,0.93,'#splitline{L=19.7 fb^{-1} #sqrt{s} = 8 TeV}{CMS Preliminary}')

		if not isSyst:
			canvas.SaveAs(outDirName+'/'+var+'.png')
		
def loadDataTemplate(varlist, inputDir, prefix):
	templPrefix = inputDir+prefix
	DataTempl = distribution('Data', [
		(templPrefix+'Data_a.root', 1),
		(templPrefix+'Data_b.root', 1),
		(templPrefix+'Data_c.root', 1),
		(templPrefix+'Data_d.root', 1),
		], varlist)
	return DataTempl

def loadQCDTemplate(varlist, inputDir, prefix):
	templPrefix = inputDir+prefix
	QCD_sf = QCDSF
	QCDTempl = distribution('QCD', [
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
	
	MCtemplates['WHIZARD'] = distribution('TTGamma'+titleSuffix, [
		(templPrefix+'TTGamma.root', TopSF*gSF*newTTgamma_xs/newTTgamma_num)
		], varList, ROOT.kRed +1, fillStyle)
	
	MCtemplates['TTJets'] = distribution('TTJets'+titleSuffix, [
		(templPrefix+'TTJets1l.root', TopSF*gSF*TTJets1l_xs/TTJets1l_num),
		(templPrefix+'TTJets2l.root', TopSF*gSF*TTJets2l_xs/TTJets2l_num),
		(templPrefix+'TTJetsHad.root', TopSF*gSF*TTJetsHad_xs/TTJetsHad_num),
		], varList ,ROOT.kRed -7, fillStyle)
	###################################
	#return MCtemplates
	###################################
	nonWJetsSF = 1.0
		
	MCtemplates['Vgamma'] = distribution('Vgamma'+titleSuffix, [
        (templPrefix+'Zgamma.root', otherMCSF*gSF*Zgamma_xs/Zgamma_num),
        (templPrefix+'Wgamma.root', otherMCSF*gSF*Wgamma_xs/Wgamma_num),
    #    (templPrefix+'WWgamma.root', gSF*WWgamma_xs/WWgamma_num),
        ], varList, ROOT.kGray, fillStyle)

	MCtemplates['SingleTop'] = distribution('SingleTop'+titleSuffix, [
		(templPrefix+'SingleT_t.root',      otherMCSF*gSF*SingTopT_xs/SingTopT_num),
        (templPrefix+'SingleT_s.root',      otherMCSF*gSF*SingTopS_xs/SingTopS_num),
        (templPrefix+'SingleT_tw.root',     otherMCSF*gSF*SingToptW_xs/SingToptW_num),
        (templPrefix+'SingleTbar_t.root',   otherMCSF*gSF*SingTopbarT_xs/SingTopbarT_num),
        (templPrefix+'SingleTbar_s.root',   otherMCSF*gSF*SingTopbarS_xs/SingTopbarS_num),
        (templPrefix+'SingleTbar_tw.root',  otherMCSF*gSF*SingTopbartW_xs/SingTopbartW_num),
		], varList, ROOT.kMagenta, fillStyle)
	
	MCtemplates['WJets'] = distribution('WJets'+titleSuffix, [
        #(templPrefix+'WJets.root', WJetsSF*gSF*WJets_xs/WJets_num),
		(templPrefix+'W3Jets.root', WJetsSF*gSF*W3Jets_xs/W3Jets_num),
		(templPrefix+'W4Jets.root', WJetsSF*gSF*W4Jets_xs/W4Jets_num),
		], varList, ROOT.kGreen -3, fillStyle)

	######## Added back in the ZJetsSF scaling ######## 
	MCtemplates['ZJets'] = distribution('ZJets'+titleSuffix, [
		(templPrefix+'ZJets.root',ZJetsSF*otherMCSF*gSF*ZJets_xs/ZJets_num)], varList, ROOT.kAzure-2, fillStyle)
	return MCtemplates

def saveAccTemplates(inputDir, outFileName):
	varList = ['MCcategory']
	AccTemplates = {}
	
	AccTemplates['TTGamma'] = distribution('TTGamma_signal', [
		(inputDir+'hist_1pho_rs_barrel_top_TTGamma.root', 1.0),
		], varList, 97)
		
	AccTemplates['TTGamma_presel'] = distribution('TTGamma_presel', [
		(inputDir+'hist_1pho_top_TTGamma.root', 1.0),
		], varList, 97)
	AccTemplates['TTJets1l'] = distribution('TTJets1l_presel', [
		(inputDir+'hist_1pho_top_TTJets1l.root', 1.0),
		], varList ,11)
	AccTemplates['TTJets2l'] = distribution('TTJets2l_presel', [
		(inputDir+'hist_1pho_top_TTJets2l.root', 1.0),
		], varList ,11)
	AccTemplates['TTJetsHad'] = distribution('TTJetsHad_presel', [
		(inputDir+'hist_1pho_top_TTJetsHad.root', 1.0),
		], varList ,11)
	
	saveTemplatesToFile(AccTemplates.values(), varList, outFileName)

def saveNoMETTemplates(inputDir, inputData, outFileName):
	varList = ['MET','MET_low','M3']
	DataTempl = loadDataTemplate(varList, inputData, 'hist_1phoNoMET_top_')
	MCTemplDict = loadMCTemplates(varList, inputDir, 'hist_1phoNoMET_top_','',1001)
	MCTempl = []
	MCTempl.append(MCTemplDict['WHIZARD'])
	MCTempl.append(MCTemplDict['TTJets'])
	MCTempl.append(MCTemplDict['Vgamma'])
	MCTempl.append(MCTemplDict['SingleTop'])
	MCTempl.append(MCTemplDict['WJets'])
	MCTempl.append(MCTemplDict['ZJets'])
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
        print "SF used : ", TopSF, WJetsSF, QCDSF, otherMCSF	
	plotTemplates( DataTempl, MCTempl, [], varList, outDirName+'/presel')
	
	
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
 	print "SF after photon selection :", TopSF ,WJetsSF, QCDSF	
	# save final templates, exactly as they are on the plots and by categories
	saveTemplatesToFile([DataTempl_b] + MCTempl_b + MCTempl_rs_b.values() + MCTempl_fe_b.values() + MCTempl_fjrb_b.values(), 
		['MET','MET_low','M3','WtransMass','MCcategory','nJets'], 
		'templates_barrel_scaled.root'
		)
	
	plotTemplates( DataTempl_b, MCTempl_b, [], shortVarList, outDirName+'/'+region+'_samples')
	
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
if systematic in ['Btag_down','Btag_up','EleEff_down','EleEff_up','JEC_down','JEC_up','JER_down','JER_up','PU_down','PU_up','elesmear_down','elesmear_up','pho_down','pho_up','toppt_down','toppt_up']:
	outSuffix = '_'+systematic
else:
	outSuffix = ''

InputHist = '../../hist_bin'+outSuffix+'/'
QCDHist = '../../QCD_bin/'
DataHist = '../../hist_bin/'

######## Added in a printout of histogram locations, for easier tracking later on ######## 

print 'Input Histogram location:', InputHist
print 'QCD Histogram location:', QCDHist
print 'Data Histogram location:', DataHist


# TTJets and TTGamma acceptance histograms
saveAccTemplates(InputHist, 'ttbar_acceptance.root')

### templates for data driven fit or closure test. No rescaling necessary
saveBarrelFitTemplates(InputHist, DataHist, 'templates_barrel.root')
templateFits.InputFilename = 'templates_barrel.root'
templateFits.fitData = False ## to do closure test
##templateFits.NpseudoExp = 3000

######## Why is this fit not run? ######## 
#phoPurity,phoPurityError = 0.657, 0.0564 #0.506, 0.078  #0.564, 0.063 #### 0.556427532887, 0.0616417156454 ## auto binsize: 0.561220079533, 0.0529980243576
phoPurity,phoPurityError,MCfrac = templateFits.doTheFit()
#exit()
# for MET fit. No rescaling
if WJetsSF == 1.0 and TopSF == 1.0:
	saveNoMETTemplates(InputHist, DataHist, 'templates_presel_nomet.root')
	saveNoMETTemplates(QCDHist, QCDHist, 'templates_presel_nomet_qcd.root')

qcd_fit.qcdMETfile = 'templates_presel_nomet_qcd.root'
qcd_fit.normMETfile = 'templates_presel_nomet.root'

QCDSF,QCDSFerror_met = qcd_fit.doQCDfit()

print "QCD SF from MET fit is :" , QCDSF
#QCD_low_SF,QCD_low_SFerror = qcd_fit.doQCD_lowfit()

# for systematics of QCD fit
if systematic == 'QCD_up':
	QCDSF *= 2
if systematic == 'QCD_down':
	QCDSF /= 2
# save templates for M3 fit
savePreselTemplates(InputHist, QCDHist, DataHist, 'templates_presel.root')

# do M3 fit, update SF for Top and WJets
qcd_fit.M3file = 'templates_presel.root'
TopSF, TopSFerror, WJetsSF, WJetsSFerror,otherMCSF,otherMCSFerror, QCDSF_m3, QCDSFerror_m3 = qcd_fit.doM3fit()
QCDSF = QCDSF *QCDSF_m3

print TopSF, WJetsSF,otherMCSF, QCDSF
#QCDSF_photon,QCDSFerror_photon = vgamma_fit.doQCDfit_photon()
#TopSF_photon, TopSFerror_photon, WJetsSF_photon, WJetsSFerror_photon = vgamma_fit.doM3fit_photon()
makeAllPlots(varList_all, InputHist, QCDHist, DataHist, 'plots')
makeQCDPlots(varList_all, QCDHist, 'QCD_plots')
######## Change the vgamma fit to return also the top fraction for use in the likelihood fit ######## 
TopSF_photon, TopSFerror_photon, WJetsSF_photon, WJetsSFerror_photon, otherMCSF_photon, otherMCSFerror_photon, m3_topFrac, m3_topFracErr = vgamma_fit.doM3fit_photon()
print '*'*80
QCDSF_photon,QCDSFerror_photon = vgamma_fit.doQCDfit_photon()
#QCD_low_SF_photon,QCD_low_SFerror_photon = vgamma_fit.doQCDlowfit_photon()
#exit()

calc_the_answer.TTJets1l_num = TTJets1l_num
calc_the_answer.TTJets2l_num = TTJets2l_num
calc_the_answer.TTJetsHad_num = TTJetsHad_num

calc_the_answer.photnPurity = phoPurity
calc_the_answer.photnPurityErr = phoPurityError
calc_the_answer.eleFakeSF = 1.5
calc_the_answer.eleFakeSFErr = 0.2
if systematic == 'EleFakeSF_up':
	calc_the_answer.eleFakeSF = 1.5 + 0.2
if systematic == 'EleFakeSF_down':
	calc_the_answer.eleFakeSF = 1.5 - 0.2

calc_the_answer.M3TopSF = TopSF
calc_the_answer.M3TopSFErr = TopSFerror
calc_the_answer.M3WJetsSF = WJetsSF
calc_the_answer.M3WJetsSFErr = WJetsSFerror
calc_the_answer.M3_photon_topFrac = m3_topFrac
calc_the_answer.M3_photon_topFracErr = m3_topFracErr

calc_the_answer.doTheCalculation()

