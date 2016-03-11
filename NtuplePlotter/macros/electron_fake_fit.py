import ROOT
import sys
import os

from Style import *
 
import CMS_lumi

drawRatio = True
padRatio = 0.25
padOverlap = 0.05
padGap = 0.03


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


useDiboson = False
useZeroB = True

if 'Diboson' in sys.argv:
	useDiboson = True
if 'WithBJet' in sys.argv:
	useZeroB = False

ROOT.gROOT.SetBatch()
openfiles = {}

def get1DHist(filename, histname):
	if filename not in openfiles:
		openfiles[filename] = ROOT.TFile(filename,'READ')
	file = openfiles[filename]
		
	hist = file.Get(histname)
	hist.SetDirectory(0)
	hist.SetFillColor(0)
	#hist.Sumw2()
	return hist


def makeFit(varname, varmin, varmax, signalHist, backgroundHist, dataHist, plotName):
	# RooFit variables
	sihihVar = ROOT.RooRealVar(varname, varname, varmin, varmax)
	sihihArgList = ROOT.RooArgList()
	sihihArgList.add(sihihVar)
	sihihArgSet = ROOT.RooArgSet()
	sihihArgSet.add(sihihVar)

	# dataHist.GetXaxis().SetTitle("e#gamma Inv. Mass (GeV)")
	# signalHist.GetXaxis().SetTitle("e#gamma Inv. Mass (GeV)")
	# backgroundHist.GetXaxis().SetTitle("e#gamma Inv. Mass (GeV)")
	# create PDFs
	signalDataHist = ROOT.RooDataHist('signalDataHist','signal RooDataHist', sihihArgList, signalHist)
	signalPdf = ROOT.RooHistPdf('signalPdf',varname+' of signal', sihihArgSet, signalDataHist)

	backgroundDataHist = ROOT.RooDataHist('backgroundDataHist','background RooDataHist', sihihArgList, backgroundHist)
	backgroundPdf = ROOT.RooHistPdf('backgroundPdf',varname+' of background', sihihArgSet, backgroundDataHist)

	# data
	dataDataHist = ROOT.RooDataHist('data '+varname, varname+' in Data', sihihArgList, dataHist)

	# signal fraction parameter
	sfname = 'signal fraction'
	if 'MET' in varname:
		sfname = 'multijet fraction'
	if 'M3' in varname:
		sfname = 'ttbar fraction'
	if 'ele1pho1Mass' in varname:
		sfname = 'e to #gamma fraction'
		
	signalFractionVar = ROOT.RooRealVar(sfname,sfname, 0.5, 0.0, 1.0)
	sumPdf = ROOT.RooAddPdf('totalPdf','signal and background', signalPdf, backgroundPdf, signalFractionVar)
	
	# fit
	sumPdf.fitTo( dataDataHist, ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.PrintLevel(-1))
	
	if plotName!='':
		# plot results
		H = 600; 
		W = 800; 

		canvas = ROOT.TCanvas('c1','c1',W,H)
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

		plotter = ROOT.RooPlot('myplot','',sihihVar,varmin,varmax,20) # nBins is dummy
		dataDataHist.plotOn(plotter, ROOT.RooFit.Name('data'))
		sumPdf.plotOn(plotter, ROOT.RooFit.Name('sum'), ROOT.RooFit.LineColor(ROOT.kRed))
		sumPdf.plotOn(plotter, ROOT.RooFit.Components('signalPdf'), ROOT.RooFit.Name('signal'), 
			ROOT.RooFit.LineColor(ROOT.kGreen))
		sumPdf.plotOn(plotter, ROOT.RooFit.Components('backgroundPdf'), ROOT.RooFit.Name('background'), 
			ROOT.RooFit.LineColor(ROOT.kBlue))
		#		sumPdf.paramOn(plotter) # fix


		


#		leg = ROOT.TLegend(.7,.99-canvas.GetTopMargin()-0.07*4,.99-canvas.GetRightMargin(),.99-canvas.GetTopMargin())
		leg = ROOT.TLegend(0.71, 1. - 0.1 - 0.07*4, 0.94, 1-0.1)
		leg.SetFillColor(ROOT.kWhite)
		leg.SetLineColor(ROOT.kWhite)
		leg.AddEntry(plotter.findObject('data'), 'Data','p')
		leg.AddEntry(plotter.findObject('sum'), 'Sum','l')
		leg.AddEntry(plotter.findObject('signal'), 'Z#rightarrowee (e to #gamma)','l')
		leg.AddEntry(plotter.findObject('background'), 'Background','l')
		
		plotter.GetYaxis().SetTitle("Events / 10 GeV")
		plotter.GetXaxis().SetTitle("M(e,#gamma) (GeV)")
		plotter.Draw()
#		plotter.GetYaxis().SetTitleOffset(1.4)
		leg.Draw()

		channelText = "e+jets"

		CMS_lumi.channelText = channelText
		CMS_lumi.writeExtraText = True
		CMS_lumi.writeChannelText = True
		CMS_lumi.CMS_lumi(canvas, 2, 11)

		canvas.Update()
		canvas.RedrawAxis();
		canvas.Print(plotName, ".png")
		canvas.Print(plotName.replace('png', 'pdf'), ".pdf")

		print sumPdf
		total = dataHist.Integral(dataHist.FindBin(varmin), dataHist.FindBin(varmax))
		signalTot = total*signalFractionVar.getVal()
		print total, signalTot
		binsX = dataHist.FindBin(varmax)-dataHist.FindBin(varmin)
		a = backgroundHist.Clone("test")
		b = signalHist.Clone("test2")
		a.Scale((1-signalFractionVar.getVal())*total/backgroundHist.Integral())
		b.Scale((signalFractionVar.getVal())*total/signalHist.Integral())
		a.SetLineColor(ROOT.kRed)
		a.Add(b)
		ratio = dataHist.Clone("temp")
		ratio.Divide(a)

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
		pad1.SetLeftMargin( L/W )
		pad1.SetRightMargin( R/W )
		pad1.SetTopMargin( T/H/(1-padRatio) )
		pad1.SetBottomMargin( (padOverlap+padGap)/(1-padRatio+padOverlap) )
		pad2.SetLeftMargin( L/W )
		pad2.SetRightMargin( R/W )
		pad2.SetTopMargin( (padOverlap)/(padRatio+padOverlap) )
		pad2.SetBottomMargin( B/H/padRatio )

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

		canvasRatio.cd()
		pad1.cd()
		plotter.GetXaxis().SetLabelSize(0)
		plotter.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()/(1.-padRatio))
		plotter.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()/(1.-padRatio))
		plotter.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleYOffset()*(1.-padRatio))
		plotter.Draw()
		pad1.Update()
#		print pad1.GetY2()
#		plotter.SetMinimum(-0.02*pad1.GetY2())
#		plotter.GetYaxis().SetRangeUser(-0.02*pad1.GetY2(), pad1.GetY2())

		leg = ROOT.TLegend(0.71, 1. - 0.1/(1.-padRatio) - 0.07/(1.-padRatio)*4, 0.94, 1-0.1/(1.-padRatio))
		leg.SetFillColor(ROOT.kWhite)
		leg.SetLineColor(ROOT.kWhite)
		leg.AddEntry(plotter.findObject('data'), 'Data','p')
		leg.AddEntry(plotter.findObject('sum'), 'Sum','l')
		leg.AddEntry(plotter.findObject('signal'), 'Z#rightarrowee (e to #gamma)','l')
		leg.AddEntry(plotter.findObject('background'), 'Background','l')

		leg.Draw()
		pad2.cd()
		oneLine = ROOT.TF1("line","1",-100000,10000)
		oneLine.SetLineColor(ROOT.kBlack)
		oneLine.SetLineWidth(1)
		oneLine.SetLineStyle(2)
		ratio.Draw()
		oneLine.Draw("same")
		ratio.GetXaxis().SetRangeUser(varmin, varmax)
		ratio.GetYaxis().SetRangeUser(0.4,1.6)
		ratio.GetYaxis().SetNdivisions(503)

		ratio.GetYaxis().SetTitle('Data/MC')
		ratio.GetYaxis().CenterTitle()
		ratio.SetTitle('')
		ratio.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()/(padRatio+padOverlap))
		ratio.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()/(padRatio+padOverlap))
		ratio.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()/(padRatio+padOverlap))
		ratio.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()/(padRatio+padOverlap))
		ratio.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleYOffset()*(padRatio+padOverlap))
		CMS_lumi.CMS_lumi(canvasRatio, 2, 11)

		canvasRatio.Update()
		canvasRatio.RedrawAxis()
		canvasRatio.Print(plotName.replace('.png', '_ratio.png'), ".png")
		canvasRatio.Print(plotName.replace('.png', '_ratio.pdf'), ".pdf")


	print 'fit returned value ',signalFractionVar.getVal(),' +- ',signalFractionVar.getError()
	return (signalFractionVar.getVal(),signalFractionVar.getError())


def getInt_Err(hist):
	err = ROOT.Double(0.0)
	integr = hist.IntegralAndError(0,hist.GetNbinsX() + 1, err)
	return (integr,err)


def addEle_Pho_contributions(sample, all, ele, pho):
	if useZeroB: fileName = 'templates_barrel_scaled_zeroB.root'
	else: fileName = 'templates_barrel_scaled.root'
	
	allHist = get1DHist(fileName, sample+'_MET')
	eleHist = get1DHist(fileName, sample+'_electron_MET')
	phoHist = get1DHist(fileName, sample+'_signal_MET')
	
	allInt,allErr = getInt_Err(allHist)
	eleInt,eleErr = getInt_Err(eleHist)
	phoInt,phoErr = getInt_Err(phoHist)
	
	if allInt <  0.0001:
		return all, ele, pho
	
	newall = (all[0] + allInt, ( (all[1])**2 + allErr**2 )**0.5 )
	
	print 'sample, all, ele, pho : ',sample, allInt,'+-', allErr, '  ', eleInt,'+-', eleErr, '  ', phoInt,'+-', phoErr
	
	newele = (ele[0] + eleInt, ( (ele[1])**2 + (eleErr)**2 )**0.5 )
	
	newpho = (pho[0] + phoInt, ( (pho[1])**2 + (phoErr)**2 )**0.5 )
	
	return (newall, newele, newpho)


def getEleMass_template(sample):
	if useZeroB: fileName = 'templates_barrel_scaled_zeroB.root'
	else: fileName = 'templates_barrel_scaled.root'
	#allHist = get1DHist(fileName, sample+'_ele1pho1Mass')
	eleHist = get1DHist(fileName, sample+'_electron_ele1pho1Mass')
	
	#normHist = get1DHist('templates_barrel_scaled_zeroB.root', sample+'_ele1pho1Mass')
	
	#eleHist.Scale( normHist.Integral() / allHist.Integral() )
	return eleHist
	

ttbar_all = (0.0, 0.0)
ttbar_ele = (0.0, 0.0)
ttbar_pho = (0.0, 0.0)

bg_all = (0.0, 0.0)
bg_ele = (0.0, 0.0)
bg_pho = (0.0, 0.0)

ttbar_all, ttbar_ele, ttbar_pho = addEle_Pho_contributions('TTJets', ttbar_all, ttbar_ele, ttbar_pho)
ttbar_all, ttbar_ele, ttbar_pho = addEle_Pho_contributions('TTGamma', ttbar_all, ttbar_ele, ttbar_pho)

bg_all, bg_ele, bg_pho = addEle_Pho_contributions('WJets', bg_all, bg_ele, bg_pho)
bg_all, bg_ele, bg_pho = addEle_Pho_contributions('ZJets', bg_all, bg_ele, bg_pho)
if useDiboson: bg_all, bg_ele, bg_pho = addEle_Pho_contributions('Diboson', bg_all, bg_ele, bg_pho)
bg_all, bg_ele, bg_pho = addEle_Pho_contributions('Wgamma', bg_all, bg_ele, bg_pho)
bg_all, bg_ele, bg_pho = addEle_Pho_contributions('Zgamma', bg_all, bg_ele, bg_pho)
bg_all, bg_ele, bg_pho = addEle_Pho_contributions('SingleTop', bg_all, bg_ele, bg_pho)
# QCD has no MC information, expect no real photon and no electron
if useZeroB: qcdInt,qcdErr = getInt_Err(get1DHist('templates_barrel_scaled_zeroB.root', 'QCD_MET'))
else: qcdInt,qcdErr = getInt_Err(get1DHist('templates_barrel_scaled.root', 'QCD_MET'))
print 'zero b-tag selection'
print 'QCD total: ',qcdInt,qcdErr
bg_all = (bg_all[0] + qcdInt, ((bg_all[1])**2 + qcdErr**2)**0.5 )

print 'ttbar all, ele, pho ',ttbar_all, ttbar_ele, ttbar_pho
print 'bg all, ele, pho ',bg_all, bg_ele, bg_pho
print
print 'ttbar ele and pho fractions',ttbar_ele[0]/ttbar_all[0], '  ', ttbar_pho[0]/ttbar_all[0]
print 'bg ele and pho fractions',bg_ele[0]/bg_all[0], '  ', bg_pho[0]/bg_all[0]


# make templates for fitting:
# Vgamma electron fakes and ZJets electron fakes as "electrons"
# everything else as "other"
electronTempl = getEleMass_template('ZJets')
electronTempl.Add( getEleMass_template('Zgamma') )

if useZeroB: fileName = 'templates_barrel_scaled_zeroB.root'
else: fileName = 'templates_barrel_scaled.root'

print fileName
otherTempl = get1DHist(fileName, 'TTJets' + '_ele1pho1Mass')
otherTempl.Add( get1DHist(fileName, 'TTGamma' + '_ele1pho1Mass') )
otherTempl.Add( get1DHist(fileName, 'WJets' + '_ele1pho1Mass') )
otherTempl.Add( get1DHist(fileName, 'ZJets' + '_ele1pho1Mass') )
if useDiboson: otherTempl.Add( get1DHist(fileName, 'Diboson' + '_ele1pho1Mass') )
otherTempl.Add( get1DHist(fileName, 'Wgamma' + '_ele1pho1Mass') )
otherTempl.Add( get1DHist(fileName, 'Zgamma' + '_ele1pho1Mass') )
otherTempl.Add( get1DHist(fileName, 'SingleTop' + '_ele1pho1Mass') )
otherTempl.Add( get1DHist(fileName, 'QCD' + '_ele1pho1Mass') )

# subtract electron component from all histograms
otherTempl.Add( electronTempl, -1)

dataTempl = get1DHist(fileName, 'Data' + '_ele1pho1Mass')

lowmass = 20.0
highmass = 180.0
(eleFrac, eleFracErr) = makeFit('ele1pho1Mass', lowmass, highmass, electronTempl, otherTempl, dataTempl, 'egammaPlots/e_gamma_mass_fit.png')

print 'Finding fit range in bins'
lowbin = dataTempl.FindBin(lowmass + 0.01)
highbin = dataTempl.FindBin(highmass - 0.01)
print 'will intergate in the bin range: ',lowbin,highbin

dataErr = ROOT.Double(0.0)
dataInt = dataTempl.IntegralAndError(lowbin,highbin,dataErr)
fittedEle = dataInt * eleFrac

mcEleErr = ROOT.Double(0.0)
MCeleInt = electronTempl.IntegralAndError(lowbin,highbin,mcEleErr)
print 'Data integral: ', dataInt, '  MC total integral ', MCeleInt + otherTempl.Integral(lowbin,highbin)
print 'MC truth electrons: ',MCeleInt,'  Fitted elecrons: ', fittedEle
ratioRelErr = ( (dataErr/dataInt)**2 + (eleFracErr/eleFrac)**2 + (mcEleErr/MCeleInt)**2 )**0.5

print 'data over MC ratio: ',fittedEle/MCeleInt, '+/-',fittedEle/MCeleInt * ratioRelErr
ele_SF = fittedEle/MCeleInt


otherErr = ROOT.Double(0.0)
otherInt = otherTempl.IntegralAndError(lowbin,highbin,dataErr)
other_SF = (dataInt-fittedEle)/otherInt
print other_SF

ttjets = get1DHist(fileName, 'TTJets' + '_ele1pho1Mass')
ttgamma =  get1DHist(fileName, 'TTGamma' + '_ele1pho1Mass') 
wjets =  get1DHist(fileName, 'WJets' + '_ele1pho1Mass')
zjets =  get1DHist(fileName, 'ZJets' + '_ele1pho1Mass')
wgamma =  get1DHist(fileName, 'Wgamma' + '_ele1pho1Mass')
zgamma =  get1DHist(fileName, 'Zgamma' + '_ele1pho1Mass')
singletop =  get1DHist(fileName, 'SingleTop' + '_ele1pho1Mass')
qcd =  get1DHist(fileName, 'QCD' + '_ele1pho1Mass')

# e_ttjets = getEleMass_template('TTJets')
# e_ttgamma =  getEleMass_template('TTGamma') 
# e_wjets =  getEleMass_template('WJets')
e_zjets =  getEleMass_template('ZJets')
#e_wgamma =  getEleMass_template('Wgamma')
e_zgamma =  getEleMass_template('Zgamma')
#e_singletop =  getEleMass_template('SingleTop')

#ttjets.Add(    e_ttjets	  , -1.)
#ttgamma.Add(   e_ttgamma  , -1.)
#wjets.Add(     e_wjets	  , -1.)
zjets.Add(     e_zjets	  , -1.)
#wgamma.Add(    e_wgamma	  , -1.)
zgamma.Add(    e_zgamma	  , -1.)
#singletop.Add( e_singletop, -1.)

ttjets.Scale(other_SF)   
ttgamma.Scale(other_SF)  	  
wjets.Scale(other_SF)    
zjets.Scale(other_SF)    
wgamma.Scale(other_SF)   
zgamma.Scale(other_SF)   
singletop.Scale(other_SF)

#ttjets.Add(    e_ttjets	  , ele_SF)
#ttgamma.Add(   e_ttgamma  , ele_SF)	  
#wjets.Add(     e_wjets	  , ele_SF)
zjets.Add(     e_zjets	  , ele_SF)
#wgamma.Add(    e_wgamma	  , ele_SF)
zgamma.Add(    e_zgamma	  , ele_SF)
#singletop.Add( e_singletop, ele_SF)


# zjets.Add(electronZjets, ele_SF-1.)
# zgamma.Add(electronZgamma, ele_SF-1.)

ttgamma.SetFillColor(ROOT.kRed+1)
ttjets.SetFillColor(ROOT.kRed-7)
wjets.SetFillColor(ROOT.kGreen-3)
zjets.SetFillColor(ROOT.kAzure-2)
zgamma.SetFillColor(ROOT.kAzure+3)
wgamma.SetFillColor(ROOT.kGray)
singletop.SetFillColor(ROOT.kMagenta)
qcd.SetFillColor(ROOT.kYellow)

ttgamma.SetLineColor(ROOT.kRed+1)
ttjets.SetLineColor(ROOT.kRed-7)
wjets.SetLineColor(ROOT.kGreen-3)
zjets.SetLineColor(ROOT.kAzure-2)
zgamma.SetLineColor(ROOT.kAzure+3)
wgamma.SetLineColor(ROOT.kGray)
singletop.SetLineColor(ROOT.kMagenta)
qcd.SetLineColor(ROOT.kYellow)

stack = ROOT.THStack()
stack.Add(qcd)
stack.Add(wjets)
stack.Add(singletop)
stack.Add(wgamma)
stack.Add(ttjets)
stack.Add(ttgamma)
stack.Add(zgamma)
stack.Add(zjets)

legend = ROOT.TLegend(0.71, 0.9 - 0.05*(8), 0.94, 0.9)
legend.SetBorderSize(0)
legend.SetFillColor(0)

legend.AddEntry(dataTempl, 'Data', 'pl')
legend.AddEntry(zjets,     'Z+jets', 'f')
legend.AddEntry(zgamma,    'Z+#gamma', 'f')
legend.AddEntry(ttgamma,   't#bar{t}+#gamma', 'f')
legend.AddEntry(ttjets,    't#bar{t}+jets', 'f')
legend.AddEntry(wgamma,    'W+#gamma', 'f')
legend.AddEntry(singletop, 'Single Top', 'f')
legend.AddEntry(wjets,     'W+jets', 'f')
legend.AddEntry(qcd,       'QCD', 'f')


H = 600; 
W = 800; 

canvas = ROOT.TCanvas('c1','c1',W,H)
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

if dataTempl.GetMaximum() > stack.GetMaximum():
	stack.SetMaximum(dataTempl.GetMaximum())

stack.Draw('HIST')
dataTempl.Draw("same")

legend.Draw()

channelText = "e+jets"

CMS_lumi.extraText = channelText
CMS_lumi.writeExtraText = True
CMS_lumi.CMS_lumi(canvas, 2, 11)

canvas.Update()
canvas.RedrawAxis();
canvas.Print("egammaPlots/ele1pho1Mass_scaled.png", ".png")
canvas.Print("egammaPlots/ele1pho1Mass_scaled.pdf", ".pdf")


