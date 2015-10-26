import ROOT
import sys
import os

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
	sumPdf.fitTo( dataDataHist, ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.PrintLevel(-1) )
	
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

		leg = ROOT.TLegend(.7,.99-canvas.GetTopMargin()-.2-0.05*4,.99-canvas.GetRightMargin(),.99-canvas.GetTopMargin()-.2)
		leg.SetFillColor(ROOT.kWhite)
		leg.SetLineColor(ROOT.kWhite)
		leg.AddEntry(plotter.findObject('data'), 'Data','p')
		leg.AddEntry(plotter.findObject('sum'), 'Sum','l')
		leg.AddEntry(plotter.findObject('signal'), 'Z#rightarrowee (e to #gamma)','l')
		leg.AddEntry(plotter.findObject('background'), 'Background','l')
		
		plotter.GetXaxis().SetTitle("M(e,#gamma) (GeV)")
		plotter.Draw()
		plotter.GetYaxis().SetTitleOffset(1.4)
		leg.Draw()

		channelText = "e+jets"

		CMS_lumi.extraText = channelText
		CMS_lumi.writeExtraText = True
		CMS_lumi.CMS_lumi(canvas, 2, 33)


		canvas.Update()
		canvas.RedrawAxis();
		canvas.Print(plotName, ".png")
		canvas.Print(plotName.replace('png', 'pdf'), ".pdf")

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


