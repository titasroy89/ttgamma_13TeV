import ROOT
import sys
import os

doZeroB = False
if 'zeroB' in sys.argv:
	doZeroB = True
doEle = False
doMu = False

import CMS_lumi

drawRatio = True
padRatio = 0.25
padOverlap = 0.05
padGap = 0.01


lep = ''
if 'mu' in sys.argv:
	doMu = True
	lep = 'mu'
	channelText = "#mu#mu"

if 'ele' in sys.argv:
	doEle = True
	lep = 'ele'
	channelText = "ee"

if (doEle and doMu) or not (doEle or doMu):
	print 'Pick one channel, electron or muon'
	sys.exit()


ROOT.gROOT.SetBatch()
openfiles = {}

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
	if lep+'1pho1Mass' in varname:
		sfname = 'e to #gamma fraction'
	if lep+'1'+lep+'2Mass' in varname:
		sfname = 'Z to '+lep+lep+' fraction'
		
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
		leg.AddEntry(plotter.findObject('signal'), 'Z+jets','l')
		leg.AddEntry(plotter.findObject('background'), 'Background','l')


		if doEle:
			plotter.GetXaxis().SetTitle("M(ee) (GeV)")
		if doMu:
			plotter.GetXaxis().SetTitle("M(#mu#mu) (GeV)")
		plotter.Draw()
		leg.Draw()
#		plotter.GetYaxis().SetTitleOffset(1.4)
		channelText = ""
		if doMu: channelText = "#mu#mu"
		if doEle: channelText = "ee"

		CMS_lumi.channelText = channelText
		CMS_lumi.writeExtraText = True
		CMS_lumi.writeChannelText = True
		CMS_lumi.channelTextLocation = 3
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
	

# make templates for fitting:
# Vgamma and ZJets as "Z to ee"
# everything else as "other"
#fileDir = 'di_ele_cross_check_zeroB/'

if doEle:
	fileDir = 'di_ele_cross_check/'
	if doZeroB:
		fileDir = 'di_ele_cross_check_zeroB/'
if doMu:
	fileDir = 'di_mu_cross_check/'
	
fileName = fileDir+'plots/templates_presel.root'
hsuffix = '_'+lep+'1'+lep+'2Mass'

electronTempl = get1DHist(fileName, 'ZJets'+hsuffix)
electronTempl.Add( get1DHist(fileName, 'Zgamma'+hsuffix) )

otherTempl = get1DHist(fileName, 'TTJets'+hsuffix)
otherTempl.Add( get1DHist(fileName, 'TTGamma'+hsuffix) )
otherTempl.Add( get1DHist(fileName, 'WJets'+hsuffix) )
otherTempl.Add( get1DHist(fileName, 'Wgamma'+hsuffix) )
otherTempl.Add( get1DHist(fileName, 'SingleTop'+hsuffix) )


dataTempl = get1DHist(fileName, 'Data'+hsuffix)

lowmass = 20.0
highmass = 180.0
(eleFrac, eleFracErr) = makeFit(lep+'1'+lep+'2Mass', lowmass, highmass, electronTempl, otherTempl, dataTempl, fileDir+lep+'_'+lep+'_fit.png')

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



zjetSF = fittedEle/MCeleInt







ttjets = get1DHist(fileName, 'TTJets' + '_'+lep+'1'+lep+'2Mass')
ttgamma =  get1DHist(fileName, 'TTGamma' + '_'+lep+'1'+lep+'2Mass') 
wjets =  get1DHist(fileName, 'WJets' + '_'+lep+'1'+lep+'2Mass')
zjets =  get1DHist(fileName, 'ZJets' + '_'+lep+'1'+lep+'2Mass')
wgamma =  get1DHist(fileName, 'Wgamma' + '_'+lep+'1'+lep+'2Mass')
zgamma =  get1DHist(fileName, 'Zgamma' + '_'+lep+'1'+lep+'2Mass')
singletop =  get1DHist(fileName, 'SingleTop' + '_'+lep+'1'+lep+'2Mass')
#qcd =  get1DHist(fileName, 'QCD' + '_'+lep+'1'+lep+'2Mass')

zjets.Scale(zjetSF)

ttgamma.SetFillColor(ROOT.kRed+1)
ttjets.SetFillColor(ROOT.kRed-7)
wjets.SetFillColor(ROOT.kGreen-3)
zjets.SetFillColor(ROOT.kAzure-2)
zgamma.SetFillColor(ROOT.kAzure+3)
wgamma.SetFillColor(ROOT.kGray)
singletop.SetFillColor(ROOT.kMagenta)
#qcd.SetFillColor(ROOT.kYellow)

ttgamma.SetLineColor(ROOT.kRed+1)
ttjets.SetLineColor(ROOT.kRed-7)
wjets.SetLineColor(ROOT.kGreen-3)
zjets.SetLineColor(ROOT.kAzure-2)
zgamma.SetLineColor(ROOT.kAzure+3)
wgamma.SetLineColor(ROOT.kGray)
singletop.SetLineColor(ROOT.kMagenta)
#qcd.SetLineColor(ROOT.kYellow)

stack = ROOT.THStack()
#stack.Add(qcd)
stack.Add(wjets)
stack.Add(singletop)
stack.Add(zgamma)
stack.Add(wgamma)
stack.Add(ttjets)
stack.Add(ttgamma)
stack.Add(zjets)

legend = ROOT.TLegend(0.71, 0.9 - 0.05*(7), 0.94, 0.9)
legend.SetBorderSize(0)
legend.SetFillColor(0)

legendR = ROOT.TLegend(0.71, 1. - 0.1/(1.-padRatio) - 0.05/(1.-padRatio)*(7), 0.94, 1-0.1/(1.-padRatio))
legendR.SetBorderSize(0)
legendR.SetFillColor(0)


legend.AddEntry(dataTempl, 'Data', 'pl')
legend.AddEntry(zjets,     'Z+jets', 'f')
legend.AddEntry(ttgamma,   't#bar{t}+#gamma', 'f')
legend.AddEntry(ttjets,    't#bar{t}+jets', 'f')
legend.AddEntry(wgamma,    'W+#gamma', 'f')
legend.AddEntry(zgamma,    'Z+#gamma', 'f')
legend.AddEntry(singletop, 'Single Top', 'f')
legend.AddEntry(wjets,     'W+jets', 'f')
#legend.AddEntry(qcd,       'QCD', 'f')

legendR.AddEntry(dataTempl, 'Data', 'pl')
legendR.AddEntry(zjets,     'Z+jets', 'f')
legendR.AddEntry(ttgamma,   't#bar{t}+#gamma', 'f')
legendR.AddEntry(ttjets,    't#bar{t}+jets', 'f')
legendR.AddEntry(wgamma,    'W+#gamma', 'f')
legendR.AddEntry(zgamma,    'Z+#gamma', 'f')
legendR.AddEntry(singletop, 'Single Top', 'f')
legendR.AddEntry(wjets,     'W+jets', 'f')
#legend.AddEntry(qcd,       'QCD', 'f')


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
stack.GetXaxis().SetRangeUser(20,180)
if doEle: stack.GetXaxis().SetTitle("M(ee) (GeV)")
if doMu: stack.GetXaxis().SetTitle("M(#mu#mu) (GeV)")
stack.GetYaxis().SetTitle("Events / 5 GeV")


legend.Draw()

#CMS_lumi.extraText = channelText
CMS_lumi.writeExtraText = True
CMS_lumi.CMS_lumi(canvas, 2, 11)

canvas.Update()
canvas.RedrawAxis();
canvas.Print(fileDir+lep+"1"+lep+"2Mass_scaled.png", ".png")
canvas.Print(fileDir+lep+"1"+lep+"2Mass_scaled.pdf", ".pdf")


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

ratio = dataTempl.Clone("temp")
ratio.Divide(stack.GetStack().Last())

pad1.cd()
stack.Draw('HIST')

pad1.Update()
y2 = pad1.GetY2()
stack.SetMinimum(-0.02*y2)
pad1.Update()
dataTempl.Draw("esame")

legendR.Draw()

stack.SetTitle('')
stack.GetXaxis().SetLabelSize(0)
stack.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()/(1.-padRatio))
stack.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()/(1.-padRatio))
stack.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleYOffset()*(1.-padRatio))

pad1.Update()

pad2.cd()
ratio.Draw()

if doEle: ratio.GetXaxis().SetTitle("M(ee) (GeV)")
if doMu: ratio.GetXaxis().SetTitle("M(#mu#mu) (GeV)")

ratio.GetYaxis().SetTitle("Data/MC")
ratio.GetYaxis().CenterTitle()

ratio.SetTitle('')
ratio.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()/(padRatio+padOverlap))
ratio.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()/(padRatio+padOverlap))
ratio.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()/(padRatio+padOverlap))
ratio.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()/(padRatio+padOverlap))
ratio.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleYOffset()*(padRatio+padOverlap))

ratio.GetXaxis().SetRangeUser(20,180)
ratio.GetYaxis().SetRangeUser(0.38,1.62)
ratio.GetYaxis().SetNdivisions(503)

ratio.SetMarkerStyle(2)		
ratio.SetLineColor(ROOT.kBlack)
ratio.SetLineWidth(1)
oneLine = ROOT.TF1("oneLine","1",-1e6,1e6)
oneLine.SetLineColor(ROOT.kBlack)
oneLine.SetLineWidth(1)
oneLine.SetLineStyle(2)

ratio.Draw()        		
oneLine.Draw("same")

pad2.Update()
CMS_lumi.CMS_lumi(canvasRatio, 2, 11)
canvasRatio.Update()
canvasRatio.RedrawAxis();
canvasRatio.Print(fileDir+lep+"1"+lep+"2Mass_scaled_ratio.png", ".png")
canvasRatio.Print(fileDir+lep+"1"+lep+"2Mass_scaled_ratio.pdf", ".pdf")
