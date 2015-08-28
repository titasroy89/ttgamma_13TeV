import ROOT
import sys

doZeroB = False
if 'zeroB' in sys.argv:
	doZeroB = True
doEle = False
doMu = False


lep = ''
if 'mu' in sys.argv:
	doMu = True
	lep = 'mu'
if 'ele' in sys.argv:
	doEle = True
	lep = 'ele'

if (doEle and doMu) or not (doEle or doMu):
	print 'Pick one channel, electron or muon'
	sys.exit()


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
		c1 = ROOT.TCanvas('c1', 'c1', 800, 600)
		plotter = ROOT.RooPlot('myplot','',sihihVar,varmin,varmax,20) # nBins is dummy
		dataDataHist.plotOn(plotter, ROOT.RooFit.Name('data'))
		sumPdf.plotOn(plotter, ROOT.RooFit.Name('sum'), ROOT.RooFit.LineColor(ROOT.kRed))
		sumPdf.plotOn(plotter, ROOT.RooFit.Components('signalPdf'), ROOT.RooFit.Name('signal'), 
			ROOT.RooFit.LineColor(ROOT.kGreen))
		sumPdf.plotOn(plotter, ROOT.RooFit.Components('backgroundPdf'), ROOT.RooFit.Name('background'), 
			ROOT.RooFit.LineColor(ROOT.kBlue))
#		sumPdf.paramOn(plotter) # fix

		leg = ROOT.TLegend(.7,.99-c1.GetTopMargin()-0.05*4,.99-c1.GetRightMargin(),.99-c1.GetTopMargin())
		leg.SetFillColor(ROOT.kWhite)
		leg.SetLineColor(ROOT.kWhite)
		leg.AddEntry(plotter.findObject('data'), 'Data','p')
		leg.AddEntry(plotter.findObject('sum'), 'Sum','l')
		leg.AddEntry(plotter.findObject('signal'), 'Z+Jets','l')
		leg.AddEntry(plotter.findObject('background'), 'Background','l')

		labelcms = ROOT.TPaveText(c1.GetLeftMargin(),1.-c1.GetTopMargin(),0.6,1.05-c1.GetTopMargin(),"NDCBR")
		labelcms.SetTextAlign(12);
		labelcms.SetTextSize(0.045);
		labelcms.SetFillColor(ROOT.kWhite);
		labelcms.SetFillStyle(0);
		labelcms.AddText("CMS Preliminary, L=19.7 fb^{-1}, #sqrt{s} = 8 TeV");
		labelcms.SetBorderSize(0);

		if doEle:
			plotter.GetXaxis().SetTitle("M(ee) (GeV)")
		if doMu:
			plotter.GetXaxis().SetTitle("M(#mu#mu) (GeV)")
		plotter.Draw()
		leg.Draw()
		labelcms.Draw()
		plotter.GetYaxis().SetTitleOffset(1.4)
		c1.SaveAs(plotName)
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


