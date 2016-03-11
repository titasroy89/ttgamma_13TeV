import ROOT
ROOT.gROOT.SetBatch()
import templateFits
import os
import sys
import CMS_lumi

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

isMuon = False
isElectron = False

channel = ''
if len(sys.argv)==2:
	if 'ele' in sys.argv[1]:
		channel = 'ele'
		isElectron = True
	if 'mu' in sys.argv[1]:
		channel = 'mu'
		isMuon = True
# ROOT.gStyle.SetOptTitle(0)
# ROOT.gStyle.SetPadLeftMargin(0.12)


templateFits.InputFilename = 'templates_barrel.root'
templateFits.fitData = False ## to do closure test
templateFits.NpseudoExp = 2000
#templateFits.NpseudoExp = 0

templateFits.savePlots = False

#templateFits.datasetsToMix = ['TTGamma','TTJets']
#templateFits.datasetsToMix = ['MGTTGamma','TTJets']

scales = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 4.5,5.0]
#scales = [0.1,1.,5.0]
graph = ROOT.TGraphErrors(len(scales))
gpoint = 0

for s in scales:
	print '*'*25
	print s
	templateFits.signalSF = s
	phoPurity,phoPurityError,MCtruth = templateFits.doTheFit()
	
	graph.SetPoint(gpoint, MCtruth, phoPurity)
	graph.SetPointError(gpoint, 0.0, phoPurityError)
	gpoint+=1
	print s,phoPurity,phoPurityError,MCtruth

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


	
graph.GetXaxis().SetTitle("MC Truth Photon Purity")
graph.GetYaxis().SetTitle("Fit Result Photon Purity")
graph.GetYaxis().SetTitleOffset(1)



graph.Draw('Ap')
graph.SetMarkerStyle(8)

_min = graph.GetXaxis().GetXmin()
_max = graph.GetXaxis().GetXmax()
line = ROOT.TLine(_min,_min,_max,_max)
line.Draw('SAME')

print canvas.GetTopMargin()
print canvas.GetLeftMargin()


channelText = ""
if isMuon: channelText = "#mu+jets"
if isElectron: channelText = "e+jets"

CMS_lumi.channelText = channelText
CMS_lumi.extraText = "Simulation"
CMS_lumi.writeExtraText = True
CMS_lumi.writeChannelText = True
		
CMS_lumi.CMS_lumi(canvas, 12, 11)
canvas.Update();
canvas.RedrawAxis();
# labelcms = ROOT.TPaveText(0.12,1.-canvas.GetTopMargin(),0.6,1.05-canvas.GetTopMargin(),"NDCBR")
# labelcms.SetTextAlign(12);
# labelcms.SetTextSize(0.045);
# labelcms.SetFillColor(ROOT.kWhite);
# labelcms.SetFillStyle(0);
# labelcms.AddText("CMS Simulation, #sqrt{s} = 8 TeV");
# labelcms.SetBorderSize(0);
# labelcms.Draw()

canvas.Print('fitplots/linearity_'+channel+'.pdf','.pdf')
canvas.Print('fitplots/linearity_'+channel+'.png','.png')

