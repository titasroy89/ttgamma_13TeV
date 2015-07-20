import ROOT
ROOT.gROOT.SetBatch()
import templateFits


ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPadLeftMargin(0.12)

templateFits.InputFilename = 'templates_barrel.root'
templateFits.fitData = False ## to do closure test
templateFits.NpseudoExp = 1000

#templateFits.datasetsToMix = ['TTGamma','TTJets']
#templateFits.datasetsToMix = ['MGTTGamma','TTJets']

scales = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 4.5,5.0]
#scales = [0.1,5.0]
graph = ROOT.TGraphErrors(len(scales))
gpoint = 0

for s in scales:
	templateFits.signalSF = s
	phoPurity,phoPurityError,MCtruth = templateFits.doTheFit()
	
	graph.SetPoint(gpoint, MCtruth, phoPurity)
	graph.SetPointError(gpoint, 0.0, phoPurityError)
	gpoint+=1
	
canvas = ROOT.TCanvas('c1','c1',640,640)

graph.GetXaxis().SetTitle("MC Truth Photon Purity")
graph.GetYaxis().SetTitle("Fit Result Photon Purity")
graph.GetYaxis().SetTitleOffset(1.4)



graph.Draw('Ap')
graph.SetMarkerStyle(8)

_min = graph.GetXaxis().GetXmin()
_max = graph.GetXaxis().GetXmax()
line = ROOT.TLine(_min,_min,_max,_max)
line.Draw('SAME')

print canvas.GetTopMargin()
print canvas.GetLeftMargin()

labelcms = ROOT.TPaveText(0.12,1.-canvas.GetTopMargin(),0.6,1.05-canvas.GetTopMargin(),"NDCBR")
labelcms.SetTextAlign(12);
labelcms.SetTextSize(0.045);
labelcms.SetFillColor(ROOT.kWhite);
labelcms.SetFillStyle(0);
labelcms.AddText("CMS Simulation, #sqrt{s} = 8 TeV");
labelcms.SetBorderSize(0);
labelcms.Draw()

canvas.SaveAs('linearityNewMG.png')
