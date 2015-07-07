import ROOT
from ROOT import *

inputFile = TFile("templates_barrel_scaled_afterPhotonM3.root",'READ')

samples = ['TTGamma','TTJets', 'Vgamma', 'WJets', 'ZJets', 'SingleTop']
samples = ['TTGamma','TTJets', 'Wgamma', 'WJets', 'Zgamma', 'ZJets', 'SingleTop']

photonGen = ['fake', 'electron', 'signal']

values = []
valuesErrs = []

for s in samples:
    tempVal = []
    tempErr = []
    for gen in photonGen:
        histName = s + "_" + gen + "_MET"
        print histName
        tempHist = inputFile.Get(histName)
        err = ROOT.Double(0.0)
        tempVal.append(tempHist.IntegralAndError(-1,-1,err))
        tempErr.append(err)
    values.append(tempVal)
    valuesErrs.append(tempErr)

tempHist = inputFile.Get("QCD_MET")
values.append([tempHist.IntegralAndError(-1,-1,err),0.,0.])
valuesErrs.append([err,0.,0.])

samples.append('QCD')
print '| samples | total | fake jet | electron | photon |'
for i in range(len(values)):
###    print "| %s | %.1f |" % (samples[i], values[i][0])
    print "| %s | %.2f | %.2f +- %.2f | %.1f +- %.2f | %.1f +- %.2f |" % (samples[i], (values[i][0]+values[i][1]+values[i][2]), values[i][0], valuesErrs[i][0], values[i][1], valuesErrs[i][1], values[i][2], valuesErrs[i][2])

