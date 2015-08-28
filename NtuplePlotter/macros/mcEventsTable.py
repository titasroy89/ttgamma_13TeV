import ROOT
from ROOT import *

latexFormat = True

inputFileName = "templates_barrel_scaled_afterPhotonM3.root"

photonGen = ['fake', 'electron', 'signal']

ttgammaSF = 1.
ttgammaSFerr = 0.

vgammaSF = 1.
vgammaSFerr = 0.

jetToPhotonSF = 1.
jetToPhotonSFerr = 0.

egammaSF = 1.47
egammaSFerr = 0.19

def printMCTable(TopSFRelErr = 0., WJetsSFRelErr = 0., WgammaSFRelErr = 0.):

    samples = ['TTGamma','TTJets', 'Wgamma', 'WJets', 'Zgamma', 'ZJets', 'SingleTop']
    inputFile = TFile(inputFileName,"READ")
    values = []
    valuesErrs = []

    for s in samples:
        tempVal = []
        tempErr = []
        for gen in photonGen:
            valSF = 1.
            histName = s + "_" + gen + "_MET"
            tempHist = inputFile.Get(histName)
            err = ROOT.Double(0.0)

            val = tempHist.IntegralAndError(-1,-1,err)

            if val > 0:
                err = (err/val)**2

            ## Apply TTgamma, Vgamma, and jetToPhoton scale factors
            if 'TTGamma' in s:
                valSF *= ttgammaSF
                err += (ttgammaSFerr/ttgammaSF)**2
            if 'Wgamma' in s or 'Zgamma' in s:
                valSF *= vgammaSF
                err += (vgammaSFerr/vgammaSF)**2
            if 'electron' in gen:
                valSF *= egammaSF
                err += (egammaSFerr/egammaSF)**2
            if 'fake' in gen:
                valSF *= jetToPhotonSF
                err += (jetToPhotonSFerr/jetToPhotonSF)**2
            if 'TTGamma' in s or 'TTJets' in s:
                err += TopSFRelErr**2
            if 'WJets' in s:
                err += WJetsSFRelErr**2
            if 'WGamma' in s:
                err += WgammaSFRelErr**2
                

            tempVal.append(val*valSF)
            tempErr.append(err**0.5*val*valSF)
        values.append(tempVal)
        valuesErrs.append(tempErr)

    tempHist = inputFile.Get("QCD_MET")
    err = ROOT.Double(0.0)
    val = tempHist.IntegralAndError(-1,-1,err)*jetToPhotonSF
    err = err/val
    values.append([val,0.,0.])
    valuesErrs.append([val*(err**2+(jetToPhotonSFerr/jetToPhotonSF)**2)**0.5,0.,0.])


    tempHist = inputFile.Get("Data_MET")
    dataErr = ROOT.Double(0.0)
    dataTot = tempHist.IntegralAndError(-1,-1,dataErr)

    inputFile.Close("R")

    totals = [0.,0.,0.,0.]
    totalsErr = [0.,0.,0.,0.]

    samples.append('QCD')
    if latexFormat:
        print '\\begin{tabular}{l c c c c}'
        print '\\hline'
        print 'Sample & Total & Fake jet & Electron & Photon \\\\'
        print '\\hline'
    else:
        print '| *Sample* | *Total* | *Fake jet* | *Electron* | *Photon* |'
    for i in range(len(values)):
        if latexFormat:
            print "%s & $%.1f \\pm %.1f$ & $%.1f \\pm %.1f$ & $%.1f \\pm %.1f$ & $%.1f \\pm %.1f$ \\\\" % (samples[i], (values[i][0]+values[i][1]+values[i][2]), (valuesErrs[i][0]**2+valuesErrs[i][1]**2+valuesErrs[i][2]**2)**0.5, values[i][0], valuesErrs[i][0], values[i][1], valuesErrs[i][1], values[i][2], valuesErrs[i][2])
        else:
            print "| %s | %.1f +- %.1f | %.1f +- %.1f | %.1f +- %.1f | %.1f +- %.1f |" % (samples[i], (values[i][0]+values[i][1]+values[i][2]), (valuesErrs[i][0]**2+valuesErrs[i][1]**2+valuesErrs[i][2]**2)**0.5, values[i][0], valuesErrs[i][0], values[i][1], valuesErrs[i][1], values[i][2], valuesErrs[i][2])

        totals[0] += values[i][0]
        totals[0] += values[i][1]
        totals[0] += values[i][2]
        totals[1] += values[i][0]
        totals[2] += values[i][1]
        totals[3] += values[i][2]

        totalsErr[0] += valuesErrs[i][0]**2
        totalsErr[0] += valuesErrs[i][1]**2
        totalsErr[0] += valuesErrs[i][2]**2
        totalsErr[1] += valuesErrs[i][0]**2
        totalsErr[2] += valuesErrs[i][1]**2
        totalsErr[3] += valuesErrs[i][2]**2

    totalsErr[0] = totalsErr[0]**0.5
    totalsErr[1] = totalsErr[1]**0.5
    totalsErr[2] = totalsErr[2]**0.5
    totalsErr[3] = totalsErr[3]**0.5
    if latexFormat:
        print '\\hline'
        print "Totals & $%.1f \\pm %.1f$ & $%.1f \\pm %.1f$ & $%.1f \\pm %.1f$ & $%.1f \\pm %.1f$ \\\\" % (totals[0], totalsErr[0], totals[1], totalsErr[1], totals[2], totalsErr[2], totals[3], totalsErr[3])
        print 'Data & $%.1f \\pm %.1f$ &  -  &  -  \\\\' % (dataTot, dataErr)
        print '\\hline'
        print '\\end{tabular}'
    else:
        print "| Totals | %.1f +- %.1f | %.1f +- %.1f | %.1f +- %.1f | %.1f +- %.1f |" % (totals[0], totalsErr[0], totals[1], totalsErr[1], totals[2], totalsErr[2], totals[3], totalsErr[3])
        print '| Data | %.1f +- %.1f |  -  |  -  | - |' % (dataTot, dataErr)

    return

#  LocalWords:  latexFormat
