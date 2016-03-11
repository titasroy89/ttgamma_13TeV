import os
import sys

SystNames = {
'Nsignal':'Likelihood fit uncertainty',
'Btag': 'b-tagging scale factor',
'EleFakeSF': 'Electron fake rate',
'JEC': 'Jet energy scale',
'JER': 'Jet energy resolution',
'PU': 'Pileup',
'QCD': 'Multijet estimate',
'ZJetsSF': '\\ZJets scale factor',
'EleEff': 'Electron efficiency',
'elesmear': 'Electron energy scale',
'otherMC': 'Background normalization',
'pho': 'Photon energy scale',
'toppt': 'Top \\pt reweighting',
'musmear': 'Muon energy scale',
'MuEff': 'Muon efficiency',
'PDF': 'PDF',
'Scale': 'Fact. and renorm. scale',
'Matching': 'ME/PS matching thresh.',
'TopMass': 'Top quark mass',
}

def getTable(args):

    directory = "/uscms/home/troy2012/TTGAMMA_trial/TTGammaSemiLep/NtuplePlotter/macros/ratioFiles/"

    takeMaxSyst = True
    twikiFormat = False #by default it prints the table in latex format

    if 'notMax' in args:
        args.remove('notMax')
        takeMaxSyst = False
    if 'twiki' in args:
        args.remove('twiki')
        twikiFormat = True

    if len(args) == 2:
        directory = args[1]


    temp = os.listdir(directory)

    fileList = list()
    
    systList = SystNames.keys()
    systList.remove('Nsignal')
    if 'ele' in directory:
        systList.remove('musmear')
        systList.remove('MuEff')
    if 'mu' in directory:
        systList.remove('elesmear')
        systList.remove('EleEff')

    systList.remove('PDF')
    # systList.remove('Scale')
    # systList.remove('Matching')
    # systList.remove('TopMass')

    for syst in systList:
        fileList.append("%s/ratio_%s_up.txt"%(directory,syst))
        fileList.append("%s/ratio_%s_down.txt"%(directory,syst))


    # for val in temp:
    #     if 'ratio_' in val and not 'nominal' in val:
    #         if 'ele' in directory and 'musmear' in val:
    #             continue
    #         elif 'ele' in directory and 'MuEff' in val:
    #             continue
    #         elif 'mu' in directory and 'elesmear' in val:
    #             continue
    #         elif 'mu' in directory and 'EleEff' in val:
    #             continue
    #         else:
    #             fileList.append(val)


    nominalVal = 0
    signalUnc = 0

    
    _file = open(directory+"ratio_nominal.txt","r")
    isResult = False
    isVisResult = False
    for line in _file:
        if 'Direct Fiducial Cross Section Value' in line:
            nominalVal = float(line.split()[-3])
            signalUnc =  float(line.split()[-1])

            
#    print nominalVal

    unc = {'Nsignal':[100*signalUnc/nominalVal,100*signalUnc/nominalVal],}
    values = {'nominal':nominalVal}
#place syst up in second spot and syst down in first spot
    upDown = {'up':1,'down':0}

    for systFile in fileList:

        systName = systFile.split('.')[0].split('_')[2:4]

        if not unc.has_key(systName[0]): 
            unc[systName[0]] = [0,0]
            values[systName[0]] = [0,0]

        _file = open(systFile,"r")
        for line in _file:
            if 'Direct Fiducial Cross Section Value' in line:
                systVal = float(line.split()[-3])

        unc[systName[0]][upDown[systName[1]]] = (systVal-nominalVal)/nominalVal*100
        values[systName[0]][upDown[systName[1]]] = systVal

    total = 0.0
    for i in unc:
        unc[i].append(abs(max(unc[i],key=abs)))
        total += max(unc[i],key=abs)**2

    table = ""
    if takeMaxSyst:
        if twikiFormat:
            table += '| * Source *   |       *Ratio Change (%)*       |\n'
        else:
            table += '\\begin{tabular}{l | c }\n'
            table += '\\hline\n'
            table += 'Source & Ratio Change (\\%) \\\\\n'
            table += '\\hline\n'
        for w in sorted(unc.items(),key=lambda x: x[1][2],reverse=True):
            if twikiFormat:
                table += "| %s  |  %.1f  |  %.1f  | \n" % (SystNames[w[0]], w[1][0][2], w[1][1][2])
            else:
                table += "%s   & %.1f \\\\ \n" % (SystNames[w[0]], w[1][2])
    
    
        if not twikiFormat:
            table += '\\hline \n'
            table += 'Total  &  %.1f  \\\\ \n' % (total**0.5)
            table += '\\hline \n'
            table += '\\end{tabular} \n'
        else:
            table += '| Total  |  %.1f  | \n' % (total**0.5)

    else:
        if twikiFormat:
            table += '|             |       *Ratio Change (%)*       ||     *Vis Ratio Change (%)*     || \n'
            table += '| * Source *  |  * Syst down *  |  * Syst up *  |  * Syst down *  |  * Syst up *  | \n'
        else:
            table += '\\begin{tabular}{l | c c | c c} \n'
            table += '\\hline \n'
            table += '& \\multicolumn{2}{c}{Ratio Change (\\%)} & \\multicolumn{2}{c}{Vis Ratio Change (\\%)} \\\\ \n'
            table += 'Source & Syst Down & Syst Up & Syst Down & Syst Up \\\\ \n'
            table += '\\hline \n'
        for w in sorted(unc.items(),key=lambda x: x[1][0][2],reverse=True):
            if twikiFormat:
                table +=  "| %s  |  %.1f  |  %.1f  |  %.1f  |  %.1f  | \n" % (w[0], w[1][0][0], w[1][0][1],  w[1][1][0], w[1][1][1])
            else:
                table +=  "%s   & %.1f & %.1f & %.1f & %.1f \\\\ \n" % (w[0], w[1][0][0], w[1][0][1],  w[1][1][0], w[1][1][1])
    
    
        if not twikiFormat:
            table +=  '\\hline \n'
            table +=  '\\end{tabular} \n'
    return table, total**0.5



if __name__=="__main__":
    
    print getTable(sys.argv)[0]
