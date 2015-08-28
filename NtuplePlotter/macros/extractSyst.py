import os
import sys


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


    for val in temp:
        if 'ratio_' in val and not 'nominal' in val:
            if 'ele' in directory and 'musmear' in val:
                continue
            elif 'ele' in directory and 'MuEff' in val:
                continue
            elif 'mu' in directory and 'elesmear' in val:
                continue
            elif 'mu' in directory and 'EleEff' in val:
                continue
            else:
                fileList.append(val)


    nominalVal = 0
    signalUnc = 0

    nominalVisVal = 0
    signalVisUnc = 0
    
    _file = file(directory+"ratio_nominal.txt","r")
    isResult = False
    isVisResult = False
    for line in _file:
        if isResult:
            nominalVal = float(line.split()[0])
            signalUnc =  float(line.split()[-1])
            isResult = False
        if isVisResult:
            nominalVisVal = float(line.split()[0])
            signalVisUnc =  float(line.split()[-1])
            isVisResult = False
        if 'final answer: cross section ratio' in line:
            isResult = True
        if 'visible cross section ratio' in line:
            isVisResult = True

            
#    print nominalVal

    unc = {'Nsignal':[[100*signalUnc/nominalVal,100*signalUnc/nominalVal],[100*signalVisUnc/nominalVisVal,100*signalVisUnc/nominalVisVal]]}
    values = {'nominal':nominalVal}
    visvalues = {'nominal':nominalVisVal}
#place syst up in second spot and syst down in first spot
    upDown = {'up':1,'down':0}

    for systFile in fileList:

        systName = systFile.split('.')[0].split('_')[1:3]

        if not unc.has_key(systName[0]): 
            unc[systName[0]] = [[0,0],[0,0]]
            values[systName[0]] = [0,0]
            visvalues[systName[0]] = [0,0]

        _file = file(directory+systFile,"r")
        isResult = False
        for line in _file:
            if isResult:
                systVal = float(line.split()[0])
                isResult = False
            if isVisResult:
                systVisVal = float(line.split()[0])
                isVisResult = False
            if 'final answer' in line:
                isResult = True
            if 'visible cross section ratio' in line:
                isVisResult = True

        unc[systName[0]][0][upDown[systName[1]]] = (systVal-nominalVal)/nominalVal*100
        unc[systName[0]][1][upDown[systName[1]]] = (systVisVal-nominalVisVal)/nominalVisVal*100
        values[systName[0]][upDown[systName[1]]] = systVal
        visvalues[systName[0]][upDown[systName[1]]] = systVisVal

    total = 0.0
    totalVis = 0.0
    for i in unc:
        #    unc[i].append(abs(max(unc[i][0],key=abs)))
        unc[i][0].append(abs(max(unc[i][0],key=abs)))
        unc[i][1].append(abs(max(unc[i][1],key=abs)))
        total += max(unc[i][0],key=abs)**2
        totalVis += max(unc[i][1],key=abs)**2

    table = ""
    if takeMaxSyst:
        if twikiFormat:
            table += '| * Source *   |       *Ratio Change (%)*       |     *Vis Ratio Change (%)*     |\n'
        else:
            table += '\\begin{tabular}{l | c c }\n'
            table += '\\hline\n'
            table += 'Source & Ratio Change (\\%) & Vis Ratio Change (\\%) \\\\\n'
            table += '\\hline\n'
        for w in sorted(unc.items(),key=lambda x: x[1][0][2],reverse=True):
            if twikiFormat:
                table += "| %s  |  %.3f  |  %.3f  | \n" % (w[0], w[1][0][2], w[1][1][2])
            else:
                table += "%s   & %.3f & %.3f \\\\ \n" % (w[0], w[1][0][2], w[1][1][2])
    
    
        if not twikiFormat:
            table += '\\hline \n'
            table += 'Total  &  %.3f  &  %.3f  \\\\ \n' % (total**0.5, totalVis**0.5)
            table += '\\hline \n'
            table += '\\end{tabular} \n'
        else:
            table += '| Total  |  %.3f  |  %.3f  | \n' % (total**0.5, totalVis**0.5)

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
                table +=  "| %s  |  %.3f  |  %.3f  |  %.3f  |  %.3f  | \n" % (w[0], w[1][0][0], w[1][0][1],  w[1][1][0], w[1][1][1])
            else:
                table +=  "%s   & %.3f & %.3f & %.3f & %.3f \\\\ \n" % (w[0], w[1][0][0], w[1][0][1],  w[1][1][0], w[1][1][1])
    
    
        if not twikiFormat:
            table +=  '\\hline \n'
            table +=  '\\end{tabular} \n'
    return table, total**0.5, totalVis**0.5



if __name__=="__main__":
    
    print getTable(sys.argv)[0]
