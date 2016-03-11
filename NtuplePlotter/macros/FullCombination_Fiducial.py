#!/usr/bin/env python

import ROOT
from math import exp
import sys
import os

import likelihoodCombination_Fiducial_Quick as likelihoodCombination_Fiducial

ROOT.gROOT.SetBatch()

def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])
    else:
        return([])


takeMax = True
twikiFormat = False

#Set directories where muon and electron files are stored, looking for both the ratio.txt files and the template root files.
e_Directory = '/uscms_data/d2/dnoonan/TTGammaElectrons/NtuplePlotter/macros/ratioFilesFiducial_ele/'
mu_Directory = '/uscms_data/d2/dnoonan/TTGammaElectrons/NtuplePlotter/macros/ratioFilesFiducial_mu/'

templatesFileName = 'templates_barrel_scaled.root'#_afterPhotonM3.root'

#list of systematics to look for
systList = ["Btag",
            "EleFakeSF",
            "JEC",
            "JER",
            "PU",
            "QCD",
            "otherMC",
            "ZJetsSF",
            "pho",
            "toppt",
            "TopMass",
            "Scale",
            "Matching",
            "EleEff",
            "elesmear",
            "MuEff",
            "musmear",
            ]

# systList = ["TopMass",
#             "Btag",
#             # "Scale",
#             # "Matching",
#             ]

#systList = []

#check if all files are there before starting the likelihood fits
for _dir in [e_Directory, mu_Directory]:
    if not os.path.exists(_dir):
        print 'Path does not exist:', _dir
        exit()
    if not os.path.exists(_dir+'/ratio_nominal.txt'):
        print 'File does not exist:', _dir+'/ratio_nominal.txt'
        exit()
    if not os.path.exists(_dir+'/'+templatesFileName.replace('.root','_nominal.root')):
        print 'File does not exist:', _dir+'/'+templatesFileName.replace('.root','_nominal.root')
        exit()
    for syst in systList[:-4]:
        if not os.path.exists(_dir+'/ratio_'+syst+'_up.txt'):
            print 'File does not exist:', _dir+'/ratio_'+syst+'_up.txt'
            exit()
        if not os.path.exists(_dir+'/'+templatesFileName.replace('.root','_'+syst+'_up.root')):
            print 'File does not exist:', _dir+'/'+templatesFileName.replace('.root','_'+syst+'_up.root')
            exit()
        if not os.path.exists(_dir+'/ratio_'+syst+'_down.txt'):
            print 'File does not exist:', _dir+'/ratio_'+syst+'_down.txt'
            exit()
        if not os.path.exists(_dir+'/'+templatesFileName.replace('.root','_'+syst+'_down.root')):
            print 'File does not exist:', _dir+'/'+templatesFileName.replace('.root','_'+syst+'_down.root')
            exit()
for syst in systList[-4:-2]:
    _dir = e_Directory
    if not os.path.exists(_dir+'/ratio_'+syst+'_up.txt'):
        print 'File does not exist:', _dir+'/ratio_'+syst+'_up.txt'
        exit()
    if not os.path.exists(_dir+'/'+templatesFileName.replace('.root','_'+syst+'_up.root')):
        print 'File does not exist:', _dir+'/'+templatesFileName.replace('.root','_'+syst+'_up.root')
        exit()
    if not os.path.exists(_dir+'/ratio_'+syst+'_down.txt'):
        print 'File does not exist:', _dir+'/ratio_'+syst+'_down.txt'
        exit()
    if not os.path.exists(_dir+'/'+templatesFileName.replace('.root','_'+syst+'_down.root')):
        print 'File does not exist:', _dir+'/'+templatesFileName.replace('.root','_'+syst+'_down.root')
        exit()
for syst in systList[-2:]:
    _dir = mu_Directory
    if not os.path.exists(_dir+'/ratio_'+syst+'_up.txt'):
        print 'File does not exist:', _dir+'/ratio_'+syst+'_up.txt'
        exit()
    if not os.path.exists(_dir+'/'+templatesFileName.replace('.root','_'+syst+'_up.root')):
        print 'File does not exist:', _dir+'/'+templatesFileName.replace('.root','_'+syst+'_up.root')
        exit()
    if not os.path.exists(_dir+'/ratio_'+syst+'_down.txt'):
        print 'File does not exist:', _dir+'/ratio_'+syst+'_down.txt'
        exit()
    if not os.path.exists(_dir+'/'+templatesFileName.replace('.root','_'+syst+'_down.root')):
        print 'File does not exist:', _dir+'/'+templatesFileName.replace('.root','_'+syst+'_down.root')
        exit()


print 'All files are present'

#list of the parameters (and errors) to look for in the ratio files
# these names are then used as keys for the dictionaries storing the values
parameters = ['photnPurity'          ,
              'M3_photon_topFrac'    ,
              'Ndata'                ,
              'photnPurityErr'       ,
              'M3_photon_topFracErr' ,
              'NdataErr'             ,
              ]

#list of the efficiencies (and errors) to look for in the ratio files
# these names are then used as keys for the dictionaries storing the values
efficiencies = ['FidTopEff',
                'FidPhoEff',
                'FidAcc',
                'nPreselSemiLepChannel',
                'nGeneratedSemiLepChannel',
                'phoAcc'           ,
                'TTGamma_topEffAcc',
                'topPreselInt'     ,
                'TTJets_topEffAcc' ,
                'phoRecoEff'       ,
                'TTGammaVis_topAcc',
                'FidTopEffErr',
                'FidPhoEffErr',
                'FidAccErr',
                'nPreselSemiLepChannelErr',
                'nGeneratedSemiLepChannelErr',
                'phoAccErr'           ,
                'TTGamma_topEffAccErr',
                'topPreselErr'        ,
                'TTJets_topEffAccErr' ,
                'phoRecoEffErr'       ,
                'TTGammaVis_topAccErr'
                ]



def findValues(e_fileName, mu_fileName):
    """
    Parses the ratio.txt files from the output of makePlots to get the parameter and efficiency values needed for the likelihood fit and cross section ratio calculation
    Returns three dictionaries:
       e_data: the parameters and uncertainties from the electron+jets channel
       mu_data: the parameters and uncertainties from the muon+jets channel
       combined_eff: the combined efficiency values for the two channels
    """

    e_ratioFile = file(e_fileName,'r')
    e_data = {}
    e_eff = {}
    
    for line in e_ratioFile:
        for key in parameters:
            if key in line: 
                if 'Err' in key and 'Err' in line:
                    e_data[key]=float(line.split()[-1])
                if not 'Err' in key and not 'Err' in line:
                    e_data[key]=float(line.split()[-1])
        for key in efficiencies:
            if key in line: 
                if 'Err' in key and 'Err' in line:
                    e_eff[key]=float(line.split()[-1])
                if not 'Err' in key and not 'Err' in line:
                    e_eff[key]=float(line.split()[-1])

    mu_ratioFile = file(mu_fileName,'r')
    
    mu_data = {}
    mu_eff = {}
    
    for line in mu_ratioFile:
        for key in parameters:
            if key in line: 
                if 'Err' in key and 'Err' in line:
                    mu_data[key]=float(line.split()[-1])
                if not 'Err' in key and not 'Err' in line:
                    mu_data[key]=float(line.split()[-1])
        for key in efficiencies:
            if key in line: 
                if 'Err' in key and 'Err' in line:
                    mu_eff[key]=float(line.split()[-1])
                if not 'Err' in key and not 'Err' in line:
                    mu_eff[key]=float(line.split()[-1])


    
                    

    combined_eff = {}

    combined_eff['topPreselInt']       = mu_eff['topPreselInt'] + e_eff['topPreselInt']
    combined_eff['topPreselErr'] = (mu_eff['topPreselErr']**2+e_eff['topPreselErr']**2)**0.5

    combined_eff['FidEff'] = (mu_eff['nPreselSemiLepChannel']*mu_eff['FidPhoEff'] + e_eff['nPreselSemiLepChannel']*e_eff['FidPhoEff'])/(mu_eff['nGeneratedSemiLepChannel'] + e_eff['nGeneratedSemiLepChannel'])

    combined_eff['FidEffErr'] = ( (mu_eff['nPreselSemiLepChannelErr']*mu_eff['FidPhoEff'])**2 + (mu_eff['nPreselSemiLepChannel']*mu_eff['FidPhoEffErr'])**2 + (e_eff['nPreselSemiLepChannelErr']*e_eff['FidPhoEff'])**2 + (e_eff['nPreselSemiLepChannel']*e_eff['FidPhoEffErr'])**2)**0.5/(mu_eff['nGeneratedSemiLepChannel'] + e_eff['nGeneratedSemiLepChannel'])
 
#    combined_eff['FidEffErr'] = (mu_eff['nPreselSemiLepChannel']*mu_eff['FidPhoEff'] + mu_eff['nPreselSemiLepChannel']*mu_eff['FidPhoEff'])/(mu_eff['nGeneratedSemiLepChannel'] + e_eff['nGeneratedSemiLepChannel'])

    combined_eff['TTJets_topEffAcc']   = mu_eff['TTJets_topEffAcc'] + e_eff['TTJets_topEffAcc']
    combined_eff['TTJets_topEffAccErr'] = (mu_eff['TTJets_topEffAccErr']**2 + e_eff['TTJets_topEffAccErr']**2)**0.5

    combined_eff['FidAcc'] = mu_eff['FidAcc'] + e_eff['FidAcc']
    combined_eff['FidAccErr'] = ( mu_eff['FidAccErr']**2 + e_eff['FidAccErr']**2 )**0.5
    
    return e_data, mu_data, combined_eff, e_eff, mu_eff

### dictionaries to store the ratios for the nominal and systematic sample combinations
ratioValues = {}
dirXSValues = {}
print 'Start'

### get the ratio for the nominal sample
### This is done last so theat the plots are saved
e_file = e_Directory + '/ratio_nominal.txt'
mu_file = mu_Directory + '/ratio_nominal.txt'
e_data, mu_data, combined_eff, e_eff, mu_eff = findValues(e_file, mu_file)
nominal_combined_eff = combined_eff
likelihoodCombination_Fiducial.e_data = e_data
likelihoodCombination_Fiducial.mu_data = mu_data
result = likelihoodCombination_Fiducial.calculateTTGamma(e_Directory+'/'+templatesFileName.replace('.root','_nominal.root'), mu_Directory+'/'+templatesFileName.replace('.root','_nominal.root'), combined_eff, saveFitPlots = True, verbose=True)

print
print 'Nominal', result
print combined_eff
# print mu_eff
# print e_eff
print 
print mu_data
print e_data

ratioValues['nominal'] = result[0]
dirXSValues['nominal'] = result[1]


# print 'Check Inverted Selection'
# e_file = e_Directory + '/ratio_nominalInverted.txt'
# mu_file = mu_Directory + '/ratio_nominalInverted.txt'
# # e_data, mu_data, combined_eff, e_eff, mu_eff = findValues(e_file, mu_file)
# # combined_eff['FidEff'] = nominal_combined_eff['FidEff']
# likelihoodCombination_Fiducial.e_data = e_data
# likelihoodCombination_Fiducial.mu_data = mu_data
# likelihoodCombination_Fiducial.ttgammaSeq = seq(0.2, 0.7, 0.01) + seq(0.705,1.15,0.005) + seq(1.16, 1.6,0.01) + seq(1.62,3,.02)
# result = likelihoodCombination_Fiducial.calculateTTGamma(e_Directory+'/'+templatesFileName.replace('.root','_nominalInverted.root'), mu_Directory+'/'+templatesFileName.replace('.root','_nominalInverted.root'), combined_eff, saveFitPlots = False, verbose=False)

# print 'Inverted Sample', result
# print combined_eff
# print 
# print mu_data
# print e_data

# ratioValues['nominalInverted'] = result[0]
# dirXSValues['nominalInverted'] = result[1]

#loop over all systematics in systlist
for syst in systList:

    print syst
    # for all systematics, the cross section ratios are stored in ratioValues as a list, with the first value being the syst_down ratio and second being syst_up
    ratioValues[syst] = [0,0]
    dirXSValues[syst] = [0,0]
    
    e_file = e_Directory + '/ratio_'+syst+'_down.txt'
    mu_file = mu_Directory + '/ratio_'+syst+'_down.txt'

    e_templateFile = e_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_down.root')
    mu_templateFile = mu_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_down.root')
    if 'elesmear' in syst or 'EleEff' in syst:
        mu_file = mu_Directory + '/ratio_nominal.txt'
        mu_templateFile = mu_Directory+'/'+templatesFileName.replace('.root','_nominal.root')        
    elif 'musmear' in syst or 'MuEff' in syst:
        e_file = e_Directory + '/ratio_nominal.txt'
        e_templateFile = e_Directory+'/'+templatesFileName.replace('.root','_nominal.root')


    print e_file, mu_file
    e_data, mu_data, combined_eff, e_eff, mu_eff = findValues(e_file, mu_file)
    likelihoodCombination_Fiducial.e_data = e_data
    likelihoodCombination_Fiducial.mu_data = mu_data

    if 'EleFake' in syst:
        likelihoodCombination_Fiducial.eleFakeSF = likelihoodCombination_Fiducial.eleFakeSF - likelihoodCombination_Fiducial.eleFakeSFErr
    else:
        likelihoodCombination_Fiducial.eleFakeSF = 1.458

    if syst in ['TopMass', 'Scale', 'Matching']:
        likelihoodCombination_Fiducial.ttgammaSeq = seq(0.705,1.15,0.005) + seq(1.16, 1.6,0.01) + seq(1.62,3,.02)
        combined_eff['FidEff'] = nominal_combined_eff['FidEff']
    else:
        likelihoodCombination_Fiducial.ttgammaSeq = seq(0.2, 0.7, 0.01) + seq(0.705,1.15,0.005) + seq(1.16, 1.6,0.01)

#    result = likelihoodCombination_Fiducial.calculateTTGamma(e_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_down.root'), mu_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_down.root'), combined_eff, saveFitPlots = False, verbose = False)
    result = likelihoodCombination_Fiducial.calculateTTGamma(e_templateFile, mu_templateFile, combined_eff, saveFitPlots = False, verbose = True)
    print
    print syst, 'Down'#, result
    print result
    print combined_eff
    print 
    print mu_data
    print e_data

    ratioValues[syst][0] = result[0]
    dirXSValues[syst][0] = result[1]

    e_file = e_Directory + '/ratio_'+syst+'_up.txt'
    mu_file = mu_Directory + '/ratio_'+syst+'_up.txt'
    e_templateFile = e_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_up.root')
    mu_templateFile = mu_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_up.root')
    if 'elesmear' in syst or 'EleEff' in syst:
        mu_file = mu_Directory + '/ratio_nominal.txt'
        mu_templateFile = mu_Directory+'/'+templatesFileName.replace('.root','_nominal.root')
    elif 'musmear' in syst or 'MuEff' in syst:
        e_file = e_Directory + '/ratio_nominal.txt'
        e_templateFile = e_Directory+'/'+templatesFileName.replace('.root','_nominal.root')
    e_data, mu_data, combined_eff, e_eff, mu_eff = findValues(e_file, mu_file)
    likelihoodCombination_Fiducial.e_data = e_data
    likelihoodCombination_Fiducial.mu_data = mu_data
    if 'EleFake' in syst:
        likelihoodCombination_Fiducial.eleFakeSF = likelihoodCombination_Fiducial.eleFakeSF + likelihoodCombination_Fiducial.eleFakeSFErr
    else:
        likelihoodCombination_Fiducial.eleFakeSF = 1.458

    if syst in ['TopMass', 'Scale', 'Matching']:
        combined_eff['FidEff'] = nominal_combined_eff['FidEff']

#    result = likelihoodCombination_Fiducial.calculateTTGamma(e_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_up.root'), mu_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_up.root'), combined_eff, saveFitPlots = False, verbose = False)
    result = likelihoodCombination_Fiducial.calculateTTGamma(e_templateFile, mu_templateFile, combined_eff, saveFitPlots = False, verbose = True)
    print
    print syst, 'Up'#, result
    print result
    print combined_eff
    print 
    print mu_data
    print e_data

    ratioValues[syst][1] = result[0]
    dirXSValues[syst][1] = result[1]
    

nominalValue = ratioValues['nominal'][0]
unc = {'Nsignal':[ratioValues['nominal'][1]/nominalValue,ratioValues['nominal'][1]/nominalValue]}

dirNomValue = dirXSValues['nominal'][0]
directUnc = {'Nsignal':[dirXSValues['nominal'][1]/dirNomValue,dirXSValues['nominal'][1]/dirNomValue]}

for syst in systList:
    unc[syst] = [0,0]
    directUnc[syst] = [0,0]

    unc[syst][0] = (ratioValues[syst][0][0]-nominalValue)/nominalValue
    unc[syst][1] = (ratioValues[syst][1][0]-nominalValue)/nominalValue

    directUnc[syst][0] = (dirXSValues[syst][0][0]-dirNomValue)/dirNomValue
    directUnc[syst][1] = (dirXSValues[syst][1][0]-dirNomValue)/dirNomValue


print unc


total = 0.0

x = []


print 'Start Systematics Dump'
for i in unc:
    print i
    print unc[i]
    
#    unc[i].append(abs(max(unc[i][0],key=abs)))
    unc[i].append(abs(max(unc[i],key=abs)))
    total += max(unc[i],key=abs)**2
    
    x.append([i,unc[i]])

print 'End Systematics Dump'

dirTotal = 0.0
print '-'*25
print 'Direct XS =',dirNomValue
print '-'*25
for i in unc:
    print i, directUnc[i]
print '-'*25

#print x

x = sorted(x,key=lambda y: y[1][2],reverse=True)

#print x

print 'Total Unc = ', total**0.5
print
print 'startTable'
if takeMax:
    if twikiFormat:
        print '| * Source *   |       *Ratio Change (%)*       |'
    else:
        print '\\begin{tabular}{l | c }'
        print '\\hline'
        print 'Source & Ratio Change (\\%) \\\\'
        print '\\hline'
    for w in x:
        if twikiFormat:
            print "| %s  |  %.3f  | " % (w[0], w[1][2])
        else:
            print "%s   & %.3f  \\\\" % (w[0], w[1][2] )
    
    
    if not twikiFormat:
        print '\\hline'
        print 'Total  &  %.3f \\\\' % (total**0.5)
        print '\\hline'
        print '\\end{tabular}'
    else:
        print '| Total  |  %.3f  |  %.3f  | ' % (total**0.5)

else:
    if twikiFormat:
        print '|             |       *Ratio Change (%)*       ||     *Vis Ratio Change (%)*     ||'
        print '| * Source *  |  * Syst down *  |  * Syst up *  |  * Syst down *  |  * Syst up *  |'
    else:
        print '\\begin{tabular}{l | c c | c c}'
        print '\\hline'
        print '& \\multicolumn{2}{c}{Ratio Change (\\%)} & \\multicolumn{2}{c}{Vis Ratio Change (\\%)} \\\\'
        print 'Source & Syst Down & Syst Up & Syst Down & Syst Up \\\\'
        print '\\hline'
    for w in x:
        if twikiFormat:
            print "| %s  |  %.3f  |  %.3f  |  %.3f  |  %.3f  |" % (w[0], w[1][0], w[1][1],  w[2][0], w[2][1])
        else:
            print "%s   & %.3f & %.3f & %.3f & %.3f \\\\" % (w[0], w[1][0], w[1][1],  w[2][0], w[2][1])
    
    
    if not twikiFormat:
        print '\\hline'
        print '\\end{tabular}'
print 'stopTable'
