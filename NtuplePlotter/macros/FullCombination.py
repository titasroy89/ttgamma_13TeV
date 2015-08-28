#!/usr/bin/env python

import ROOT
from math import exp
import sys
import os

import likelihoodCombination

ROOT.gROOT.SetBatch()

takeMax = True
twikiFormat = False

#Set directories where muon and electron files are stored, looking for both the ratio.txt files and the template root files.
e_Directory = '/uscms_data/d2/dnoonan/TTGammaElectrons/NtuplePlotter/macros/ratioFiles_ele/'
mu_Directory = '/uscms_data/d2/dnoonan/TTGammaElectrons/NtuplePlotter/macros/ratioFiles_mu/'

templatesFileName = 'templates_barrel_scaled_afterPhotonM3.root'

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
            "EleEff",
            "elesmear",
            "MuEff",
            "musmear",
            ]


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
efficiencies = ['phoAcc'           ,
                'TTGamma_topEffAcc',
                'topPreselInt'     ,
                'TTJets_topEffAcc' ,
                'phoRecoEff'       ,
                'TTGammaVis_topAcc',
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
    combined_eff['TTgammaPhoEffAcc']   = mu_eff['TTGamma_topEffAcc']*mu_eff['phoAcc'] + e_eff['TTGamma_topEffAcc']*e_eff['phoAcc']
    combined_eff['TTJets_topEffAcc']   = mu_eff['TTJets_topEffAcc'] + e_eff['TTJets_topEffAcc']
    combined_eff['TTgammaPhoEffAccVis'] = mu_eff['TTGammaVis_topAcc']*mu_eff['phoRecoEff'] + e_eff['TTGammaVis_topAcc']*e_eff['phoRecoEff']
    
    combined_eff['topPreselErr'] = (mu_eff['topPreselErr']**2+e_eff['topPreselErr']**2)**0.5
    combined_eff['TTgammaPhoEffAccErr'] = ((mu_eff['TTGamma_topEffAccErr']/mu_eff['TTGamma_topEffAcc'])**2 + 
                                           (mu_eff['phoAccErr']/mu_eff['phoAcc'])**2 + 
                                           (e_eff['TTGamma_topEffAccErr']/e_eff['TTGamma_topEffAcc'])**2 + 
                                           (e_eff['phoAccErr']/e_eff['phoAcc'])**2)**0.5 * combined_eff['TTgammaPhoEffAcc']
    combined_eff['TTJets_topEffAccErr'] = (mu_eff['TTJets_topEffAccErr']**2 + e_eff['TTJets_topEffAccErr']**2)**0.5
    combined_eff['TTgammaPhoEffAccVisErr'] =  ((mu_eff['TTGammaVis_topAccErr']/mu_eff['TTGammaVis_topAcc'])**2 + 
                                               (mu_eff['phoRecoEffErr']/mu_eff['phoRecoEff'])**2 + 
                                               (e_eff['TTGammaVis_topAccErr']/e_eff['TTGammaVis_topAcc'])**2 + 
                                               (e_eff['phoRecoEffErr']/e_eff['phoRecoEff'])**2)**0.5 * combined_eff['TTgammaPhoEffAccVis']
    
    return e_data, mu_data, combined_eff

### dictionaries to store the ratios for the nominal and systematic sample combinations
ratioValues = {}
vis_ratioValues = {}

print 'Start'

### get the ratio for the nominal sample
### This is done last so theat the plots are saved
e_file = e_Directory + '/ratio_nominal.txt'
mu_file = mu_Directory + '/ratio_nominal.txt'
e_data, mu_data, combined_eff = findValues(e_file, mu_file)
likelihoodCombination.e_data = e_data
likelihoodCombination.mu_data = mu_data
result = likelihoodCombination.calculateTTGamma(e_Directory+'/'+templatesFileName.replace('.root','_nominal.root'), mu_Directory+'/'+templatesFileName.replace('.root','_nominal.root'), combined_eff, saveFitPlots = True, verbose=False)

print
print 'Nominal',  result
print combined_eff

ratioValues['nominal'] = result[0]
vis_ratioValues['nominal'] = result[1]


#loop over all systematics in systlist
for syst in systList:

    print syst
    # for all systematics, the cross section ratios are stored in ratioValues as a list, with the first value being the syst_down ratio and second being syst_up
    ratioValues[syst] = [0,0]
    vis_ratioValues[syst] = [0,0]
    
    e_file = e_Directory + '/ratio_'+syst+'_down.txt'
    mu_file = mu_Directory + '/ratio_'+syst+'_down.txt'
    if 'elesmear' in syst or 'EleEff' in syst:
        mu_file = mu_Directory + '/ratio_nominal.txt'
        
    elif 'musmear' in syst or 'MuEff' in syst:
        e_file = e_Directory + '/ratio_nominal.txt'

    e_data, mu_data, combined_eff = findValues(e_file, mu_file)
    likelihoodCombination.e_data = e_data
    likelihoodCombination.mu_data = mu_data
    if 'EleFake' in syst:
        likelihoodCombination.eleFakeSF = likelihoodCombination.eleFakeSF - likelihoodCombination.eleFakeSFErr
    else:
        likelihoodCombination.eleFakeSF = 1.5

    result = likelihoodCombination.calculateTTGamma(e_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_down.root'), mu_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_down.root'), combined_eff, saveFitPlots = False, verbose = False)

    print
    print 'Down'#, result

    ratioValues[syst][0] = result[0]
    vis_ratioValues[syst][0] = result[1]

    e_file = e_Directory + '/ratio_'+syst+'_up.txt'
    mu_file = mu_Directory + '/ratio_'+syst+'_up.txt'
    if 'elesmear' in syst or 'EleEff' in syst:
        mu_file = mu_Directory + '/ratio_nominal.txt'
    elif 'musmear' in syst or 'MuEff' in syst:
        e_file = e_Directory + '/ratio_nominal.txt'
    e_data, mu_data, combined_eff = findValues(e_file, mu_file)
    likelihoodCombination.e_data = e_data
    likelihoodCombination.mu_data = mu_data
    if 'EleFake' in syst:
        likelihoodCombination.eleFakeSF = likelihoodCombination.eleFakeSF + likelihoodCombination.eleFakeSFErr
    else:
        likelihoodCombination.eleFakeSF = 1.5

    result = likelihoodCombination.calculateTTGamma(e_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_up.root'), mu_Directory+'/'+templatesFileName.replace('.root','_'+syst+'_up.root'), combined_eff, saveFitPlots = False, verbose = False)
    print
    print 'Up'#, result
    ratioValues[syst][1] = result[0]
    vis_ratioValues[syst][1] = result[1]
    

nominalValue = ratioValues['nominal'][0]
unc = {'Nsignal':[ratioValues['nominal'][1]/nominalValue,ratioValues['nominal'][1]/nominalValue]}

vis_nominalValue = vis_ratioValues['nominal'][0]
visUnc = {'Nsignal':[vis_ratioValues['nominal'][1]/vis_nominalValue,vis_ratioValues['nominal'][1]/vis_nominalValue]}

for syst in systList:
    unc[syst] = [0,0]
    visUnc[syst] = [0,0]
    
    unc[syst][0] = (ratioValues[syst][0][0]-nominalValue)/nominalValue
    unc[syst][1] = (ratioValues[syst][1][0]-nominalValue)/nominalValue

    visUnc[syst][0] = (vis_ratioValues[syst][0][0]-vis_nominalValue)/vis_nominalValue
    visUnc[syst][1] = (vis_ratioValues[syst][1][0]-vis_nominalValue)/vis_nominalValue


print unc

print visUnc


total = 0.0
totalVis = 0.0

x = []

for i in unc:
    print i
    print unc[i]
    
#    unc[i].append(abs(max(unc[i][0],key=abs)))
    unc[i].append(abs(max(unc[i],key=abs)))
    visUnc[i].append(abs(max(visUnc[i],key=abs)))
    total += max(unc[i],key=abs)**2
    totalVis += max(visUnc[i],key=abs)**2
    
    x.append([i,unc[i],visUnc[i]])

#print x

x = sorted(x,key=lambda y: y[1][2],reverse=True)

#print x

print 'Total Unc = ', total**0.5
print 'Total Vis Unc = ', totalVis**0.5
print
print 'startTable'
if takeMax:
    if twikiFormat:
        print '| * Source *   |       *Ratio Change (%)*       |     *Vis Ratio Change (%)*     |'
    else:
        print '\\begin{tabular}{l | c c }'
        print '\\hline'
        print 'Source & Ratio Change (\\%) & Vis Ratio Change (\\%) \\\\'
        print '\\hline'
    for w in x:
        if twikiFormat:
            print "| %s  |  %.3f  |  %.3f  | " % (w[0], w[1][2], w[2][2])
        else:
            print "%s   & %.3f & %.3f \\\\" % (w[0], w[1][2], w[2][2])
    
    
    if not twikiFormat:
        print '\\hline'
        print 'Total  &  %.3f  &  %.3f  \\\\' % (total**0.5, totalVis**0.5)
        print '\\hline'
        print '\\end{tabular}'
    else:
        print '| Total  |  %.3f  |  %.3f  | ' % (total**0.5, totalVis**0.5)

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
