#! /usr/bin/env python
import os
import glob
import math

from optparse import OptionParser

parser = OptionParser()


############################################
#            Job steering                  #
############################################

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 
parser.add_option('--inputFiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='inputFiles',
                  help='Input files')

parser.add_option('--txtfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='txtfiles',
                  help='Input txt files')

parser.add_option("--onDcache", action='store_true',
                  default=True,
                  dest="onDcache",
                  help="onDcache(1), onDcache(0)")

# Output name to use. 
parser.add_option('--outputFile', metavar='F', type='string', action='store',
                  default='shyft_fwlite.root',
                  dest='outputFile',
                  help='output file name')

# Using MC info or not. For MC, truth information is accessed.
parser.add_option('--doMC', metavar='F', action='store_true',
                  default=False,
                  dest='doMC',
                  help='Check MC Information')

# Which lepton type to use
parser.add_option('--lepType', metavar='F', type='int', action='store',
                  default=0,
                  dest='lepType',
                  help='Lepton type. Options are 0 = muons, 1 = electrons')

# Invert MET cut?
parser.add_option('--invertMET', action='store_true',
                  default=False,
                  dest='invertMET',
                  help='Invert MET cut')

# Invert PF isolation cut?
parser.add_option('--invertPFIso', action='store_true',
                  default=False,
                  dest='invertPFIso',
                  help='Invert PF isolation cut')


(options, args) = parser.parse_args()

argv = []

# Import everything from ROOT
import ROOT
ROOT.gROOT.Macro("rootlogon.C")

# Import stuff from FWLite
import sys
from DataFormats.FWLite import Events, Handle

#infile = open( options.inputFiles )
#infileStr = infile.read().rstrip()

#print 'Getting files from this dir: ' + infileStr

# Get the file list. 
#files = glob.glob( infileStr )

# Get the file list.
if options.inputFiles:
    files = glob.glob( options.files )
    print 'getting files', files
elif options.txtfiles:
    files = []
    with open(options.txtfiles, 'r') as input_:
        for line in input_:
            files.append(line.strip())
else:
    files = []

print 'getting files: ', files

if options.onDcache:
        files = ["dcap://" + x for x in files]
        #print 'new files', *files, sep='\n'
        #print('new files', files[0], files[1], ..., sep='\n')

fname = options.txtfiles
fileN = fname[fname.rfind('/')+1:]

print files


# Create the output file. 
f = ROOT.TFile(options.outputFile, "recreate")
f.cd()


# Make histograms
print "Creating histograms"
secvtxMassHist = ROOT.TH1F('secvtxMassHist', "Secondary Vertex Mass", 150, 0., 5.0)
secvtxMassHistB = ROOT.TH1F('secvtxMassHistB', "Secondary Vertex Mass, b jets", 150, 0., 5.0)
secvtxMassHistC = ROOT.TH1F('secvtxMassHistC', "Secondary Vertex Mass, c jets", 150, 0., 5.0)
secvtxMassHistL = ROOT.TH1F('secvtxMassHistL', "Secondary Vertex Mass, udsg jets", 150, 0., 5.0)
metVsIso = ROOT.TH2F('metVsIso', 'MET Versus PFIsolation', 15, 0., 150., 125, 0., 2.5)
jetPtHist = ROOT.TH1F('jetPtHist', 'Jet p_{T}', 150, 0., 600.)
m3Hist = ROOT.TH1F('m3Hist', 'M3 Histogram', 150, 0., 600.)
muPtHist = ROOT.TH1F('muPtHist', 'Mu p_{T}', 150, 0., 600.)
njetHist = ROOT.TH1F('njetHist', 'Number of Jets', 10, 0., 10.)
HtHist = ROOT.TH1F('HtHist','H_{T}', 150, 0., 1500.)
METHist = ROOT.TH1F('METHist','MET', 150, 0., 600.)
WmtHist = ROOT.TH1F('WmtHist','W Transverse Mass', 150, 0., 600.)

############################################
# Physics level parameters for systematics #
############################################

# Kinematic cuts:
jetPtMin = 30.0
leadJetPtMin = 30.0
isoMax = 0.2
ssvheCut = 1.74
minJets = 4

if options.lepType == 0 :
    muonPtMin = 45.0
    electronPtMin = 20.0
    metMin = 20.0
    lepStr = 'Mu'
else:
    muonPtMin = 20.0
    electronPtMin = 35.0
    metMin = 20.0
    lepStr = 'Ele'




events = Events (files)

# Make the entirety of the handles required for the
# analysis. 
postfix = ""
if options.invertPFIso :
    postfix = "Loose"

puHandle         = Handle( "std::vector<float>" )
puLabel    = ( "PUNtupleDumper",   "PUweightNominalUpDown" )


jetPtHandle         = Handle( "std::vector<float>" )
jetPtLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "pt" )
jetEtaHandle         = Handle( "std::vector<float>" )
jetEtaLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "eta" )
jetPhiHandle         = Handle( "std::vector<float>" )
jetPhiLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "phi" )
jetMassHandle         = Handle( "std::vector<float>" )
jetMassLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "mass" )
jetSecvtxMassHandle         = Handle( "std::vector<float>" )
jetSecvtxMassLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "secvtxMass" )
jetSSVHEHandle         = Handle( "std::vector<float>" )
jetSSVHELabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "ssvhe" )
jetFlavorHandle         = Handle( "std::vector<float>" )
jetFlavorLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "flavor" )

muonPtHandle         = Handle( "std::vector<float>" )
muonPtLabel    = ( "pfShyftTupleMuons"+  postfix,   "pt" )
muonEtaHandle         = Handle( "std::vector<float>" )
muonEtaLabel    = ( "pfShyftTupleMuons"+  postfix,   "eta" )
muonPhiHandle         = Handle( "std::vector<float>" )
muonPhiLabel    = ( "pfShyftTupleMuons"+  postfix,   "phi" )
muonNhIsoHandle         = Handle( "std::vector<float>" )
muonNhIsoLabel    = ( "pfShyftTupleMuons"+  postfix,   "nhIso" )
muonChIsoHandle         = Handle( "std::vector<float>" )
muonChIsoLabel    = ( "pfShyftTupleMuons"+  postfix,   "chIso" )
muonPhIsoHandle         = Handle( "std::vector<float>" )
muonPhIsoLabel    = ( "pfShyftTupleMuons"+  postfix,   "phIso" )
muonPuIsoHandle         = Handle( "std::vector<float>" )
muonPuIsoLabel    = ( "pfShyftTupleMuons"+  postfix,   "puIso" )

electronPtHandle         = Handle( "std::vector<float>" )
electronPtLabel    = ( "pfShyftTupleElectrons"+  postfix,   "pt" )
electronEtaHandle         = Handle( "std::vector<float>" )
electronEtaLabel    = ( "pfShyftTupleElectrons"+  postfix,   "eta" )
electronPhiHandle         = Handle( "std::vector<float>" )
electronPhiLabel    = ( "pfShyftTupleElectrons"+  postfix,   "phi" )
electronNhIsoHandle         = Handle( "std::vector<float>" )
electronNhIsoLabel    = ( "pfShyftTupleElectrons"+  postfix,   "nhIso" )
electronChIsoHandle         = Handle( "std::vector<float>" )
electronChIsoLabel    = ( "pfShyftTupleElectrons"+  postfix,   "chIso" )
electronPhIsoHandle         = Handle( "std::vector<float>" )
electronPhIsoLabel    = ( "pfShyftTupleElectrons"+  postfix,   "phIso" )
electronPuIsoHandle         = Handle( "std::vector<float>" )
electronPuIsoLabel    = ( "pfShyftTupleElectrons"+  postfix,   "puIso" )


metHandle = Handle( "std::vector<float>" )
metLabel = ("pfShyftTupleMET" + lepStr +  postfix,   "pt" )
metPhiHandle = Handle( "std::vector<float>" )
metPhiLabel = ("pfShyftTupleMET" + lepStr +  postfix,   "phi" )


# Keep some timing information
nEventsAnalyzed = 0
nEventsPassed4Jets = 0
nEventsPassed1Tag = 0
timer = ROOT.TStopwatch()
timer.Start()

pairs = []

# loop over events
count = 0
ntotal = events.size()
percentDone = 0.0
ipercentDone = 0
ipercentDoneLast = -1
print "Start looping"
for event in events:
    nEventsAnalyzed += 1
    ipercentDone = int(percentDone)
    if ipercentDone != ipercentDoneLast :
        ipercentDoneLast = ipercentDone
        print 'Processing {0:10.0f}/{1:10.0f} : {2:5.0f}%'.format(
            count, ntotal, ipercentDone )
    count = count + 1
    percentDone = float(count) / float(ntotal) * 100.0

    #define and initialize Ht
    Ht = 0.0

    ################################################
    #   Retrieve the jet four-vector
    #   ------------------------------------
    #      The jet 4-vectors are large and hence
    #      take a long time to read out. If you don't
    #      need the other products (eta,phi,mass of jet)
    #      then don't read them out. 
    ################################################

    event.getByLabel( jetPtLabel, jetPtHandle )
    if not jetPtHandle.isValid():
        jetPts = None
    else :
        jetPts = jetPtHandle.product()

    if jetPts is None :
        continue
    event.getByLabel( jetEtaLabel, jetEtaHandle )
    jetEtas = jetEtaHandle.product()
    event.getByLabel( jetPhiLabel, jetPhiHandle )
    jetPhis = jetPhiHandle.product()
    event.getByLabel( jetMassLabel, jetMassHandle )
    jetMasses = jetMassHandle.product()


    # Find the njet bin we're in
    njets = -1
    event.getByLabel( jetPtLabel, jetPtHandle )
    jetPts = jetPtHandle.product()
    njets = 0
    for ijet in range( 0, len( jetPts ) ) :
        if jetPts[ijet] > 30.0 :
            njets += 1

    # We're not interested in <=4 jets
    if njets < minJets :
        continue
    nEventsPassed4Jets = nEventsPassed4Jets + 1



    ################################################
    #   Retrieve the jet vertex mass and plot the
    #   secondary vertex mass for tagged jets. 
    ################################################
    event.getByLabel (jetSSVHELabel, jetSSVHEHandle)
    jetSSVHEs = jetSSVHEHandle.product()
    
    # Now loop over the jets, and count tags
    ntags = 0
    for ijet in range(0,len(jetPts) ) :
        jetPt = jetPts[ijet]
        Ht += jetPt
        if jetPt < 30.0 :
            continue
        #jetEta = jetEtas[ijet]
        #jetPhi = jetPhis[ijet]
        #jetMass = jetMasses[ijet]
        jetSSVHE = jetSSVHEs[ijet]
        # plot secondary vertex mass for tagged jets
        if jetSSVHE >= ssvheCut :
            ntags = ntags + 1                    
    if ntags > 0:
        nEventsPassed1Tag = nEventsPassed1Tag + 1
        pairs.append( [event.object().id().run(),
                        event.object().id().luminosityBlock(),
                        event.object().id().event(),
                        njets,
                        ntags] )
    else :
        continue


    ################################################
    #   Require exactly one lepton (e or mu)
    #   ------------------------------------
    #      Our ntuples have both muon and electron
    #      events, and hence we must select events
    #      based on one or the other type. 
    #      To accomplish this we check the products
    #      for the type we're currently plotting
    #      (Mu or Ele), and check if the product is
    #      present. 
    ################################################
    muonPts = None
    electronPts = None
    if options.lepType == 0 :
        event.getByLabel (muonPtLabel, muonPtHandle)
        if not muonPtHandle.isValid():
            muonPts = None
        else :
            muonPts = muonPtHandle.product()

    elif options.lepType == 1 :
        event.getByLabel (electronPtLabel, electronPtHandle)
        if not electronPtHandle.isValid():
            electronPts = None
        else :
            electronPts = electronPtHandle.product()

    # If neither muons nor electrons are found, skip
    if muonPts is None and electronPts is None :
        continue
    # If we are looking for muons but none are found, skip
    if options.lepType == 0 and muonPts is None :
        continue
    # If we are looking for electrons but none are found, skip
    if options.lepType == 1 and electronPts is None :
        continue

    # keep leptons with certain pt threshold in the event
    if options.lepType == 0 and muonPts[0] <= muonPtMin:
        continue
    if options.lepType == 1 and electronPts[0] <= electronPtMin:
        continue
    
    #Now get muon eta
    event.getByLabel( muonEtaLabel, muonEtaHandle)
    muonEtas = muonEtaHandle.product()

    #Now get muon phi
    event.getByLabel( muonPhiLabel, muonPhiHandle)
    muonPhis = muonPhiHandle.product()

    #Now get the MET phi
    event.getByLabel( metPhiLabel, metPhiHandle)
    metPhi = metPhiHandle.product()[0]

    # Now get the MET
    event.getByLabel( metLabel, metHandle )
    metRaw = metHandle.product()[0]

    # Now get the PF isolation
    lepIso = -1.0
    if options.lepType == 0 and muonPts is not None :
        events.getByLabel( muonNhIsoLabel, muonNhIsoHandle )
        events.getByLabel( muonChIsoLabel, muonChIsoHandle )
        events.getByLabel( muonPhIsoLabel, muonPhIsoHandle )
        events.getByLabel( muonPuIsoLabel, muonPuIsoHandle )
        nhIso = muonNhIsoHandle.product()[0]
        chIso = muonChIsoHandle.product()[0]
        phIso = muonPhIsoHandle.product()[0] 
        puIso = muonPuIsoHandle.product()[0]  
        lepIso = (chIso + max(0.0, nhIso + phIso - 0.5*puIso)) / muonPts[0]
    if options.lepType == 1 and electronPts is not None :
        events.getByLabel( electronNhIsoLabel, electronNhIsoHandle )
        events.getByLabel( electronChIsoLabel, electronChIsoHandle )
        events.getByLabel( electronPhIsoLabel, electronPhIsoHandle )
        events.getByLabel( electronPuIsoLabel, electronPuIsoHandle )
        nhIso = electronNhIsoHandle.product()[0]
        chIso = electronChIsoHandle.product()[0]
        phIso = electronPhIsoHandle.product()[0]  
        puIso = electronPuIsoHandle.product()[0]
        lepIso = (chIso + max(0.0, nhIso + phIso - 0.5*puIso)) / electronPts[0]

    # Make a plot of the MET versus ISO for normalization purposes
    metVsIso.Fill( metRaw, lepIso )

    # If the MET is lower than our cut, skip, unless we want it inverted
    if not options.invertMET :
        if metRaw < metMin :
            continue
    else :
        if metRaw > metMin :
            continue

    # If the ISO is higher than our cut, skip, unless we want it inverted
    if not options.invertPFIso :
        if lepIso > isoMax :
            continue
    else :
        if lepIso < isoMax :
            continue


    event.getByLabel (jetSecvtxMassLabel, jetSecvtxMassHandle)
    jetSecvtxMasses = jetSecvtxMassHandle.product()
    if options.doMC :
        event.getByLabel( jetFlavorLabel, jetFlavorHandle )
        jetFlavors = jetFlavorHandle.product()
    # Now loop over the jets, and store the secondary vertex mass.
    ntags = 0
    jets_p4 = []
    for ijet in range(0,len(jetPts) ) :
        jetPt = jetPts[ijet]
        if jetPt < 30.0 :
            continue
        ijetP4 = ROOT.TLorentzVector()
        ijetP4.SetPtEtaPhiM( jetPts[ijet], jetEtas[ijet], jetPhis[ijet], jetMasses[ijet] )
        jets_p4.append( ijetP4 )
        jetSSVHE = jetSSVHEs[ijet]
        jetSecvtxMass = jetSecvtxMasses[ijet]
        jetPtHist.Fill( jetPt )
        # plot secondary vertex mass for tagged jets
        if jetSSVHE >= ssvheCut :
            ntags = ntags + 1
            secvtxMassHist.Fill( jetSecvtxMass )
            if options.doMC :
                if abs(jetFlavors[ijet]) == 5 :
                    secvtxMassHistB.Fill( jetSecvtxMass )
                elif abs(jetFlavors[ijet]) == 4 :
                    secvtxMassHistC.Fill( jetSecvtxMass )
                else :
                    secvtxMassHistL.Fill( jetSecvtxMass )

    # Now compute m3
    maxPt = -1.0
    m3 = -1.0
    for ijet in range(0, len(jets_p4) ) :
        for jjet in range(ijet + 1, len(jets_p4) ) :
            for kjet in range(jjet + 1, len(jets_p4) ) :
                sumP4 = jets_p4[ijet] + jets_p4[jjet] + jets_p4[kjet]
                if sumP4.Perp() > maxPt :
                    maxPt = sumP4.Perp()
                    m3 = sumP4.M()
    if maxPt > 0.0 :
        m3Hist.Fill( m3 )

    # Now plot mu pt
    muPtHist.Fill( muonPts[0] )
    
    #Now plot number of jets
    njetHist.Fill( njets )

    #Now plot Ht
    HtHist.Fill( Ht )

    #Now plot MET
    METHist.Fill(metRaw)

    #Now compute and plot W transerve mass
    #vpMu = ROOT.TLorentzVector()
    #vpMu.SetPtEtaPhiE( muonPts[0],muonEtas[0],muonPhis[0],muonPts[0])
    #vpNu = ROOT.TLorentzVector()
    #vpNu.SetPtEtaPhiE( metRaw, 0., metPhi, metRaw)
    #vpW = ROOT.TLorentzVector()
    #vpW = vpMu + vpNu
    #WmtHist.Fill(vpW.Mt())
    
    

# Stop our timer
timer.Stop()

# Print out our timing information
rtime = timer.RealTime(); # Real time (or "wall time")
ctime = timer.CpuTime(); # CPU time
print("Analyzed events: {0:6d}").format(nEventsAnalyzed)
print(">=4 jet events : {0:6d}").format(nEventsPassed4Jets)
print(">=1 tag events : {0:6d}").format(nEventsPassed1Tag)
print("RealTime={0:6.2f} seconds, CpuTime={1:6.2f} seconds").format(rtime,ctime)
print("{0:4.2f} events / RealTime second .").format( nEventsAnalyzed/rtime)
print("{0:4.2f} events / CpuTime second .").format( nEventsAnalyzed/ctime)



f.cd()
f.Write()
f.Close()

