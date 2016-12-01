#!/usr/bin/env python
import ROOT
from ROOT import *
import re
from array import array

c1 = TCanvas('c1', 'Plots', 1000, 500)
c1.SetFillColor(0)
c1.SetGrid()

h1 = []
f = ROOT.TFile("./sync.root")

t = f.Get("EventTree")
h1.Fill(
