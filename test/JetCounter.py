#!/usr/bin/env python 

import ROOT
from ROOT import TChain, TCanvas, JetCounter, vector

rootfile = "../../../../../BetterTTrees/TTbar_Summer13.root"

# instance of selelector
jetCounter = JetCounter()

# taggers to analyze 
jetCounter.addTagger("combinedInclusiveSecondaryVertexBJetTags", 0.0, 1.0, 100)

# pt bins
ptBins = vector('double')()
ptBins += [0.0, 100.0, 200.0]
jetCounter.setPtBins(ptBins)
# eta bins
etaBins = vector('double')()
etaBins += [-2.4, 2.4]
jetCounter.setEtaBins(etaBins)

# tchain ("tree" is the TTree name)
tchain = TChain("tree")
tchain.Add( rootfile )

tchain.Process(jetCounter)

tcanvas = TCanvas()
jetCounter.all_jets.Draw("COLZ TEXT")
tcanvas.Print("all_jets.pdf")



