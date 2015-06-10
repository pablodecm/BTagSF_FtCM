#!/usr/bin/env python 

import ROOT
from ROOT import TChain, TCanvas, JetCounter, vector

rootfile = "../../../../../BetterTTrees/TTbar_Summer13.root"

# instance of selelector
jetCounter = JetCounter()

# taggers to analyze 
jetCounter.addTagger("combinedInclusiveSecondaryVertexBJetTags", 0.0, 1.0, 10)

# pt bins
ptBins = vector('double')()
ptBins += [0.0, 50.0, 100.0, 2000.0]
jetCounter.setPtBins(ptBins)
# eta bins
etaBins = vector('double')()
etaBins += [-2.4, 2.4]
jetCounter.setEtaBins(etaBins)

# tchain ("tree" is the TTree name)
tchain = TChain("tree")
tchain.Add( rootfile )

tchain.Process(jetCounter)

print "Pass events:  " + str(jetCounter.jetRegistry_.nEventPass_)

tcanvas = TCanvas()
jetCounter.all_good_jets.Draw("COLZ TEXT")
tcanvas.Print("good_jets.pdf")

tcanvas = TCanvas()
jetCounter.jetRegistry_.good_b_jets_.Draw("COLZ TEXT")
tcanvas.Print("good_b_jets.pdf")

# test is JetRegistry can be written to a file
tfile = ROOT.TFile("test.root","RECREATE")
tfile.WriteObject(jetCounter.jetRegistry_,"jetRegistry")
tfile.close()




