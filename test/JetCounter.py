#!/usr/bin/env python 

import ROOT
from ROOT import TChain, TCanvas, JetCounter, vector

rootfile = "../../../../../BetterTTrees/TTbar_Summer13.root"

# instance of selelector
jetCounter = JetCounter()

# taggers to analyze 
jetCounter.addTagger("combinedInclusiveSecondaryVertexBJetTags", 0.0, 1.0, 20)
jetCounter.addTagger("combinedMVABJetTags", 0.0, 1.0, 10)

# pt bins
ptBins = vector('double')()
ptBins += [0.0, 10000.0]
jetCounter.setPtBins(ptBins)
# eta bins
etaBins = vector('double')()
etaBins += [-2.4, 2.4]
jetCounter.setEtaBins(etaBins)

# tchain ("tree" is the TTree name)
tchain = TChain("tree")
tchain.Add( rootfile )

tchain.Process(jetCounter)

# test serialization of results to JSON
jetCounter.serialize("mc_output.json")

# test is JetRegistry can be written to a file
tfile = ROOT.TFile("test.root","RECREATE")
tfile.WriteObject(jetCounter.jetRegistry_,"jetRegistry")




