
from mut_framework.BTagSF_FtCM.JetCounting import JetCountingManager

mc_samples = [ "../../../../../BetterTTrees/TTbar_Summer13.root",
               "../../../../../BetterTTrees/T_t_Summer13.root",
               "../../../../../BetterTTrees/T_tW_Summer13.root",
               "../../../../../BetterTTrees/WJets_Summer13.root",
               "../../../../../BetterTTrees/ZJets_Summer13.root",
               "../../../../../BetterTTrees/Tbar_t_Summer13.root",
               "../../../../../BetterTTrees/Tbar_tW_Summer13.root"]

data_samples = ["../../../../../BetterTTrees/Data_2012ABCD_Winter13_ReReco.root"]
ptBins = [0.0, 100000.0]
etaBins = [-2.4, 2.4]

taggers = {"combinedInclusiveSecondaryVertexBJetTags": [0.3, 0.5, 0.7],
           "combinedMVABJetTags": [0.4, 0.8]} 


jetCM = JetCountingManager(ptBins, etaBins, taggers)  

jetCM.process(mc_samples, isData = False)
jetCM.process(data_samples, isData = True)
