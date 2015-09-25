
from mut_framework.BTagSF_FtCM.JetCounting import JetCountingManager
from ROOT import TFile
import numpy as np


file_dir = "../../../../../BetterTTrees/"

data = ["Data_2012ABCD_Winter13_ReReco.root"]

mc_bn = ["TTbar_{}Summer13.root",
         "T_t_{}Summer13.root",
         "Tbar_t_{}Summer13.root",
         "T_tW_{}Summer13.root",
         "Tbar_tW_{}Summer13.root",
         "WJets_{}Summer13.root",
         "ZJets_{}Summer13.root"]

nom_samples = [s.format("") for s in mc_bn]
ttbar_s = [s.format("SUP_") for s in mc_bn[0:1]] + [s.format("SDOWN_") for s in mc_bn[0:1]]
ttbar_m= [s.format("MUP_") for s in mc_bn[0:1]] + [s.format("MDOWN_") for s in mc_bn[0:1]]
jec_samples = [s.format("JEC_UP_") for s in mc_bn] + [s.format("JEC_DOWN_") for s in mc_bn]
jer_samples = [s.format("JER_UP_") for s in mc_bn] + [s.format("JER_DOWN_") for s in mc_bn]

mc_samples = [(file_dir + s) for s in (nom_samples+ttbar_s+ttbar_m+jec_samples+jer_samples)]
data_samples = [(file_dir + s) for s in data]
ptBins = [0, 10000]
etaBins = [-2.4, 2.4]


wp_list = np.linspace(0.0,1.0,20).tolist() +  [0.244, 0.679, 0.898]
taggers = {"combinedSecondaryVertexBJetTags": wp_list}

out_dir = "../output/CSV_wp_scan/"

jetCM = JetCountingManager(ptBins, etaBins, taggers, oldMuonSF = True, out_dir = out_dir)  

jetCM.process(mc_samples, isData = False )
jetCM.process(data_samples, isData = True)