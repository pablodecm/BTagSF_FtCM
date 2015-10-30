
from mut_framework.BTagSF_FtCM.JetCounting import JetCountingManager
from ROOT import TFile


file_dir = "../../../../../BetterTTrees/13TeV_50ns/"

data_names = [ "SingleMuon_PR",
              "SingleMuon_RR",
              "SingleMuon_PRC_50ns"]

mc_names = ["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_50ns",
            "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_50ns",
            "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_50ns",
            "ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_50ns",
            "ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1_50ns",
            "ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1_50ns",
            "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_50ns",
            "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_50ns",
            "WW_TuneCUETP8M1_13TeV-pythia8_50ns",
            "WZ_TuneCUETP8M1_13TeV-pythia8_50ns" ]

mc_samples = [(file_dir + s + ".root") for s in mc_names]
data_samples = [(file_dir + s + ".root") for s in data_names]
ptBins = [30,10000]
etaBins = [-2.4, 2.4]

taggers = {"pfCombinedInclusiveSecondaryVertexV2BJetTags": [0.604, 0.890, 0.970]}

out_dir = "../output/13Tev_50ns_inclusive/" 

jetCM = JetCountingManager(ptBins, etaBins, taggers, oldMuonSF = False, out_dir=out_dir)  

jetCM.process(mc_samples, isData = False )
jetCM.process(data_samples, isData = True)
