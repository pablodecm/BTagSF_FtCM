
from mut_framework.BTagSF_FtCM.JetCounting import JetCountingManager
from ROOT import TFile


file_dir = "../../../../../BetterTTrees/13TeV_25ns/"

data_names = ["SingleMuon_RRC_25ns",
              "SingleMuon_RRD_25ns",
              "SingleMuon_PRD_25ns"]

mc_names = ["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_25ns",
            "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_25ns",
            "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_25ns",
            "ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_25ns",
            "ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_25ns",
            "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_25ns",
            "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_25ns",
            "WW_TuneCUETP8M1_13TeV-pythia8_25ns",
            "WZ_TuneCUETP8M1_13TeV-pythia8_25ns",
            "TT_TuneCUETP8M1_13TeV-powheg-pythia8"]

qcd_names = ["QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_25ns",
             "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_25ns",
             "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_25ns",
             "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_25ns",
             "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_25ns",
             "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_25ns",
             "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_25ns",
             "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_25ns",
             "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_25ns"]


mc_samples = [(file_dir + s + ".root") for s in (mc_names + qcd_names)]
data_samples = [(file_dir + s + ".root") for s in data_names]
ptBins = [30, 10000]
etaBins = [-2.4, 2.4]

taggers = {"pfCombinedInclusiveSecondaryVertexV2BJetTags": [0.604, 0.890, 0.970]}

out_dir = "../output/13Tev_25ns_inclusive/" 

jetCM = JetCountingManager(ptBins, etaBins, taggers, oldMuonSF = False, out_dir=out_dir)  

jetCM.process(mc_samples, isData = False )
jetCM.process(data_samples, isData = True)
