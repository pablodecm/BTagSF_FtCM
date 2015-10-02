
from mut_framework.BTagSF_FtCM.JetCounting import JetCountingManager
from ROOT import TFile


file_dir = "../../../../../BetterTTrees/13TeV/"

data = [ "Single_el_Prompt.root",
        "Single_el_Rereco.root",
        "Single_mu_Prompt.root",
        "Single_mu_Rereco.root"]

qcd = [ "QCD_1000to1400.root",
       "QCD_120to170.root",
       "QCD_1400to1800.root",
       "QCD_170to300.root",
       "QCD_1800to2400.root",
       "QCD_2400to3200.root",
       "QCD_300to470.root",
       "QCD_30to50.root",
       "QCD_3200toInf.root",
       "QCD_470to600.root",
       "QCD_50to80.root",
       "QCD_600to800.root",
       "QCD_800to1000.root",
       "QCD_80to120.root"]
top = ["T_t_atop.root",
       "TTbar.root",
       "T_t_top.root",
       "T_tW_atop.root",
       "TTW.root",
       "T_tW_top.root",
       "TTZll.root",
       "TTZqq.root"]
ewk = [ "WJets.root",
       "WW.root",
       "WZ.root",
       "ZJets.root",
       "ZZ.root"]

mc_samples = [(file_dir + s) for s in (qcd+top+ewk)]
data_samples = [(file_dir + s) for s in data]
ptBins = [30, 50, 80, 120, 160, 320]
etaBins = [-2.4, 2.4]

taggers = {"pfCombinedInclusiveSecondaryVertexV2BJetTags": [0.604, 0.890, 0.970]}

out_dir = "../output/kin_cat_runII/" 

jetCM = JetCountingManager(ptBins, etaBins, taggers, oldMuonSF = False, out_dir=out_dir)  

jetCM.process(mc_samples, isData = False )
jetCM.process(data_samples, isData = True)
