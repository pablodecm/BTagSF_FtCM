
from mut_framework.BTagSF_FtCM.DataManager import Sample, DataManager

ttbar = Sample("../output/TTbar_Summer13.json", "TTbar", 6923750, 240.0)

print "PreTagging Efficiency: " + str(ttbar.pretag_eff())
tagger = "combinedMVABJetTags"
workPoint = 0.8
print "tagMultiplicity: " + str(ttbar.tag_multiplicity(tagger, workPoint))
cat_fractions = ttbar.cat_fractions() 
print "catFractions: " + str(cat_fractions)
print "sum catFractions: " + str(sum(cat_fractions.values()))
print "ttbar formula 0 btag: " + str(ttbar.formula(0))
print "ttbar formula 1 btag: " + str(ttbar.formula(1))


samples = [Sample("../output/TTbar_Summer13.json", "TTbar", 6923750, 240.0),
           Sample("../output/WJets_Summer13.json", "WJets", 76041475, 36257.2),
           Sample("../output/T_tW_Summer13.json", "T_tW", 497658, 11.15),
           Sample("../output/Tbar_tW_Summer13.json", "Tbar_tW", 493460, 11.15),
           Sample("../output/Tbar_t_Summer13.json", "Tbar_t", 1935072, 30.0),
           Sample("../output/ZJets_Summer13.json", "ZJets", 30459503, 3503.0),
           Sample("../output/T_t_Summer13.json", "T_t", 3758227, 55.531)]

lumi = 19789.0

dm = DataManager(samples)
dm.save_formula_list("../output/f_list.pickle")

