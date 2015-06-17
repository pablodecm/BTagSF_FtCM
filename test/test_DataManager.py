
from mut_framework.BTagSF_FtCM.DataManager import DataManager

ttbar = DataManager("../output/TTbar_Summer13.json", "TTbar", 6923750, 240.0)

print "PreTagging Efficiency: " + str(ttbar.pretag_eff())
tagger = "combinedMVABJetTags"
workPoint = 0.8
print "tagMultiplicity: " + str(ttbar.tag_multiplicity(tagger, workPoint))
cat_fractions = ttbar.cat_fractions() 
print "catFractions: " + str(cat_fractions)
print "sum catFractions: " + str(sum(cat_fractions.values()))
print "ttbar formula 0 btag: " + str(ttbar.formula(0))
print "ttbar formula 1 btag: " + str(ttbar.formula(1))

