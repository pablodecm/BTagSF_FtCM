from ROOT import Model, ModelPdf 
from ROOT import RooArgSet 
from ROOT import RooUniform, RooSimultaneous, RooExtendPdf 
from ROOT import RooFit
from ROOT import vector

mc_samples = [["../output/TTbar_Summer13.json", 6923750, 240.0, 0 ],
              ["../output/WJets_Summer13.json", 76041475, 36257.2, 1],
              ["../output/T_tW_Summer13.json", 497658, 11.15, 1],
              ["../output/Tbar_tW_Summer13.json", 493460, 11.15, 1],
              ["../output/Tbar_t_Summer13.json", 1935072, 30.0, 1],
              ["../output/ZJets_Summer13.json", 30459503, 3503.0, 1],
              ["../output/T_t_Summer13.json", 3758227, 55.531,1]]
data_samples = [["../output/Data_2012ABCD_Winter13_ReReco.json"]]

lumi = 19789.0
max_n_tags = 7 

m = Model(lumi)
for mc_s in mc_samples:
    m.add_mc_component(*mc_s)
for data_s in data_samples:
    m.add_data_component(*data_s)

cat_mapping = vector(vector("int"))()    
for i in range(4):
    cat_mapping.push_back(vector("int")())
    cat_mapping[i] += [i]

m.set_category_mapping(cat_mapping)
m.set_tag_wp("combinedMVABJetTags", 0.4)
for i in range(3):
    m.tag_effs_[i].setConstant()
m.set_pdfs(max_n_tags)
data = m.get_data_hist(max_n_tags)
m.sim_pdf_.fitTo(data, RooFit.Extended()) 
