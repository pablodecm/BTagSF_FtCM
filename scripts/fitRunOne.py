#!/usr/bin/env python

import argparse
from ROOT import Model, ModelPdf 
from ROOT import vector
from ROOT import RooFit, RooAddition, RooMinuit, RooArgList

parser = argparse.ArgumentParser(description='NLL fit with FtCM.')
parser.add_argument('tagger', metavar='tag', type=str, 
                           help='the name of the tagger to perform the fit')
parser.add_argument('workpoint', metavar='wp', type=float,
                           help='the wp of the tagger to perform the fit')
args = parser.parse_args()


mc_samples = [["../output/TTbar_Summer13.json", 6923750, 240.0, 0 ],
              ["../output/WJets_Summer13.json", 76041475, 36257.2, 1],
              ["../output/T_tW_Summer13.json", 497658, 11.15, 1],
              ["../output/Tbar_tW_Summer13.json", 493460, 11.15, 1],
              ["../output/Tbar_t_Summer13.json", 1935072, 30.0, 1],
              ["../output/ZJets_Summer13.json", 30459503, 3503.0, 1],
              ["../output/T_t_Summer13.json", 3758227, 55.531,1]]
data_samples = [["../output/Data_2012ABCD_Winter13_ReReco.json"]]

lumi = 19789.0
max_n_tags = 6 

tagger = args.tagger
wp = args.workpoint
result_dir = "./"

print "--- Fitting Run I Data ---"
print "  - tagger : {}".format(tagger)
print "  - workpoint : {}".format(wp)

m = Model(lumi)
for mc_s in mc_samples:
    m.add_mc_component(*mc_s)
for data_s in data_samples:
    m.add_data_component(*data_s)

m.set_tag_wp(tagger, wp)
m.set_pdfs(max_n_tags)

h_ftcm = m.get_data_hist(max_n_tags)
h_kin = m.get_data_kin_hist()
nll_ftcm = m.sim_pdf_.createNLL(h_ftcm, RooFit.Extended(), RooFit.NumCPU(6))
nll_kin = m.sim_kin_pdf_.createNLL(h_kin, RooFit.Extended(), RooFit.NumCPU(4))
nll = RooAddition("nll","nll",RooArgList(nll_ftcm,nll_kin))
minuit = RooMinuit(nll)
fit_result = minuit.fit("r")

fit_result.SaveAs("fit_result_{}_{}.root".format(tagger,wp))


