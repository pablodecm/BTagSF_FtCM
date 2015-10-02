#!/usr/bin/env python

import argparse
from ROOT import KinFtCM, RooArgSet 
from ROOT import vector, RooRealVar, RooLognormal, RooProdPdf
from ROOT import RooFit, RooAddition, RooMinuit, RooArgList
import os

parser = argparse.ArgumentParser(description='NLL fit with FtCM.')
parser.add_argument('tagger', metavar='tag', type=str, 
                           help='the name of the tagger to perform the fit')
parser.add_argument('workpoint', metavar='wp', type=float,
                           help='the wp of the tagger to perform the fit')
args = parser.parse_args()

json_dir = "../output/kin_cat_runII/"

mc_samples = [[json_dir+"TTbar.json", 4992231, 831.76, 0 ],
              [json_dir+"WJets.json", 24089991, 61526.7, 1],
              [json_dir+"ZJets.json", 19925500, 6025.2, 1],
              [json_dir+"WW.json", 989608, 118.7, 1],
              [json_dir+"WZ.json", 996920, 66.1, 1],
              [json_dir+"ZZ.json", 998848, 19.3, 1],
              [json_dir+"T_t_top.json", 1273800, 136.02, 1],
              [json_dir+"T_t_atop.json", 681900, 80.95,1], 
              [json_dir+"T_tW_top.json", 998400, 35.6, 1],
              [json_dir+"T_tW_atop.json", 1000000, 35.6,1],
              [json_dir+"TTW.json", 1000000, 11.1, 1],
              [json_dir+"TTZll.json", 1000000, 11.1, 1],
              [json_dir+"TTZqq.json", 1000000, 11.1, 1] ]

data_samples = [[json_dir+"Single_mu_Prompt.json"],
                 [json_dir+"Single_mu_Rereco.json"]]

lumi = 41.6 

result_dir = "./KinFtCM_17092015/"
tagger = args.tagger
wp = args.workpoint


print "--- Fitting Run I Data ---"
print "  - tagger : {}".format(tagger)
print "  - workpoint : {}".format(wp)
print "--------------------------"


m = KinFtCM.Builder(tagger, wp, lumi)
for mc_s in mc_samples:
    m.add_mc_component(*mc_s)
for data_s in data_samples:
    m.add_data_component(*data_s)

# min data counts per pretag category
cats = m.add_all_categories(10)
m.set_mc_jet_tag_effs()

# get and fit data
data = m.get_data_hist()
#fit_result = m.sim_kin_pdf_.fitTo(data, RooFit.Extended(1), RooFit.NumCPU(8), RooFit.Save(1))
#fit_result.SaveAs(result_dir+"fit_result_{}_{}_first.root".format(tagger,wp))
#
## set mc effs and fit again
#m.b_jet_tag_effs_.Print("v")
#m.c_jet_tag_effs_.Print("v")
#m.l_jet_tag_effs_.Print("v")
#print "Set MC tag effs"
#m.set_mc_jet_tag_effs()
#m.b_jet_tag_effs_.Print("v")
#m.c_jet_tag_effs_.Print("v")
#m.l_jet_tag_effs_.Print("v")
#
# # fit data
#fit_result = m.sim_kin_pdf_.fitTo(data, RooFit.Extended(1), RooFit.NumCPU(8), RooFit.Save(1))
#fit_result.SaveAs(result_dir+"fit_result_{}_{}_second.root".format(tagger,wp))
#
## set mc effs and fit again
#m.b_jet_tag_effs_.Print("v")
#m.c_jet_tag_effs_.Print("v")
#m.l_jet_tag_effs_.Print("v")
#print "Set MC tag effs"
#m.set_mc_jet_tag_effs()
#m.b_jet_tag_effs_.Print("v")
#m.c_jet_tag_effs_.Print("v")
#m.l_jet_tag_effs_.Print("v")
#
#
#



