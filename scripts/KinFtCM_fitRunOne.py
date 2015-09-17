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

json_dir = "../output/kin_cat/"

mc_samples = [[json_dir+"TTbar_Summer13.json", 6923750, 240.0, 0 ],
              [json_dir+"WJets_Summer13.json", 76041475, 36257.2, 1],
              [json_dir+"T_tW_Summer13.json", 497658, 11.15, 1],
              [json_dir+"Tbar_tW_Summer13.json", 493460, 11.15, 1],
              [json_dir+"Tbar_t_Summer13.json", 1935072, 30.0, 1],
              [json_dir+"ZJets_Summer13.json", 30459503, 3503.0, 1],
              [json_dir+"T_t_Summer13.json", 3758227, 55.531,1]]
data_samples = [[json_dir+"Data_2012ABCD_Winter13_ReReco.json"]]

lumi = 19789.0

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
m.add_all_categories(100)
m.set_mc_jet_tag_effs()

# get and fit data
data = m.get_data_hist()
fit_result = m.sim_kin_pdf_.fitTo(data, RooFit.Extended(1), RooFit.NumCPU(8), RooFit.Save(1))
fit_result.SaveAs(result_dir+"fit_result_{}_{}_first.root".format(tagger,wp))

# set mc effs and fit again
m.b_jet_tag_effs_.Print("v")
m.c_jet_tag_effs_.Print("v")
m.l_jet_tag_effs_.Print("v")
print "Set MC tag effs"
m.set_mc_jet_tag_effs()
m.b_jet_tag_effs_.Print("v")
m.c_jet_tag_effs_.Print("v")
m.l_jet_tag_effs_.Print("v")

 # fit data
fit_result = m.sim_kin_pdf_.fitTo(data, RooFit.Extended(1), RooFit.NumCPU(8), RooFit.Save(1))
fit_result.SaveAs(result_dir+"fit_result_{}_{}_second.root".format(tagger,wp))

# set mc effs and fit again
m.b_jet_tag_effs_.Print("v")
m.c_jet_tag_effs_.Print("v")
m.l_jet_tag_effs_.Print("v")
print "Set MC tag effs"
m.set_mc_jet_tag_effs()
m.b_jet_tag_effs_.Print("v")
m.c_jet_tag_effs_.Print("v")
m.l_jet_tag_effs_.Print("v")






