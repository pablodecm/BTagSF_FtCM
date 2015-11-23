#!/usr/bin/env python

from mut_framework.LJMET_converter.samples_25ns import data_samples,mc_samples
import argparse
from ROOT import Builder, RooArgSet 
from ROOT import vector, RooRealVar, RooLognormal, RooProdPdf
from ROOT import RooFit, RooAddition, RooMinuit, RooArgList
import os

default_input_dir = "../output/13TeV_25ns/" 
default_output_dir = "./13TeV_25ns_inclusive/"

parser = argparse.ArgumentParser(description='NLL fit with FtCM.')
parser.add_argument('tagger', metavar='tag', type=str, 
                           help='the name of the tagger to perform the fit')
parser.add_argument('workpoint', metavar='wp', type=float,
                           help='the wp of the tagger to perform the fit')
parser.add_argument('--input_folder', default=default_input_dir, type=str)
parser.add_argument('--output_folder', default=default_output_dir, type=str)
args = parser.parse_args()


json_dir = args.input_folder 

data_names = ["SingleMuon_PRD",
              "SingleMuon_RRD"]


mc_names = [["TT_powheg",0],
            ["WJetsToLNu",1],
            ["DYJetsToLL_M-50", 1],
            ["ST_t-channel", 1],
            ["ST_s-channel", 1],
            ["ST_tW_top",1],
            ["ST_tW_antitop", 1],
            ["QCD_Pt-1000toInf_MuEnrichedPt5",1],
            ["QCD_Pt-120to170_MuEnrichedPt5",1],
            ["QCD_Pt-170to300_MuEnrichedPt5",1],
            ["QCD_Pt-300to470_MuEnrichedPt5",1],
            ["QCD_Pt-470to600_MuEnrichedPt5",1],
            ["QCD_Pt-50to80_MuEnrichedPt5",1],
            ["QCD_Pt-600to800_MuEnrichedPt5", 1],
            ["QCD_Pt-800to1000_MuEnrichedPt5",1],
            ["QCD_Pt-80to120_MuEnrichedPt5", 1],
            ["WW", 1],
            ["WZ", 1],
            ["ZZ", 1]]


data_comps = []
lumi = 0.0
for s in data_names:
    data_comps.append([json_dir+data_samples[s]["full_name"]+".json"])
    lumi += data_samples[s]["lumi"]
            
mc_comps = []
for s in mc_names:
    mc_comps.append([json_dir+mc_samples[s[0]]["full_name"]+".json",
                     mc_samples[s[0]]["gen_events"],
                     mc_samples[s[0]]["xs"],
                     s[1]])

result_dir =  args.output_folder
tagger = args.tagger
wp = args.workpoint


print "--- Fitting Run II 25ns Data ---"
print "  - tagger : {}".format(tagger)
print "  - workpoint : {}".format(wp)
print "--------------------------"


m = Builder(tagger, wp, lumi)
for mc_s in mc_comps:
    m.add_mc_component(*mc_s)
for data_s in data_comps:
    m.add_data_component(*data_s)

# builder options
m.zeroNegativeFracs = False
m.useDataPretagNorm = False 

if (m.useDataPretagNorm):
    m.kappa_.setConstant(True)

# min data counts per pretag category
cats = m.add_all_categories(True, 200)
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

