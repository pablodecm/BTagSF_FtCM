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

lumi = 71.52 


json_dir = "../output/13Tev_50ns_inclusive/" 

data_names = ["SingleMuon_PR",
              "SingleMuon_RR",
              "SingleMuon_PRC_50ns"]

mc_names = [["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_50ns",4994250*0.331907,831.76,0], 
            ["WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_50ns",24089991*0.683927, 61526.7,1],
            ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_50ns", 19925500*0.670032, 6025.2, 1],
            ["ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_50ns",3999910*0.215131,70.69,1],
            ["ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_50ns", 998400, 35.6,1],
            ["ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_50ns",1000000, 35.6, 1]]

data_samples = []
for s in data_names:
    data_samples.append([json_dir+s+".json"])
            
mc_samples = []
for s in mc_names:
    mc_samples.append([json_dir+s[0]+".json",s[1],s[2],s[3]])

result_dir = "./13TeV_50ns_inclusive/"
tagger = args.tagger
wp = args.workpoint


print "--- Fitting Run II Data ---"
print "  - tagger : {}".format(tagger)
print "  - workpoint : {}".format(wp)
print "--------------------------"


m = KinFtCM.Builder(tagger, wp, lumi)
for mc_s in mc_samples:
    m.add_mc_component(*mc_s)
for data_s in data_samples:
    m.add_data_component(*data_s)

# builder options
m.zeroNegativeFracs = False
m.useDataPretagNorm = False 

if (m.useDataPretagNorm):
    m.kappa_.setConstant(True)

# min data counts per pretag category
cats = m.add_all_categories(True, 10)
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






