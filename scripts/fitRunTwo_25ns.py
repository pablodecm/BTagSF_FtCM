#!/usr/bin/env python

from mut_framework.LJMET_converter.samples_25ns import data_samples,mc_samples
from mut_framework.BTagSF_FtCM.ScaleFactorReader import ScaleFactorReader 
import argparse
from ROOT import Builder, RooArgSet, RooGaussian 
from ROOT import vector, RooRealVar, RooLognormal, RooProdPdf
from ROOT import RooFit, RooAddition, RooMinuit, RooArgList
from ROOT import SetOwnership
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


mc_names = [["TT_powheg","ttbar_norm"],
            ["WJetsToLNu","bkg_norm"],
            ["DYJetsToLL_M-50", "bkg_norm"],
            ["ST_t-channel", "bkg_norm"],
            ["ST_s-channel", "bkg_norm"],
            ["ST_tW_top","bkg_norm"],
            ["ST_tW_antitop", "bkg_norm"],
            ["QCD_Pt-1000toInf_MuEnrichedPt5","bkg_norm"],
            ["QCD_Pt-120to170_MuEnrichedPt5","bkg_norm"],
            ["QCD_Pt-170to300_MuEnrichedPt5","bkg_norm"],
            ["QCD_Pt-300to470_MuEnrichedPt5","bkg_norm"],
            ["QCD_Pt-470to600_MuEnrichedPt5","bkg_norm"],
            ["QCD_Pt-50to80_MuEnrichedPt5","bkg_norm"],
            ["QCD_Pt-600to800_MuEnrichedPt5", "bkg_norm"],
            ["QCD_Pt-800to1000_MuEnrichedPt5","bkg_norm"],
            ["QCD_Pt-80to120_MuEnrichedPt5", "bkg_norm"],
            ["WW", "bkg_norm"],
            ["WZ", "bkg_norm"],
            ["ZZ", "bkg_norm"]]




data_comps = []
lumi = 0.0
for s in data_names:
    data_comps.append([json_dir+data_samples[s]["full_name"]+".json"])
    lumi += data_samples[s]["lumi"]
            
mc_comps = []
for s in mc_names:
    mc_comps.append([s[0],json_dir+mc_samples[s[0]]["full_name"]+".json",
                     mc_samples[s[0]]["gen_events"], mc_samples[s[0]]["xs"],
                     s[1]])

result_dir =  args.output_folder
tagger = args.tagger
wp = args.workpoint

if tagger == "pfCombinedInclusiveSecondaryVertexV2BJetTags": 
    s_tag = "CSVv2"
    op = {0.604 : 0, 0.890: 1, 0.970 : 2}[wp]


print "--- Fitting Run II 25ns Data ---"
print "  - tagger : {}".format(tagger)
print "  - workpoint : {}".format(wp)
print "    {} - OP {}".format(s_tag, op)
print "--------------------------"


m = Builder(tagger, wp, lumi)
for mc_s in mc_comps:
    m.add_mc_component(*mc_s)
for data_s in data_comps:
    m.add_data_component(*data_s)

# builder options
m.zeroNegativeFracs = False

# min data counts per pretag category
cats = m.add_all_categories(True, 200)
m.init()

n_cat = m.get_n_cat()
bin_centers = [[0.0, 200]]
csv_file = "../data/ext_SF/RunIISpring15MiniAODv2/CSVv2.csv"
jet_tag_sfs_constraints = [["free", [1.0]*n_cat ], # b jet sf is floating free
                           ["around", m.jet_tag_sfs_[0], [0.1]*n_cat],
                           ["meas", "comb", bin_centers]]    
                          # ["fixed", [1.0]*n_cat]]    

for t, const in enumerate(jet_tag_sfs_constraints):
    if const[0] == "free":
        for cat in range(n_cat):
            m.jet_tag_sfs_[t][cat].setVal(const[1][cat])
            m.jet_tag_sfs_[t][cat].setConstant(False)
    elif const[0] == "fixed":        
            m.jet_tag_sfs_[t][cat].setVal(const[1][cat])
            m.jet_tag_sfs_[t][cat].setConstant(True)
    elif const[0] == "around":        
        for cat in range(n_cat):
            x_a_r = m.jet_tag_sfs_[t][cat]
            mean_a_r = const[1][cat] 
            sigma_a = RooRealVar(x_a_r.GetName()+"_sigma",
                                 x_a_r.GetName()+"_sigma",
                                 const[2][cat])
            m.aux_list_.addOwned(sigma_a)
            sigma_a_i = m.aux_list_.getSize()-1
            SetOwnership(sigma_a, False)
            sigma_a_r = m.aux_list_[sigma_a_i]                              
            c_name = x_a_r.GetName()+"_around_"+mean_a_r.GetName()
            gauss_a = RooGaussian(c_name, c_name, x_a_r, mean_a_r, sigma_a_r)
            m.c_list_.addOwned(gauss_a)
            SetOwnership(gauss_a, False)
    elif const[0] == "meas":      
        for cat in range(n_cat):
            x_m_r = m.jet_tag_sfs_[t][cat]
            r = ScaleFactorReader(s_tag,csv_file, op, const[1]) 
            eta = bin_centers[cat][0]
            pt = bin_centers[cat][1]
            mean_m = RooRealVar(x_m_r.GetName()+"_mean",
                                x_m_r.GetName()+"_mean",
                                r.get_sf("central",t, eta, pt))
            m.aux_list_.addOwned(mean_m)
            mean_m_i = m.aux_list_.getSize()-1
            SetOwnership(mean_m, False)
            sigma_m = RooRealVar(x_m_r.GetName()+"_sigma",
                                 x_m_r.GetName()+"_sigma",
                                 r.get_sf("up",t, eta, pt)/2. -
                                 r.get_sf("down",t,eta,pt)/2.)
            m.aux_list_.addOwned(sigma_m)
            sigma_m_i = m.aux_list_.getSize()-1
            SetOwnership(sigma_m, False)
            mean_m_r = m.aux_list_[mean_m_i]
            sigma_m_r = m.aux_list_[sigma_m_i]
            c_name = x_m_r.GetName()+"_measured"
            gauss_m = RooGaussian(c_name, c_name, x_m_r,
                                  mean_m_r, sigma_m_r)
            m.c_list_.addOwned(gauss_m)
            SetOwnership(gauss_m, False)

m.set_constrained_pdf()
# get and fit data
data = m.get_data_hist()
fit_result = m.fit_pdf_.fitTo(data, RooFit.Extended(1),RooFit.Constrained(),
                                    RooFit.NumCPU(8), RooFit.Save(1), RooFit.Minos(True))


