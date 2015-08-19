#!/usr/bin/env python

from ROOT import Model, ModelPdf 
from ROOT import vector, RooRealVar, RooLognormal, RooProdPdf
from ROOT import RooFit, RooAddition, RooMinuit, RooArgList
import numpy as np
import os

json_dir = "../output/CSV_wp_scan/"

mc_samples = [[json_dir+"TTbar_Summer13.json", 6923750, 240.0, 0 ],
              [json_dir+"WJets_Summer13.json", 76041475, 36257.2, 1],
              [json_dir+"T_tW_Summer13.json", 497658, 11.15, 1],
              [json_dir+"Tbar_tW_Summer13.json", 493460, 11.15, 1],
              [json_dir+"Tbar_t_Summer13.json", 1935072, 30.0, 1],
              [json_dir+"ZJets_Summer13.json", 30459503, 3503.0, 1],
              [json_dir+"T_t_Summer13.json", 3758227, 55.531,1]]
data_samples = [[json_dir+"Data_2012ABCD_Winter13_ReReco.json"]]

lumi = 19789.0
min_n_tags = 0 
max_n_tags = 4 

#wp_list = np.linspace(0.0,1.0,20).tolist()
wp_list = [0.244, 0.679, 0.898]
tagger = "combinedSecondaryVertexBJetTags"
result_dir = "./MC_study_CVS_wp_scan_w_0_bin/"

if not os.path.exists(result_dir):
        os.makedirs(result_dir)

for  wp in wp_list:

    print "--- Fitting Run I Data ---"
    print "  - tagger : {}".format(tagger)
    print "  - workpoint : {}".format(wp)
    print "--------------------------"


    m = Model(lumi)
    for mc_s in mc_samples:
        m.add_mc_component(*mc_s)
    for data_s in data_samples:
        m.add_data_component(*data_s)

    m.set_tag_wp(tagger, wp)
    m.set_pdfs(min_n_tags, max_n_tags)

    ln_b_m0 = RooRealVar("m0","m0",1.0)
    ln_b_k = RooRealVar("k","k",1.15)
    ln_b =  RooLognormal("log_normal_b","log_normal_b",
                         m.kappa_,
                         ln_b_m0,
                         ln_b_k)
    simulpdf = RooProdPdf("simulpdf","simulpdf", m.sim_pdf_, ln_b)
 
    h_ftcm = m.get_mc_hist(min_n_tags, max_n_tags)
    nll_ftcm = simulpdf.createNLL(h_ftcm, RooFit.Extended(), RooFit.NumCPU(6))
    minuit = RooMinuit(nll_ftcm)
    fit_result = minuit.fit("r")

    fit_result.SaveAs(result_dir+"fit_result_{}_{}.root".format(tagger,wp))


