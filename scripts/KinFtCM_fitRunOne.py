#!/usr/bin/env python

from ROOT import KinFtCM, RooArgSet 
from ROOT import vector, RooRealVar, RooLognormal, RooProdPdf
from ROOT import RooFit, RooAddition, RooMinuit, RooArgList
import numpy as np
import os

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

wp_list = [0.679]
tagger = "combinedSecondaryVertexBJetTags"
result_dir = "./KinFtCM_tests/"

if not os.path.exists(result_dir):
        os.makedirs(result_dir)

for  wp in wp_list:

    print "--- Fitting Run I Data ---"
    print "  - tagger : {}".format(tagger)
    print "  - workpoint : {}".format(wp)
    print "--------------------------"


    m = KinFtCM.Builder(tagger, wp, lumi)
    for mc_s in mc_samples:
        m.add_mc_component(*mc_s)
    for data_s in data_samples:
        m.add_data_component(*data_s)

    #m.add_category("00112", "00000")
    #m.add_category("00112", "00100")
    m.add_all_categories()

    ext_pdf_0 = m.get_extended_pdf_ptr("00112", "00000")
    ext_pdf_1 = m.get_extended_pdf_ptr("00112", "00100")
    m.get_extended_pdf_ptr("00130", "00000")
    #ext_pdf_1 = m.get_extended_pdf_ptr("00220", "00100")
    m.set_mc_jet_tag_effs()
    print ext_pdf_0.expectedEvents(RooArgSet())
    print ext_pdf_1.expectedEvents(RooArgSet())



