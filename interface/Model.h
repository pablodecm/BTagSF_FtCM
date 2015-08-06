
#pragma once

#include <vector>
#include <string>
#include <memory>

#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooUniform.h>
#include <RooExtendPdf.h>
#include <RooSimultaneous.h>
#include <RooDataHist.h>

#include "../interface/Component.h"
#include "../interface/ModelPdf.h"
#include "../interface/PtBinPdf.h"

class Model {
  public:

    std::vector<Component> mc_comps_;
    std::vector<Norm> mc_norms_;
    std::vector<Component> data_comps_;
    
    RooCategory mul_tag_; // n_tag category 
    RooCategory kin_cat_; // kinematic category 

    RooRealVar lumi_; 
    RooRealVar kappa_; // background normalization
    RooArgList pretag_effs_;
    RooArgList xsecs_;
    RooArgList tag_effs_;
    RooArgList c_jet_tag_effs_;
    RooArgList l_jet_tag_effs_;

    std::vector<std::unique_ptr<RooArgList>> mean_b_jet_muls_;
    std::vector<std::unique_ptr<RooArgList>> mean_c_jet_muls_;
    std::vector<std::unique_ptr<RooArgList>> mean_l_jet_muls_;

    RooArgList pdf_norms_;
    RooArgList kin_pdf_norms_;
    RooArgList ext_pdfs_;
    RooArgList ext_kin_pdfs_;
    RooUniform uni_;
    RooSimultaneous sim_pdf_;
    RooSimultaneous sim_kin_pdf_;

    Model(double lumi) : 
       mul_tag_("mul_tag", "Event Tag Multiplicity"),
       kin_cat_("kin_cat", "Kinematic Category"),
       lumi_("lumi","lumi",lumi),
       kappa_("kappa","kappa",1.0,0.0,2.0),
       pretag_effs_("pretag_effs"),
       xsecs_("xsecs"),
       tag_effs_("tag_effs"),
       c_jet_tag_effs_("c_jet_tag_effs"),
       l_jet_tag_effs_("l_jet_tag_effs"),
       pdf_norms_("pdf_norms"),
       kin_pdf_norms_("kin_pdf_norms"),
       ext_pdfs_("ext_pdfs"),
       ext_kin_pdfs_("ext_kin_pdfs"),
       uni_("uni","uniform",RooArgSet()),
       sim_pdf_("sim_pdf","sim_pdf",mul_tag_),
       sim_kin_pdf_("sim_kin_pdf","sim_kin_pdf",kin_cat_)  {}
    ~Model() {}

    void add_mc_component(std::string filename, double nEventGen,
                          double xSec, Norm n); 
    void add_data_component(std::string filename); 

    void set_category_mapping( std::vector<std::vector<int>> cat_mapping) ;
    void set_tag_wp(std::string tag, double wp);
    void set_pdfs(int min_n_tag = 0, int max_n_tag = 4);

    std::vector<double> get_mc_tag_effs() const;
    std::vector<double> get_b_jet_tag_effs() const;
    std::vector<double> get_c_jet_tag_effs() const;
    std::vector<double> get_l_jet_tag_effs() const;
    std::vector<double> get_data_tag_multiplicity() const;
    std::vector<double> get_data_kin_categories() const;

    ModelPdf get_n_tag_pdf(unsigned n_tag);
    ModelPdf * get_n_tag_pdf_ptr(unsigned n_tag);

    PtBinPdf * get_pt_bin_pdf_ptr(unsigned n_bin);

    RooDataHist get_data_hist(int min_n_tag = 0, int max_n_tag = 4);
    RooDataHist get_data_kin_hist();
    

    
};

