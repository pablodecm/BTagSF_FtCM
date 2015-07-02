
#pragma once

#include <vector>
#include <string>

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

class Model {
  public:

    std::vector<Component> mc_comps_;
    std::vector<Norm> mc_norms_;
    std::vector<Component> data_comps_;
    
    RooCategory mul_tag_; // n_tag category 

    RooRealVar lumi_; 
    RooRealVar kappa_; // background normalization
    RooArgList pretag_effs_;
    RooArgList xsecs_;
    RooArgList tag_effs_;

    RooArgList pdf_norms_;
    RooArgList ext_pdfs_;
    RooUniform uni_;
    RooSimultaneous sim_pdf_;

    Model(double lumi) : 
       mul_tag_("mul_tag", "Event Tag Multiplicity"),
       lumi_("lumi","lumi",lumi),
       kappa_("kappa","kappa",1.0,0.0,2.0),
       pretag_effs_("pretag_effs"),
       xsecs_("xsecs"),
       tag_effs_("tag_effs"),
       pdf_norms_("pdf_norms"),
       ext_pdfs_("ext_pdfs"),
       uni_("uni","uniform",RooArgSet()),
       sim_pdf_("sim_pdf","sim_pdf",mul_tag_) {}
    ~Model() {}

    void add_mc_component(std::string filename, double nEventGen,
                          double xSec, Norm n); 
    void add_data_component(std::string filename); 

    void set_category_mapping( std::vector<std::vector<int>> cat_mapping) ;
    void set_tag_wp(std::string tag, double wp);
    void set_pdfs(int max_n_tag = 4);

    std::vector<double> get_mc_tag_effs() const;
    std::vector<double> get_data_tag_multiplicity() const;

    ModelPdf get_n_tag_pdf(unsigned n_tag);
    ModelPdf * get_n_tag_pdf_ptr(unsigned n_tag);

    RooDataHist get_data_hist(int max_n_tag = 4);
    

    
};

