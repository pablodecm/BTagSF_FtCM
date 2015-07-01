
#pragma once

#include <vector>
#include <string>

#include <RooArgList.h>
#include <RooRealVar.h>

#include "../interface/Component.h"
#include "../interface/ModelPdf.h"

class Model {
  public:

    std::vector<Component> mc_comps_;
    std::vector<Norm> mc_norms_;
    std::vector<Component> data_comps_;
    
  
    RooRealVar lumi_; 
    RooRealVar kappa_; // background normalization
    RooArgList pretag_effs_;
    RooArgList xsecs_;
    RooArgList tag_effs_;

    Model(double lumi) : 
       lumi_("lumi","lumi",lumi),
       kappa_("kappa","kappa",1.0,0.0,2.0),
       pretag_effs_("pretag_effs"),
       xsecs_("xsecs"),
       tag_effs_("tag_effs") {}
    ~Model() {}

    void add_mc_component(std::string filename, double nEventGen,
                          double xSec, Norm n); 
    void add_data_component(std::string filename); 

    void set_category_mapping( std::vector<std::vector<int>> cat_mapping) ;
    void set_tag_wp(std::string tag, double wp);

    std::vector<double> get_mc_tag_effs() const;
    std::vector<double> get_data_tag_multiplicity() const;

    ModelPdf get_n_tag_pdf(unsigned n_tag);

    
};

