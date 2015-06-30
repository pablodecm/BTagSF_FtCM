
#pragma once

#include <vector>
#include <string>

#include <RooArgList.h>
#include <RooRealVar.h>

#include "../interface/Component.h"

class Model {
  public:

    std::vector<Component> mc_comps_;
    std::vector<Component> data_comps_;
  
    RooRealVar lumi_; 
    RooRealVar kappa_; // background normalization
    RooArgList pre_tag_effs_;
    RooArgList xsecs_;
    RooArgList tag_effs_;

    Model(double lumi) : 
       lumi_("lumi","lumi",lumi),
       kappa_("kappa","kappa",1.0,0.0,2.0),
       pre_tag_effs_("pre_tag_effs"),
       xsecs_("xsecs"),
       tag_effs_("tag_effs") {}
    ~Model() {}

    void add_mc_component(std::string filename, double nEventGen,
                          double xSec, Normalization n = FIXED); 
    void add_data_component(std::string filename); 

    void set_category_mapping( std::vector<std::vector<int>> cat_mapping) ;
    void set_tag_wp(std::string tag, double wp);

    std::vector<double> get_mc_tag_eff() const;
    std::vector<double> get_data_tag_multiplicity() const;

    
};

