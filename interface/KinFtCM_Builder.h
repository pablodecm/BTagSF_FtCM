
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

#include "KinFtCM_Component.h"
#include "KinFtCM_ExtendedPdf.h"

namespace KinFtCM {

  class Builder {
    public:

    std::string tagger_;  
    double workPoint_;  

    std::vector<Component> mc_comps_;
    std::vector<Norm> mc_norms_;
    std::vector<Component> data_comps_;
    
    // pretag:tag category (one for each fitted bin)
    RooCategory kin_cat_;  

    RooRealVar lumi_; 
    RooRealVar kappa_; // background normalization
    RooArgList pretag_effs_;
    RooArgList xsecs_;
    // b-tagging efficiencies
    RooArgList b_jet_tag_effs_;
    RooArgList c_jet_tag_effs_;
    RooArgList l_jet_tag_effs_;

    RooArgList kin_bin_pdfs_;
    RooSimultaneous sim_kin_pdf_;

    std::set<std::pair<std::string,std::string>> cat_set_;


    Builder(std::string tagger, double workPoint, double lumi) : 
       tagger_(tagger),
       workPoint_(workPoint),
       kin_cat_("kin_cat", "pretag:tag kinematic category"),
       lumi_("lumi","lumi",lumi),
       kappa_("kappa","kappa",1.0,0.0,2.0),
       pretag_effs_("pretag_effs"),
       xsecs_("xsecs"),
       b_jet_tag_effs_("b_jet_tag_effs"),
       c_jet_tag_effs_("c_jet_tag_effs"),
       l_jet_tag_effs_("l_jet_tag_effs"),
       kin_bin_pdfs_("ext_kin_pdfs"),
       sim_kin_pdf_("sim_kin_pdf","sim_kin_pdf",kin_cat_)  {}
    ~Builder() {}

    void add_mc_component(std::string filename, double nEventGen,
                          double xSec, Norm n); 
    void add_data_component(std::string filename); 

    void add_category( const std::string & pretag_cat, const std::string & tag_cat);
      
    std::vector<double> get_mc_jet_tag_effs( const std::vector<int> & type) const;
    void set_mc_jet_tag_effs();
    void add_all_categories(double min_counts_pretag = -1.0,
                            double min_counts_tag = -1.0 );
    void add_pretag_category(const std::string & pretag_cat);

    ExtendedPdf * get_extended_pdf_ptr(const std::string & pretag_cat, const std::string & tag_cat);
    double get_data_tag_counts(const std::string & pretag_cat, const std::string & tag_cat);
    
    RooDataHist get_data_hist();
    
  };
}

