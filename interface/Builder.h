
#pragma once

#include <vector>
#include <string>
#include <memory>
#include <iterator>

#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooUniform.h>
#include <RooProduct.h>
#include <RooExtendPdf.h>
#include <RooSimultaneous.h>
#include <RooDataHist.h>

#include "Component.h"
#include "PretagTagPdf.h"
#include "SimBTagEff.h"


class Builder {
  public:

  std::string tagger_;  
  double workPoint_;  

  std::vector<std::string> jet_types_ = {"b", "c", "l" };

  // nominal MC components 
  std::vector<Component> mc_comps_;
  // systematic shifted components
  // i (different samples)
  // j (different syst) 
  // k ( 0 is down & 1 is up)
  std::vector<std::vector<std::vector<Component>>> mc_systs_;
  // data components
  std::vector<Component> data_comps_;

  // names of samples and systematics
  std::vector<std::string> mc_names_; 
  std::vector<std::string> syst_names_; 
  
  // pretag:tag category (one for each fitted bin)
  RooCategory kin_cat_;  

  // luminosity
  RooRealVar nom_lumi_; 
  RooRealVar unc_lumi_;
  RooProduct lumi_; 

  // cross section variables
  RooArgList th_xsecs_;
  RooArgList kappas_;
  RooArgList xsecs_;
  // MC normalization lumi*xsecs
  RooArgList mc_norms_;

  // b-tagging efficiency (0 is b, 1 is c, 2 is l) 
  // each element of the list correspong to a kin. cat.
  std::vector<RooArgList> mc_jet_tag_effs_;  
  std::vector<RooArgList> jet_tag_sfs_;
  std::vector<RooArgList> jet_tag_effs_;

  // shape syst nuissance parameters
  RooArgList shape_nuis_pars_;

  RooArgList kin_bin_pdfs_;
  RooSimultaneous sim_kin_pdf_;

  // number of categories, samples and jet types
  std::size_t n_cat_ = 0;
  std::size_t n_sam_ = 0;
  std::size_t n_jty_ = 0;
  std::size_t n_sys_ = 0;

  std::vector<std::vector<int>> cat_mapping_ = {{0}, {1}, {2,3}};

  std::set<std::pair<std::string,std::string>> cat_set_;

  // configuration options
  bool zeroNegativeFracs = false;

  Builder(std::string tagger, double workPoint, double lumi) : 
     tagger_(tagger),
     workPoint_(workPoint),
     kin_cat_("kin_cat", "pretag:tag kinematic category"),
     nom_lumi_("nom_lumi", "Nominal luminosity", lumi),
     unc_lumi_("unc_lumi", "Luminosity uncertainty nuissance", 1.0),
     lumi_("lumi","Luminosity", RooArgList(nom_lumi_, unc_lumi_)),
     th_xsecs_("th_xsecs"),
     kappas_("kappas"),
     xsecs_("xsecs"),
     mc_norms_("mc_norms"),
     kin_bin_pdfs_("ext_kin_pdfs"),
     sim_kin_pdf_("sim_kin_pdf","sim_kin_pdf",kin_cat_)  {

    // construct RooArgList 
    for (const auto & jet_type : jet_types_) {
      mc_jet_tag_effs_.emplace_back((jet_type+"_mc_jet_tag_effs").c_str());
      jet_tag_sfs_.emplace_back((jet_type+"_jet_tag_sfs").c_str());
      jet_tag_effs_.emplace_back((jet_type+"_jet_tag_effs").c_str());
    }

    // set names again (were deleted, probably a ROOT bug)
    for (std::size_t i_jty = 0; i_jty < jet_types_.size(); i_jty++) {
      const auto & jet_type = jet_types_.at(i_jty);
      mc_jet_tag_effs_.at(i_jty).setName((jet_type+"_mc_jet_tag_effs").c_str());
      jet_tag_sfs_.at(i_jty).setName((jet_type+"_mc_jet_tag_sfs").c_str());
      jet_tag_effs_.at(i_jty).setName((jet_type+"_jet_tag_effs").c_str());
    }
  }

  ~Builder() { }

  void add_mc_component(std::string mc_name, std::string filename,
                        double nEventGen, double theory_xs, std::string k_name); 
  void add_mc_systematic(std::string mc_name, std::string syst_name,
                         std::string filename_m, std::string filename_p,
                         double nEventGen_m, double nEventGen_p);
  void add_data_component(std::string filename); 

  void set_mc_jet_tag_effs();
  void set_jet_tag_sfs(std::vector<bool> floating = {true, false, false});
  void set_jet_tag_effs();

  void add_category( const std::string & pretag_cat, const std::string & tag_cat);

  std::vector<std::string> get_mcs_names() const;
    
  std::vector<double> get_mc_jet_tag_effs( const std::vector<int> & type) const;
  void init();
  std::vector<std::string> add_all_categories(bool all_tag_cats = true, double min_counts_pretag = -1.0,
                                              double min_counts_tag = -1.0 );
  void add_pretag_category(const std::string & pretag_cat);
  std::vector<std::string> get_tag_categories(const std::string & pretag_cat) const;

  PretagTagPdf * get_extended_pdf_ptr(const std::string & pretag_cat, const std::string & tag_cat);
  double get_data_tag_counts(const std::string & pretag_cat, const std::string & tag_cat) const;
  double get_data_pretag_counts(const std::string & pretag_cat) const;
  double get_expected_pretag_counts(const std::string & pretag_cat) const;
  double get_expected_tag_counts(const std::string & pretag_cat, const std::string & tag_cat) const;
  std::vector<double> get_mcs_pretag_counts(const std::string & pretag_cat) const;
  std::vector<double> get_mcs_tag_counts(const std::string & pretag_cat, const std::string & tag_cat) const;

  std::size_t get_n_cat() const {return n_cat_; }
  std::size_t get_n_sam() const {return n_sam_; }
  std::size_t get_n_jty() const {return n_jty_; }

  RooDataHist get_data_hist();
  
};

