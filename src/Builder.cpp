
#include "../interface/Builder.h"
#include "mut_framework/mut_utils/interface/prettyprint.hpp"

void Builder::add_mc_component(std::string mc_name, std::string filename,
                               double nEventGen, double theory_xs,
                               std::string k_name) {

  std::size_t mc_index = mc_names_.size();
  auto mc_iter =  std::find(mc_names_.begin(), mc_names_.end(), mc_name);
  if ( mc_iter == mc_names_.end() ) {
    mc_names_.emplace_back(mc_name);
    mc_comps_.emplace_back(filename, tagger_, workPoint_, nEventGen);
    th_xsecs_.addOwned(*new RooRealVar(("xs_"+mc_name).c_str(),
                                       ("xs_"+mc_name).c_str(),
                                       theory_xs));
    int k_index = kappas_.index(k_name.c_str());  
    if (k_index < 0) {  
       k_index = kappas_.getSize();
       kappas_.addOwned(*new RooRealVar( k_name.c_str(),
                                         k_name.c_str(),
                                         0.5, 1.5));
    }
    mc_norms_.addOwned(*new RooProduct( ("n_"+mc_name).c_str(),
                                        ("n_"+mc_name).c_str(),
                                         RooArgList(th_xsecs_[mc_index],
                                         kappas_[k_index], lumi_)));


  } else {
    mc_index = std::distance(mc_names_.begin(), mc_iter);
    std::cout << "WARNING: overwriting nominal " << mc_name << " component" << std::endl;
    mc_comps_.at(mc_index) = Component(filename, tagger_, workPoint_, nEventGen);
  }

}

void Builder::add_mc_systematic(std::string mc_name, std::string syst_name,
                                std::string filename_m, std::string filename_p,
                                double nEventGen_m, double nEventGen_p)  {
  std::size_t mc_index = mc_names_.size();
  auto mc_iter =  std::find(mc_names_.begin(), mc_names_.end(), mc_name);
  if ( mc_iter == mc_names_.end() ) {
    std::cout << "ERROR: " << mc_name << " nominal MC sample not added" << std::endl;
  } else {
    mc_index = std::distance(mc_names_.begin(), mc_iter);
    mc_systs_.at(mc_index).emplace_back();
    // get systematic index  
    std::size_t syst_index = syst_names_.size();
    auto syst_iter = std::find(syst_names_.begin(), syst_names_.end(), syst_name);
    if ( syst_iter == syst_names_.end() ) {
       // add systematic if not used yet
       syst_names_.emplace_back(syst_name);
       for (auto & mc_syst : mc_systs_) {
         mc_syst.emplace_back();
       }
    } else {
       // get index if systematic already known
       syst_index = std::distance(syst_names_.begin(), syst_iter);
    }
    auto & syst_comp = mc_systs_.at(mc_index).at(syst_index);
    syst_comp.emplace_back(filename_m, tagger_, workPoint_, nEventGen_m);
    syst_comp.emplace_back(filename_p, tagger_, workPoint_, nEventGen_p);
  }
}

void Builder::add_data_component(std::string filename) {
  data_comps_.emplace_back(filename, tagger_, workPoint_);
}


void Builder::set_mc_jet_tag_effs() {

  // get number of tagged and good each for each sample
  std::vector<std::vector<std::vector<double>>> tag_jets(n_sam_);
  std::vector<std::vector<std::vector<double>>> good_jets(n_sam_);
  for (std::size_t i_sam = 0; i_sam < n_sam_; i_sam++) {
    const auto & mc_comp = mc_comps_.at(i_sam);
    for (std::size_t i_jty = 0; i_jty < n_jty_; i_jty++) {
      const auto & type = cat_mapping_.at(i_jty);
      tag_jets.at(i_sam).emplace_back(mc_comp.get_tag_jets(type, cat_set_));
      good_jets.at(i_sam).emplace_back(mc_comp.get_good_jets(type, cat_set_));
    }
  }

  // create a SimBTagEff objects
  for (std::size_t i_jty = 0; i_jty < n_jty_; i_jty++) {
    for (std::size_t i_cat = 0; i_cat < n_cat_; i_cat++) {
      // get tagged and good jets for each cat and jty
      std::vector<double> tag_jet_jc(n_sam_, 0.0);
      std::vector<double> good_jet_jc(n_sam_, 0.0);
      for (std::size_t i_sam = 0; i_sam < n_sam_; i_sam++) {
        tag_jet_jc.at(i_sam) = tag_jets.at(i_sam).at(i_jty).at(i_cat);
        good_jet_jc.at(i_sam) = good_jets.at(i_sam).at(i_jty).at(i_cat);
      }
      std::string var_name = std::string(mc_jet_tag_effs_.at(i_jty).GetName())
                              + "_" + std::to_string(i_cat);

      mc_jet_tag_effs_.at(i_jty).addOwned(*new SimBTagEff(var_name.c_str(),
                                                          var_name.c_str(),
                                                          mc_norms_,
                                                          tag_jet_jc,
                                                          good_jet_jc));
    }
  }
}

void Builder::set_jet_tag_sfs(std::vector<bool> floating) {
  
  for (std::size_t i_jty = 0; i_jty < n_jty_; i_jty++) {
    for (std::size_t i_cat = 0; i_cat < n_cat_; i_cat++) {
       std::string var_name = std::string(jet_tag_sfs_.at(i_jty).GetName())
                               + "_" + std::to_string(i_cat);
      if (floating.at(i_jty)) {
        jet_tag_sfs_.at(i_jty).addOwned(*new RooRealVar(var_name.c_str(),
                                         var_name.c_str(), 0.0, 2.0));
      } else {
        jet_tag_sfs_.at(i_jty).addOwned(*new RooRealVar(var_name.c_str(),
                                         var_name.c_str(), 1.0));
      }
    }
  }
}

void Builder::set_jet_tag_effs() {
  
  for (std::size_t i_jty = 0; i_jty < n_jty_; i_jty++) {
    for (std::size_t i_cat = 0; i_cat < n_cat_; i_cat++) {
       std::string var_name = std::string(jet_tag_effs_.at(i_jty).GetName())
                               + "_" + std::to_string(i_cat);
       jet_tag_effs_.at(i_jty).addOwned(*new RooProduct(var_name.c_str(),
                                           var_name.c_str(), RooArgList(
                                           mc_jet_tag_effs_.at(i_jty)[i_cat],
                                           jet_tag_sfs_.at(i_jty)[i_cat])));
    }
  }
}


void Builder::add_category( const std::string & pretag_cat, const std::string & tag_cat) { 
  cat_set_.insert(std::make_pair(pretag_cat, tag_cat)); 
}


std::vector<std::string> Builder::get_mcs_names() const {
  return mc_names_;
}

void Builder::init() {

  // set values for useful varibles
  n_cat_ = mc_comps_.at(0).get_n_cat();
  n_sam_ = mc_names_.size();
  n_jty_ = cat_mapping_.size(); 
  n_sys_ = syst_names_.size();

  // setup MC b-tagging efficiencies
  set_mc_jet_tag_effs();
  set_jet_tag_sfs();  
  set_jet_tag_effs();

  // create pdfs for all categories
  for (const auto & cat : cat_set_) {
    const auto & pretag_cat = cat.first;
    const auto & tag_cat = cat.second;
    std::string cat_name = pretag_cat + ":" + tag_cat;
    kin_cat_.defineType(cat_name.c_str());
    kin_bin_pdfs_.addOwned(*get_extended_pdf_ptr(pretag_cat, tag_cat)); 
    sim_kin_pdf_.addPdf(dynamic_cast<RooAbsPdf &>(
                        kin_bin_pdfs_[kin_bin_pdfs_.getSize()-1]),
                        cat_name.c_str());
  }
}

void Builder::set_constrained_pdf() {
  list_pdfs_.add(sim_kin_pdf_);
  list_pdfs_.add(c_list_);
  fit_pdf_ = new RooProdPdf("fit_pdf","fit_pdf",list_pdfs_);
}

std::vector<std::string> Builder::add_all_categories( bool all_tag_cats, double min_counts_pretag,
                                  double min_counts_tag) {

  std::set<std::pair<std::string,std::string>> unique_cats;
  std::set<std::pair<std::string,std::string>> selected_cats;

  std::vector<std::pair<std::string, double>> pretag_pairs;

  // get all unique categories in all data files
  for (std::size_t n_s = 0; n_s < data_comps_.size(); n_s++) { 
    for ( const auto & pretag_cat : data_comps_.at(n_s).tag_cat_counts_) {
      for ( const auto & tag_cat : pretag_cat.second) {
        unique_cats.insert(std::make_pair(pretag_cat.first, tag_cat.first));
      }
    }
  }

  
  // for each unique category pair check if verify conditions
  for ( const auto & cat : unique_cats) {
    double counts_pretag = 0.0;
    double counts_tag = 0.0;
    for (std::size_t n_s = 0; n_s < data_comps_.size(); n_s++) { 
      if (data_comps_.at(n_s).tag_cat_counts_.count(cat.first) > 0) {
        const auto & pretag_map = data_comps_.at(n_s).tag_cat_counts_.at(cat.first);
        for ( const auto & tag_cat : pretag_map) {
          counts_pretag += tag_cat.second.at(0);
        }
        if (pretag_map.count(cat.second) > 0) {
          counts_tag += pretag_map.at(cat.second).at(0);
        }
      }
    }

    if ( counts_pretag > min_counts_pretag) {
      if (std::none_of(pretag_pairs.cbegin(), pretag_pairs.cend(),
          [&](const std::pair<std::string,double> & element) { return (element.first == cat.first); }))
      {
        pretag_pairs.emplace_back(cat.first, counts_pretag);
      }
      if (all_tag_cats) {
        std::vector<std::string> tag_cats = FtCM::all_tag_cats(cat.first);
        for ( const auto & tag_cat : tag_cats) {
          selected_cats.insert(std::make_pair(cat.first, tag_cat));
        }
      } else if ( counts_tag > min_counts_tag) {
        // translate tag cat name to short format
        std::string short_cat(data_comps_.at(0).get_n_cat(),'0');
        for (std::size_t b = 0; b < data_comps_.at(0).get_n_cat(); b++) {
          int jet_sum = 0;
          for (std::size_t t = 0; t < 4; t++) {
            jet_sum += int(cat.second.at(4*b+t)-'0');
          }
          short_cat.at(b) = char(jet_sum) + '0';
        }     
        selected_cats.insert(std::make_pair(cat.first, short_cat));
      }
    }
  }

  // add all selected category pair
  for ( const auto & cat : selected_cats) {
    add_category(cat.first, cat.second);
  }

  std::sort(pretag_pairs.begin(), pretag_pairs.end(), 
      [](const std::pair<std::string,double> & lhs, const std::pair<std::string,double> & rhs ) 
      {
        return lhs.second > rhs.second;
      });

  std::vector<std::string> pretag_names;
  for ( const auto & pretag_pair : pretag_pairs) {
    pretag_names.emplace_back(pretag_pair.first);
  }
  
  return pretag_names;
}

void Builder::add_pretag_category( const std::string & pretag_cat) {

  std::set<std::pair<std::string,std::string>> unique_cats;
  for (std::size_t n_s = 0; n_s < data_comps_.size(); n_s++) { 
    if ( data_comps_.at(n_s).pretag_jet_counts_.count(pretag_cat) > 0) {
      for ( const auto & tag_cat : data_comps_.at(n_s).pretag_jet_counts_.at(pretag_cat)) {
        unique_cats.insert(std::make_pair(pretag_cat, tag_cat.first));
      }
    }
  }

  for ( const auto & cat : unique_cats) {
    add_category(cat.first, cat.second);
  }
 
}

std::vector<std::string> Builder::get_tag_categories( const std::string & pretag_cat) const {

  std::vector<std::pair<std::string, double>> tag_pairs;
  for (const auto & cat : cat_set_) {
    if (cat.first == pretag_cat) {
      tag_pairs.emplace_back(cat.second, get_data_tag_counts(cat.first, cat.second));
    }
  } 

  std::sort(tag_pairs.begin(), tag_pairs.end(), 
      [](const std::pair<std::string,double> & lhs, const std::pair<std::string,double> & rhs ) 
      {
        return lhs.second > rhs.second;
      });

  std::vector<std::string> tag_names;
  for ( const auto & tag_pair : tag_pairs) {
    tag_names.emplace_back(tag_pair.first);
  }
  
  return tag_names;
}



PretagTagPdf * Builder::get_extended_pdf_ptr(const std::string & pretag_cat,
                                            const std::string & tag_cat)
{
  std::string name = "PretagTagPdf-"+pretag_cat+":"+tag_cat;
  PretagTagPdf * ext_pdf =  new PretagTagPdf(name.c_str(), name.c_str(),
                                             mc_norms_,
                                             jet_tag_effs_);   
  ext_pdf->set_category(pretag_cat, tag_cat);
  
  std::set<std::string> unique_cats;
  std::vector<std::map<std::string, std::vector<double>>> flav_fracs; 

  for (std::size_t n_s = 0; n_s < mc_comps_.size(); n_s++) { 
    flav_fracs.emplace_back(mc_comps_.at(n_s).get_flav_frac(pretag_cat));
    for ( const auto & c : flav_fracs.at(n_s)) {
      unique_cats.insert(c.first);
    }
  }
 
  std::vector<std::vector<std::string>> cat;
  std::vector<std::vector<double>> frac;

  for (const auto & unique_cat : unique_cats) {
    cat.emplace_back();
    frac.emplace_back();
    for (std::size_t n_s = 0; n_s < flav_fracs.size(); n_s++) { 
      cat.back().emplace_back(unique_cat);
      if (flav_fracs.at(n_s).count(unique_cat) > 0) {
        if (zeroNegativeFracs && (flav_fracs.at(n_s).at(unique_cat).at(0) < 0.0 )) {
          frac.back().emplace_back(0.0);
        } else {
          frac.back().emplace_back(flav_fracs.at(n_s).at(unique_cat).at(0));
        }
      } else {
        frac.back().emplace_back(0.0);
      } 
    }
  }
  
  ext_pdf->set_category(pretag_cat, tag_cat);
  ext_pdf->set_cat_frac(cat, frac);

  return ext_pdf;
  
}

double Builder::get_data_tag_counts(const std::string & pretag_cat, const std::string & tag_cat) const {

  double counts = 0.0;
  for (const auto & data_comp : data_comps_) {
    counts += data_comp.get_tag_counts(pretag_cat, tag_cat);
  }
  return counts;
}

RooDataHist Builder::get_data_hist() {
  RooDataHist data_hist("data_hist","data_hist",RooArgSet(kin_cat_));
  for (const auto & cat : cat_set_) {
    std::string cat_name = cat.first + ":" + cat.second;
    kin_cat_.setLabel(cat_name.c_str());
    data_hist.add(RooArgSet(kin_cat_),
                  get_data_tag_counts(cat.first, cat.second));

  }  
  return data_hist;
}


double Builder::get_data_pretag_counts(const std::string & pretag_cat) const {

  double counts = 0.0;
  for (const auto & data_comp : data_comps_) {
    counts += data_comp.get_pretag_counts(pretag_cat);
  }
  return counts;
}

double Builder::get_expected_pretag_counts(const std::string & pretag_cat) const {
  double counts = 0.0;
  for (std::size_t n_s = 0; n_s < mc_comps_.size(); n_s++) {
    double eff  = mc_comps_.at(n_s).get_pretag_counts(pretag_cat)/mc_comps_.at(n_s).nEventGen_;
    double factor = lumi_.getVal()*dynamic_cast<RooAbsReal&>(xsecs_[n_s]).getVal();
    counts += factor*eff;
  }
  return counts;
}
double Builder::get_expected_tag_counts(const std::string & pretag_cat, const std::string & tag_cat) const {

  std::string pdf_name = "PretagTagPdf-"+pretag_cat+":"+tag_cat;
  int pos = kin_bin_pdfs_.index(pdf_name.c_str());


  if (pos > -1) { // found in RooArgList
    PretagTagPdf & ext_pdf = dynamic_cast<PretagTagPdf &>(kin_bin_pdfs_[pos]);
    return ext_pdf.expectedEvents(RooArgSet());  
  } 

  std::cout << " Warning - Pretag category not found in RooArgList " << std::endl;

  return -1.0;
}

std::vector<double> Builder::get_mcs_pretag_counts(const std::string & pretag_cat) const {
  std::vector<double> mcs_counts;
  for (std::size_t n_s = 0; n_s < mc_comps_.size(); n_s++) {
    double eff  = mc_comps_.at(n_s).get_pretag_counts(pretag_cat)/mc_comps_.at(n_s).nEventGen_;
    double factor = lumi_.getVal()*dynamic_cast<RooAbsReal&>(xsecs_[n_s]).getVal();
    mcs_counts.emplace_back(factor*eff);
  }
  return mcs_counts;
}

std::vector<double> Builder::get_mcs_tag_counts(const std::string & pretag_cat, const std::string & tag_cat) const {

  std::string pdf_name = "PretagTagPdf-"+pretag_cat+":"+tag_cat;
  int pos = kin_bin_pdfs_.index(pdf_name.c_str());

  if (pos > -1) { // found in RooArgList
    PretagTagPdf & ext_pdf = dynamic_cast<PretagTagPdf &>(kin_bin_pdfs_[pos]);
    return ext_pdf.get_mcs_tag_counts();  
  } 

  std::cout << " Warning - Pretag category not found in RooArgList " << std::endl;

  return {};
}

