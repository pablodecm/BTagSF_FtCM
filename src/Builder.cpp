
#include "../interface/Builder.h"
#include "mut_framework/mut_utils/interface/prettyprint.hpp"

void Builder::add_mc_component(std::string filename, double nEventGen,
                             double xSec, Norm n) {

  mc_comps_.emplace_back(filename, tagger_, workPoint_, nEventGen, xSec, n);
 
  // init tag_eff_lists
  if (mc_comps_.size() == 1) { 
    for (std::size_t i=0; i< mc_comps_.back().get_n_cat(); i++) {
      std::string b_tag_eff = "b_tag_eff_"+std::to_string(i);
      b_jet_tag_effs_.addOwned(*new RooRealVar(b_tag_eff.c_str(),
                                           b_tag_eff.c_str(), 0.0, 1.0));
      std::string c_tag_eff = "c_tag_eff_"+std::to_string(i);
      c_jet_tag_effs_.addOwned(*new RooRealVar(c_tag_eff.c_str(),
                                           c_tag_eff.c_str(), 0.5, 0.0, 1.0));
      dynamic_cast<RooRealVar &>(c_jet_tag_effs_[i]).setConstant();
      std::string l_tag_eff = "l_tag_eff_"+std::to_string(i);
      l_jet_tag_effs_.addOwned(*new RooRealVar(l_tag_eff.c_str(),
                                           l_tag_eff.c_str(), 0.5, 0.0, 1.0));
      dynamic_cast<RooRealVar &>(l_jet_tag_effs_[i]).setConstant();
    }
  }

  // pretag_eff and cross section
  const Component & c = mc_comps_.back();
  std::string n_pretag_eff = "pre_tag_eff_"+c.get_name();
  pretag_effs_.addOwned(*new RooRealVar(n_pretag_eff.c_str(),
                                        n_pretag_eff.c_str(),
                                        c.get_pretag_eff()[0]));
  std::string n_xsec = "xsec_"+c.get_name();
  if ( n == BKG ) { 
    xsecs_.addOwned(*new RooRealVar(n_xsec.c_str(), n_xsec.c_str(), xSec));
  } else {
    xsecs_.addOwned(*new RooRealVar(n_xsec.c_str(), n_xsec.c_str(),
                                    xSec, xSec*0.5, xSec*1.5));   
  }
  // BKG or SIGNAL
  mc_norms_.emplace_back(n);

}

void Builder::add_data_component(std::string filename) {
  data_comps_.emplace_back(filename, tagger_,workPoint_, 1.0, 1.0, DATA);
}

void Builder::add_category( const std::string & pretag_cat, const std::string & tag_cat) { 
  cat_set_.insert(std::make_pair(pretag_cat, tag_cat)); 
  std::string cat_name = pretag_cat + ":" + tag_cat;
  kin_cat_.defineType(cat_name.c_str());
  kin_bin_pdfs_.addOwned(*get_extended_pdf_ptr(pretag_cat, tag_cat)); 
  sim_kin_pdf_.addPdf(dynamic_cast<RooAbsPdf &>(kin_bin_pdfs_[kin_bin_pdfs_.getSize()-1]),
      cat_name.c_str());
}

std::vector<std::string> Builder::get_mcs_names() const {
  std::vector<std::string> mc_names;
  for (const auto & mc_comp : mc_comps_) {
    mc_names.emplace_back(mc_comp.get_name());
  }
  return mc_names;
}

std::vector<double> Builder::get_mc_jet_tag_effs(const std::vector<int> & type) const {

  std::size_t n_cat = mc_comps_.back().get_n_cat();
  std::vector<double> tag_jets(n_cat, 0.0);
  std::vector<double> good_jets(n_cat, 0.0);
  std::vector<double> jet_tag_effs(n_cat, 0.0);

  // sum all good and tagged jets
  for (std::size_t n_s = 0; n_s < mc_comps_.size(); n_s++) {
    std::vector<double> c_tag_jets = mc_comps_.at(n_s).get_tag_jets(type, cat_set_);
    std::vector<double> c_good_jets = mc_comps_.at(n_s).get_good_jets(type, cat_set_);
    for (std::size_t i_cat = 0; i_cat < n_cat; i_cat++) {
      double factor = lumi_.getVal()*dynamic_cast<RooAbsReal&>(xsecs_[n_s]).getVal();
      if (mc_norms_.at(n_s) == SIGNAL ) { 
        tag_jets.at(i_cat) += factor*c_tag_jets.at(i_cat); 
        good_jets.at(i_cat) += factor*c_good_jets.at(i_cat); 
      } else if (mc_norms_.at(n_s) == BKG) {
        tag_jets.at(i_cat) += kappa_.getVal()*factor*c_tag_jets.at(i_cat); 
        good_jets.at(i_cat) += kappa_.getVal()*factor*c_good_jets.at(i_cat); 
      }
    } 
  }

  // compute mc efficiencies
   for (std::size_t i_cat = 0; i_cat < n_cat; i_cat++) {
     std::cout << "category " << i_cat << " - " << tag_jets.at(i_cat) << " " << good_jets.at(i_cat) << std::endl;
     jet_tag_effs.at(i_cat) = tag_jets.at(i_cat) / good_jets.at(i_cat);
   }

  return jet_tag_effs;

}

void Builder::set_mc_jet_tag_effs() {
  std::vector<double> b_jet_tag_effs = get_mc_jet_tag_effs({0});
  for ( std::size_t i_c = 0; i_c < b_jet_tag_effs.size(); i_c++) {
    dynamic_cast<RooRealVar &>(b_jet_tag_effs_[i_c]).setVal(b_jet_tag_effs.at(i_c));
  }
  std::vector<double> c_jet_tag_effs = get_mc_jet_tag_effs({1});
  for ( std::size_t i_c = 0; i_c < c_jet_tag_effs.size(); i_c++) {
    dynamic_cast<RooRealVar &>(c_jet_tag_effs_[i_c]).setVal(c_jet_tag_effs.at(i_c));
  }
  std::vector<double> l_jet_tag_effs = get_mc_jet_tag_effs({2,3});
  for ( std::size_t i_c = 0; i_c < l_jet_tag_effs.size(); i_c++) {
    dynamic_cast<RooRealVar &>(l_jet_tag_effs_[i_c]).setVal(l_jet_tag_effs.at(i_c));
  }
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
  double pretag_ev_data = - 1.0;
  if (useDataPretagNorm ) {
    pretag_ev_data = get_data_pretag_counts(pretag_cat); 
  }
  PretagTagPdf * ext_pdf =  new PretagTagPdf(name.c_str(), name.c_str(),
                                           lumi_, kappa_,
                                           pretag_effs_,
                                           xsecs_,
                                           b_jet_tag_effs_,   
                                           c_jet_tag_effs_,   
                                           l_jet_tag_effs_,
                                           pretag_ev_data);   
  ext_pdf->set_norms(mc_norms_);
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
  
  ext_pdf->set_norms(mc_norms_);
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

