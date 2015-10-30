
#include "Riostream.h" 
#include "../interface/KinFtCM_ExtendedPdf.h" 
#include "mut_framework/mut_utils/interface/prettyprint.hpp"

ClassImp(KinFtCM::ExtendedPdf) 

namespace KinFtCM {

 ExtendedPdf::ExtendedPdf(const char *name, const char *title, 
                    RooAbsReal& lumi,
                    RooAbsReal& kappa,
                    const RooArgList& pretag_effs,
                    const RooArgList& xsecs,
                    const RooArgList& b_jet_tag_effs,
                    const RooArgList& c_jet_tag_effs,
                    const RooArgList& l_jet_tag_effs,
                    double pretag_ev_data) :
   RooAbsPdf(name,title), 
   lumi_("lumi","lumi",this,lumi),
   kappa_("kappa","kappa",this,kappa),
   pretag_effs_("pretag_effs","pretag_effs",this),
   xsecs_("xsecs","xsecs",this),
   b_jet_tag_effs_("tag_effs","tag_effs",this),
   c_jet_tag_effs_("tag_effs","tag_effs",this),
   l_jet_tag_effs_("tag_effs","tag_effs",this),
   cat_(),
   frac_(),
   pretag_ev_data_(pretag_ev_data),
   pretag_cat_effs_(xsecs.getSize(), 0.0)
 { 
   pretag_effs_.add(pretag_effs);
   xsecs_.add(xsecs);
   b_jet_tag_effs_.add(b_jet_tag_effs);
   c_jet_tag_effs_.add(c_jet_tag_effs);
   l_jet_tag_effs_.add(l_jet_tag_effs);
  
 } 


 ExtendedPdf::ExtendedPdf(const ExtendedPdf & other, const char* name) :  
   RooAbsPdf(other, name), 
   lumi_("lumi", this, other.lumi_),
   kappa_("kappa", this, other.kappa_),
   pretag_effs_("pretag_effs", this, other.pretag_effs_),
   xsecs_("xsecs", this, other.xsecs_),
   b_jet_tag_effs_("b_jet_tag_effs", this, other.b_jet_tag_effs_),
   c_jet_tag_effs_("c_jet_tag_effs", this, other.c_jet_tag_effs_),
   l_jet_tag_effs_("l_jet_tag_effs", this, other.l_jet_tag_effs_),
   pretag_cat_(other.pretag_cat_),
   tag_cat_(other.tag_cat_),
   norms_(other.norms_),
   cat_(other.cat_),
   frac_(other.frac_),
   pretag_ev_data_(other.pretag_ev_data_),
   pretag_cat_effs_(other.pretag_cat_effs_)
 { 
 } 



 Double_t ExtendedPdf::evaluate() const 
 { 
   return 1.0;
 } 

Double_t ExtendedPdf::expectedEvents(const RooArgSet* nset) const {

  std::size_t n_cat = b_jet_tag_effs_.getSize();
  std::size_t n_sam = xsecs_.getSize();
  double value = 0.0;
  double lumi = double(lumi_);
  double kappa = double(kappa_);
  double w_sum = 0.0;

  std::vector<std::vector<double>> jet_tag_effs( n_cat, std::vector<double>(3, 0.0));  
  std::vector<double> c_jet_tag_effs(c_jet_tag_effs_.getSize(),0.0);  
  std::vector<double> l_jet_tag_effs(l_jet_tag_effs_.getSize(),0.0);  
  for (std::size_t j_i=0; j_i < n_cat; j_i++) { // for each jet kin bin
    jet_tag_effs.at(j_i).at(0) = dynamic_cast<RooAbsReal&>(b_jet_tag_effs_[j_i]).getVal();
    jet_tag_effs.at(j_i).at(1) = dynamic_cast<RooAbsReal&>(c_jet_tag_effs_[j_i]).getVal();
    jet_tag_effs.at(j_i).at(2) = dynamic_cast<RooAbsReal&>(l_jet_tag_effs_[j_i]).getVal();
   }

  std::vector<double> xsecs(xsecs_.getSize(),0.0);  
  std::vector<double> pretag_effs(pretag_effs_.getSize(),0.0);  
  for (std::size_t s_i=0; s_i < n_sam ; s_i++) { // for each sample 
     xsecs.at(s_i) = dynamic_cast<RooAbsReal&>(xsecs_[s_i]).getVal();
     pretag_effs.at(s_i) = dynamic_cast<RooAbsReal&>(pretag_effs_[s_i]).getVal();
     if (norms_.at(s_i) == SIGNAL)  { 
        w_sum += lumi*xsecs.at(s_i)*pretag_effs.at(s_i)*pretag_cat_effs_.at(s_i);
     } else if (norms_.at(s_i) == BKG) {
        w_sum += kappa*lumi*xsecs.at(s_i)*pretag_effs.at(s_i)*pretag_cat_effs_.at(s_i);
     } 
  }

  for (std::size_t c_i=0; c_i < cat_.size(); c_i++) { // for each kinematic category
    for (std::size_t s_i=0; s_i < cat_.at(c_i).size(); s_i++) { // for each sample
      const std::string & cat = cat_.at(c_i).at(s_i);
      const double & frac = frac_.at(c_i).at(s_i);

      std::vector<std::string> pre_cats;
      std::vector<std::vector<std::string>> sub_cats;
      for (std::size_t j_i=0; j_i < n_cat; j_i++) { // for each jet kin bin
        pre_cats.emplace_back(cat.substr(3*j_i, 3));
        sub_cats.emplace_back();  
        for (const auto & sub_cat : FtCM::submultiset(pre_cats.back(), tag_cat_.at(j_i) - '0')) {
          sub_cats.back().emplace_back(sub_cat);
        }
      }

      // cartesian product to get all possiblities
      std::vector<std::vector<std::string>> pos_cats = FtCM::cart_product(sub_cats);

      for ( const auto & pos_cat : pos_cats) { // for each posibility
        double pos_prob = 1;
        for ( std::size_t j_i = 0; j_i < pos_cat.size(); j_i++) { // for each jet kin 
          for ( std::size_t j_t = 0; j_t < 3; j_t++) {
            int i = pre_cats.at(j_i).at(j_t)-'0'; 
            int i_p = pos_cat.at(j_i).at(j_t)-'0';
            pos_prob *= TMath::Binomial(i, i_p)*std::pow(jet_tag_effs.at(j_i).at(j_t),i_p)*
                       std::pow(1.0-jet_tag_effs.at(j_i).at(j_t), i-i_p); 
          }
        }
        // update value sum
        if (norms_.at(s_i) == SIGNAL)  { 
          value += lumi*xsecs.at(s_i)*pretag_effs.at(s_i)*frac*pos_prob;
        } else if (norms_.at(s_i) == BKG) {
          value += kappa*lumi*xsecs.at(s_i)*pretag_effs.at(s_i)*frac*pos_prob;
        } 
      }  
    }
  }

  // use pretag data counts as normalization
  if (pretag_ev_data_ > 0) {
    if (w_sum > 1e-50) { // avoid null division
      value *= pretag_ev_data_/w_sum;
    } else {
      value = 1e-50;
    }
  }
  // avoid negative expected events (negative MC weights)
  if (value < 0.0) {
    value = 0.0;
  }
  return value;
}

std::vector<double> ExtendedPdf::get_mcs_tag_counts() const {

  std::size_t n_cat = b_jet_tag_effs_.getSize();
  std::size_t n_sam = xsecs_.getSize();
  std::vector<double> mcs_tag_counts(n_sam, 0.0);
  double lumi = double(lumi_);
  double kappa = double(kappa_);
  double w_sum = 0.0;

  std::vector<std::vector<double>> jet_tag_effs( n_cat, std::vector<double>(3, 0.0));  
  std::vector<double> c_jet_tag_effs(c_jet_tag_effs_.getSize(),0.0);  
  std::vector<double> l_jet_tag_effs(l_jet_tag_effs_.getSize(),0.0);  
  for (std::size_t j_i=0; j_i < n_cat; j_i++) { // for each jet kin bin
    jet_tag_effs.at(j_i).at(0) = dynamic_cast<RooAbsReal&>(b_jet_tag_effs_[j_i]).getVal();
    jet_tag_effs.at(j_i).at(1) = dynamic_cast<RooAbsReal&>(c_jet_tag_effs_[j_i]).getVal();
    jet_tag_effs.at(j_i).at(2) = dynamic_cast<RooAbsReal&>(l_jet_tag_effs_[j_i]).getVal();
   }

  std::vector<double> xsecs(xsecs_.getSize(),0.0);  
  std::vector<double> pretag_effs(pretag_effs_.getSize(),0.0);  
  for (std::size_t s_i=0; s_i < n_sam ; s_i++) { // for each sample 
     xsecs.at(s_i) = dynamic_cast<RooAbsReal&>(xsecs_[s_i]).getVal();
     pretag_effs.at(s_i) = dynamic_cast<RooAbsReal&>(pretag_effs_[s_i]).getVal();
     if (norms_.at(s_i) == SIGNAL)  { 
        w_sum += lumi*xsecs.at(s_i)*pretag_effs.at(s_i)*pretag_cat_effs_.at(s_i);

     } else if (norms_.at(s_i) == BKG) {
        w_sum += kappa*lumi*xsecs.at(s_i)*pretag_effs.at(s_i)*pretag_cat_effs_.at(s_i);
     } 
  }


  for (std::size_t c_i=0; c_i < cat_.size(); c_i++) { // for each kinematic category
    for (std::size_t s_i=0; s_i < cat_.at(c_i).size(); s_i++) { // for each sample
      const std::string & cat = cat_.at(c_i).at(s_i);
      const double & frac = frac_.at(c_i).at(s_i);

      std::vector<std::string> pre_cats;
      std::vector<std::vector<std::string>> sub_cats;
      for (std::size_t j_i=0; j_i < n_cat; j_i++) { // for each jet kin bin
        pre_cats.emplace_back(cat.substr(3*j_i, 3));
        sub_cats.emplace_back();  
        for (const auto & sub_cat : FtCM::submultiset(pre_cats.back(), tag_cat_.at(j_i) - '0')) {
          sub_cats.back().emplace_back(sub_cat);
        }
      }

      // cartesian product to get all possiblities
      std::vector<std::vector<std::string>> pos_cats = FtCM::cart_product(sub_cats);

      for ( const auto & pos_cat : pos_cats) { // for each posibility
        double pos_prob = 1;
        for ( std::size_t j_i = 0; j_i < pos_cat.size(); j_i++) { // for each jet kin 
          for ( std::size_t j_t = 0; j_t < 3; j_t++) {
            int i = pre_cats.at(j_i).at(j_t)-'0'; 
            int i_p = pos_cat.at(j_i).at(j_t)-'0';
            pos_prob *= TMath::Binomial(i, i_p)*std::pow(jet_tag_effs.at(j_i).at(j_t),i_p)*
                       std::pow(1.0-jet_tag_effs.at(j_i).at(j_t), i-i_p); 
          }
        }
  
        // update value sum
        if (norms_.at(s_i) == SIGNAL)  { 
          mcs_tag_counts.at(s_i) += lumi*xsecs.at(s_i)*pretag_effs.at(s_i)*frac*pos_prob;
        } else if (norms_.at(s_i) == BKG) {
          mcs_tag_counts.at(s_i) += kappa*lumi*xsecs.at(s_i)*pretag_effs.at(s_i)*frac*pos_prob;
        } 
      }  
    }
  }

  for (auto & mc_tag_counts : mcs_tag_counts) {
    // use pretag data counts as normalization
    if (pretag_ev_data_ > 0) {
      if (w_sum > 1e-5) { // avoid null division
        mc_tag_counts *= pretag_ev_data_/w_sum;
      } else {
        mc_tag_counts = 0.0;
      }
    }
    // avoid negative expected events (negative MC weights)
    if (mc_tag_counts < 0.0) {
      mc_tag_counts = 0.0;
    }
  } 

  return mcs_tag_counts;
}


}


