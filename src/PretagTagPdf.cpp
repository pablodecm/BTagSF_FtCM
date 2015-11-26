
#include "Riostream.h" 
#include "../interface/PretagTagPdf.h" 

ClassImp(PretagTagPdf) 


 PretagTagPdf::PretagTagPdf(const char *name, const char *title, 
                            const RooArgList& mc_norms,
                            const std::vector<RooArgList> & jet_tag_effs) :
   RooAbsPdf(name,title), 
   mc_norms_("mc_norms","mc_norms", this),
   jet_tag_effs_()
 { 
   mc_norms_.add(mc_norms);
   for (std::size_t c=0; c < jet_tag_effs.size(); c++) {
     jet_tag_effs_.push_back(std::unique_ptr<RooListProxy>(
                             new RooListProxy(
                               jet_tag_effs.at(c).GetName(),
                               jet_tag_effs.at(c).GetTitle(),
                               this)));
     jet_tag_effs_.at(c)->add(jet_tag_effs.at(c));
   }
  
 } 


 PretagTagPdf::PretagTagPdf(const PretagTagPdf & other, const char* name) :  
   RooAbsPdf(other, name), 
   mc_norms_("mc_norms", this, other.mc_norms_),
   jet_tag_effs_(),
   pretag_cat_(other.pretag_cat_),
   tag_cat_(other.tag_cat_),
   cat_(other.cat_),
   frac_(other.frac_),
   syst_frac_(other.syst_frac_)
 { 
   for (const auto & proxy : other.jet_tag_effs_) {
     jet_tag_effs_.push_back(std::unique_ptr<RooListProxy>(
                             new RooListProxy(
                               proxy->GetName(),
                               this,
                               *proxy)));
   }
 } 



 Double_t PretagTagPdf::evaluate() const 
 { 
   return 1.0;
 } 

Double_t PretagTagPdf::expectedEvents(const RooArgSet* nset) const {

  std::size_t n_jty = jet_tag_effs_.size();
  std::size_t n_cat = jet_tag_effs_.at(0)->getSize();
  std::size_t n_sam = mc_norms_.getSize();
  double value = 0.0;

  std::vector<std::vector<double>> jet_tag_effs( n_cat, std::vector<double>(n_jty, 0.0));  
  for (std::size_t j_i=0; j_i < n_cat; j_i++) { // for each jet kin bin
    for (std::size_t t_i=0; t_i < n_jty; t_i++) { // for each jet type
    jet_tag_effs.at(j_i).at(t_i) = dynamic_cast<RooAbsReal&>((*jet_tag_effs_.at(t_i))[j_i]).getVal();
    }
   }

  std::vector<double> mc_norms(n_sam,0.0);  
  for (std::size_t s_i=0; s_i < n_sam ; s_i++) { // for each sample 
     mc_norms.at(s_i) = dynamic_cast<RooAbsReal&>(mc_norms_[s_i]).getVal();
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
          for ( std::size_t j_t = 0; j_t < 3; j_t++) { // for each jet
            int i = pre_cats.at(j_i).at(j_t)-'0'; 
            int i_p = pos_cat.at(j_i).at(j_t)-'0';
            pos_prob *= TMath::Binomial(i, i_p)*std::pow(jet_tag_effs.at(j_i).at(j_t),i_p)*
                       std::pow(1.0-jet_tag_effs.at(j_i).at(j_t), i-i_p); 
          }
        }
        value += mc_norms.at(s_i)*frac*pos_prob;
      }  
    }
  }

  // avoid negative expected events (negative MC weights)
  if (value < 0.0) {
    value = 0.0;
  }
  return value;
}

std::vector<double> PretagTagPdf::get_mcs_tag_counts() const {

  std::size_t n_jty = jet_tag_effs_.size();
  std::size_t n_cat = jet_tag_effs_.at(0)->getSize();
  std::size_t n_sam = mc_norms_.getSize();
  std::vector<double> mcs_tag_counts(n_sam, 0.0);


  std::vector<std::vector<double>> jet_tag_effs( n_cat, std::vector<double>(n_jty, 0.0));  
  for (std::size_t j_i=0; j_i < n_cat; j_i++) { // for each jet kin bin
    for (std::size_t t_i=0; t_i < n_jty; t_i++) { // for each jet type
    jet_tag_effs.at(j_i).at(t_i) = dynamic_cast<RooAbsReal&>((*jet_tag_effs_.at(t_i))[j_i]).getVal();
    }
   }

  std::vector<double> mc_norms(n_sam,0.0);  
  for (std::size_t s_i=0; s_i < n_sam ; s_i++) { // for each sample 
     mc_norms.at(s_i) = dynamic_cast<RooAbsReal&>(mc_norms_[s_i]).getVal();
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
  
        mcs_tag_counts.at(s_i) += mc_norms.at(s_i)*frac*pos_prob;
      }  
    }
  }

  for (auto & mc_tag_counts : mcs_tag_counts) {
    // avoid negative expected events (negative MC weights)
    if (mc_tag_counts < 0.0) {
      mc_tag_counts = 0.0;
    }
  } 

  return mcs_tag_counts;
}

