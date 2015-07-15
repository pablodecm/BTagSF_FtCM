
#include "Riostream.h" 
#include "../interface/ModelPdf.h" 

ClassImp(ModelPdf) 

 ModelPdf::ModelPdf(const char *name, const char *title, 
                    RooAbsReal& lumi,
                    RooAbsReal& kappa,
                    const RooArgList& pretag_effs,
                    const RooArgList& xsecs,
                    const RooArgList& tag_effs) :
   RooAbsPdf(name,title), 
   n_tag_(0),
   lumi_("lumi","lumi",this,lumi),
   kappa_("kappa","kappa",this,kappa),
   pretag_effs_("pretag_effs","pretag_effs",this),
   xsecs_("xsecs","xsecs",this),
   tag_effs_("tag_effs","tag_effs",this),
   cat_(),
   frac_()
 { 
   pretag_effs_.add(pretag_effs);
   xsecs_.add(xsecs);
   tag_effs_.add(tag_effs);
 } 


 ModelPdf::ModelPdf(const ModelPdf& other, const char* name) :  
   RooAbsPdf(other, name), 
   n_tag_(other.n_tag_),
   lumi_("lumi", this, other.lumi_),
   kappa_("kappa", this, other.kappa_),
   pretag_effs_("pretag_effs", this, other.pretag_effs_),
   xsecs_("xsecs", this, other.xsecs_),
   tag_effs_("tag_effs", this, other.tag_effs_),
   norms_(other.norms_),
   cat_(other.cat_),
   frac_(other.frac_)
 { 
 } 



 Double_t ModelPdf::evaluate() const 
 { 
   // cache some parameters (avoid access overhead)
   double value = 0.0;
   double lumi = double(lumi_);
   double kappa = double(kappa_);
   std::vector<double> tag_effs(tag_effs_.getSize(),0.0);  
   for (std::size_t j_i=0; j_i < std::size_t(tag_effs_.getSize()); j_i++) 
     tag_effs[j_i] = dynamic_cast<RooAbsReal&>(tag_effs_[j_i]).getVal();
   for (std::size_t s_i=0; s_i < cat_.size(); s_i++) { // for each sample
    double pretag_eff =  dynamic_cast<RooAbsReal&>(pretag_effs_[s_i]).getVal();
    double xsec =  dynamic_cast<RooAbsReal&>(xsecs_[s_i]).getVal();
    for (std::size_t c_i=0; c_i < cat_.at(s_i).size(); c_i++) { // for each category
      const std::string & cat = cat_.at(s_i).at(c_i);
      const double & frac = frac_.at(s_i).at(c_i);
      // for each combination
      for (const std::string & comb : FtCM::submultiset(cat, n_tag_)) { 
        double prod = 1.0;
        for (std::size_t j_i=0; j_i < cat.size(); j_i++) {
          int i = cat.at(j_i)-'0'; 
          int i_p = comb.at(j_i)-'0';
          prod *= TMath::Binomial(i, i_p)*std::pow(tag_effs[j_i],i_p)*
            std::pow(1.0-tag_effs[i], i-i_p);
        } 
        // update value sum
        if (norms_.at(s_i) == SIGNAL)  { 
          value += lumi*xsec*pretag_eff*frac*prod;
        } else if (norms_.at(s_i) == BKG) {
          value += kappa*lumi*xsec*pretag_eff*frac*prod;
        } 
      }
    } 
   }
   return value;
 } 



