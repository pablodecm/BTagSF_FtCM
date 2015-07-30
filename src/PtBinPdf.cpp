
#include "Riostream.h" 
#include "../interface/PtBinPdf.h" 

ClassImp(PtBinPdf) 

 PtBinPdf::PtBinPdf(const char *name, const char *title,
                    unsigned n_bin,
	                  RooAbsReal& lumi,
	                  RooAbsReal& kappa,
	                  const RooArgList& pretag_effs,
	                  const RooArgList& xsecs,
	                  const RooArgList& tag_effs,
	                  const RooArgList& mean_b_jet_muls,
	                  const RooArgList& mean_c_jet_muls,
	                  const RooArgList& mean_l_jet_muls,
                    RooAbsReal& c_jet_tag_eff,
                    RooAbsReal& l_jet_tag_eff) :
   RooAbsPdf(name,title), 
   n_bin_(n_bin),
   lumi_("lumi","lumi",this,lumi),
   kappa_("kappa","kappa",this,kappa),
   pretag_effs_("pretag_effs","pretag_effs",this),
   xsecs_("xsecs","xsecs",this),
   tag_effs_("tag_effs","tag_effs",this),
   mean_b_jet_muls_("mean_b_jet_muls","mean_b_jet_muls",this),
   mean_c_jet_muls_("mean_c_jet_muls","mean_c_jet_muls",this),
   mean_l_jet_muls_("mean_l_jet_muls","mean_l_jet_muls",this),
   c_jet_tag_eff_("c_jet_tag_eff","c_jet_tag_eff",this,c_jet_tag_eff),
   l_jet_tag_eff_("l_jet_tag_eff","l_jet_tag_eff",this,l_jet_tag_eff)
 { 
   pretag_effs_.add(pretag_effs);
   xsecs_.add(xsecs);
   tag_effs_.add(tag_effs);
   mean_b_jet_muls_.add(mean_b_jet_muls);
   mean_c_jet_muls_.add(mean_c_jet_muls);
   mean_l_jet_muls_.add(mean_l_jet_muls);
 } 


 PtBinPdf::PtBinPdf(const PtBinPdf& other, const char* name) :  
   RooAbsPdf(other, name), 
   n_bin_(other.n_bin_),
   lumi_("lumi", this, other.lumi_),
   kappa_("kappa", this, other.kappa_),
   pretag_effs_("pretag_effs", this, other.pretag_effs_),
   xsecs_("xsecs", this, other.xsecs_),
   tag_effs_("tag_effs", this, other.tag_effs_),
   norms_(other.norms_),
   mean_b_jet_muls_("mean_b_jet_muls",this, other.mean_b_jet_muls_),
   mean_c_jet_muls_("mean_c_jet_muls",this, other.mean_c_jet_muls_),
   mean_l_jet_muls_("mean_l_jet_muls",this, other.mean_l_jet_muls_),
   c_jet_tag_eff_("c_jet_tag_eff",this,other.c_jet_tag_eff_),
   l_jet_tag_eff_("l_jet_tag_eff",this,other.l_jet_tag_eff_)
 { 
 } 



 Double_t PtBinPdf::evaluate() const 
 { 
   // cache some parameters (avoid access overhead)
   double value = 0.0;
   double lumi = double(lumi_);
   double kappa = double(kappa_);
   double b_jet_tag_eff = dynamic_cast<RooAbsReal&>(tag_effs_[n_bin_+2]).getVal();
   double c_jet_tag_eff = double(c_jet_tag_eff_);
   double l_jet_tag_eff = double(l_jet_tag_eff_);
   for (std::size_t s_i=0; s_i < std::size_t(xsecs_.getSize()); s_i++) {
    double pretag_eff =  dynamic_cast<RooAbsReal&>(pretag_effs_[s_i]).getVal();
    double xsec =  dynamic_cast<RooAbsReal&>(xsecs_[s_i]).getVal();
    double mean_b_jet_mul =  dynamic_cast<RooAbsReal&>(mean_b_jet_muls_[s_i]).getVal();
    double mean_c_jet_mul =  dynamic_cast<RooAbsReal&>(mean_c_jet_muls_[s_i]).getVal();
    double mean_l_jet_mul =  dynamic_cast<RooAbsReal&>(mean_l_jet_muls_[s_i]).getVal();
      if (norms_.at(s_i) == SIGNAL)  { 
        value += lumi*xsec*pretag_eff*( mean_b_jet_mul*b_jet_tag_eff +
                                        mean_c_jet_mul*c_jet_tag_eff +
                                        mean_l_jet_mul*l_jet_tag_eff );
      } else if (norms_.at(s_i) == BKG) {
        value += kappa*lumi*xsec*pretag_eff*( mean_b_jet_mul*b_jet_tag_eff +
                                              mean_c_jet_mul*c_jet_tag_eff +
                                              mean_l_jet_mul*l_jet_tag_eff );
      } 
   } 
   return value;
 } 


