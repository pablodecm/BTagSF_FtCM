
#pragma once

#include <vector>
#include <cmath> 

#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "TMath.h" 

#include "../interface/FtCM_formula.h"
#include "../interface/KinFtCM_Component.h"

namespace KinFtCM {
class ExtendedPdf : public RooAbsPdf {
public:
  ExtendedPdf() {} ; 
  ExtendedPdf(const char *name, const char *title,
	         RooAbsReal& lumi,
	         RooAbsReal& kappa,
	         const RooArgList& pretag_effs,
	         const RooArgList& xsecs,
	         const RooArgList& b_jet_tag_effs,
	         const RooArgList& c_jet_tag_effs,
	         const RooArgList& l_jet_tag_effs);
  ExtendedPdf(const ExtendedPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new ExtendedPdf(*this,newname); }
  inline virtual ~ExtendedPdf() { }

  void set_cat_frac(const std::vector<std::vector<std::string>> & cat,
                    const std::vector<std::vector<double>> & frac) 
  {
    cat_ = cat;
    frac_ = frac; 
  }

  void set_norms(std::vector<Norm> norms) { norms_ =norms;}
  void set_category(const std::string & pretag_cat, const std::string & tag_cat)
  { pretag_cat_ = pretag_cat;
    tag_cat_ = tag_cat_; }
  virtual ExtendMode extendMode() const { return CanBeExtended ; }
  virtual Double_t expectedEvents(const RooArgSet* nset) const;
  virtual Double_t expectedEvents(const RooArgSet& nset) const { return expectedEvents(&nset); };

protected:

  RooRealProxy lumi_ ;
  RooRealProxy kappa_ ;
  RooListProxy pretag_effs_ ;
  RooListProxy xsecs_ ;
  RooListProxy tag_effs_ ;
  RooListProxy b_jet_tag_effs_;
  RooListProxy c_jet_tag_effs_;
  RooListProxy l_jet_tag_effs_;

  std::string pretag_cat_;
  std::string tag_cat_;

  std::vector<Norm> norms_;

  std::vector<std::vector<std::string>> cat_;
  std::vector<std::vector<double>> frac_;

  Double_t evaluate() const ;

private:

  ClassDef(ExtendedPdf,1) 
};

}
 
