
#pragma once

#include <vector>
#include <cmath> 

#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsCategory.h"
#include "TMath.h" 

#include "CombinatorialFormulas.h"
#include "Component.h"

#include "mut_framework/mut_utils/interface/prettyprint.hpp"

class PretagTagPdf : public RooAbsPdf {
public:
  PretagTagPdf() {} ; 
  PretagTagPdf(const char *name, const char *title,
	             const RooArgList& mc_norms,
	             const std::vector<RooArgList> & jet_tag_effs);
  PretagTagPdf(const PretagTagPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new PretagTagPdf(*this,newname); }
  inline virtual ~PretagTagPdf() { }

  void set_cat_frac(const std::vector<std::vector<std::string>> & cat,
                    const std::vector<std::vector<double>> & frac,
                    const std::vector<std::vector<std::vector<std::vector<double>>>> & syst_frac = {}) 
  {
    cat_ = cat;
    frac_ = frac; 
    syst_frac_ = syst_frac;
  }

  void set_category(const std::string & pretag_cat, const std::string & tag_cat)
  { pretag_cat_ = pretag_cat;
    tag_cat_ = tag_cat;
  }

  virtual ExtendMode extendMode() const { return CanBeExtended ; }
  virtual Double_t expectedEvents(const RooArgSet* nset) const;
  virtual Double_t expectedEvents(const RooArgSet& nset) const { return expectedEvents(&nset); };
  std::vector<double> get_mcs_tag_counts() const;

protected:

  RooListProxy mc_norms_ ;
  std::vector<std::unique_ptr<RooListProxy>> jet_tag_effs_;

  std::string pretag_cat_;
  std::string tag_cat_;

  std::vector<std::vector<std::string>> cat_;
  std::vector<std::vector<double>> frac_;
  std::vector<std::vector<std::vector<std::vector<double>>>> syst_frac_; 

  Double_t evaluate() const ;

private:

  ClassDef(PretagTagPdf,1) 
};

 
