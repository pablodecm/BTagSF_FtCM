
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
#include "../interface/Component.h"

class ModelPdf : public RooAbsPdf {
public:
  ModelPdf() {} ; 
  ModelPdf(const char *name, const char *title,
	         RooAbsReal& lumi,
	         RooAbsReal& kappa,
	         const RooArgList& pretag_effs,
	         const RooArgList& xsecs,
	         const RooArgList& tag_effs);
  ModelPdf(const ModelPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new ModelPdf(*this,newname); }
  inline virtual ~ModelPdf() { }

  void set_cat_frac(const std::vector<std::vector<std::string>> & cat,
                    const std::vector<std::vector<double>> & frac) 
  {
    cat_ = cat;
    frac_ = frac; 
  }

  void set_n_tag(unsigned n_tag) { n_tag_ = n_tag; }
  void set_norms(std::vector<Norm> norms) { norms_ =norms;}

protected:

  unsigned n_tag_; 

  RooRealProxy lumi_ ;
  RooRealProxy kappa_ ;
  RooListProxy pretag_effs_ ;
  RooListProxy xsecs_ ;
  RooListProxy tag_effs_ ;
  std::vector<Norm> norms_;

  std::vector<std::vector<std::string>> cat_;
  std::vector<std::vector<double>> frac_;

  Double_t evaluate() const ;

private:

  ClassDef(ModelPdf,1) 
};
 
