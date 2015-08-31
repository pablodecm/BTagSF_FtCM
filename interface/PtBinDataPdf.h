
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

class PtBinDataPdf : public RooAbsPdf {
public:
  PtBinDataPdf() {} ; 
  PtBinDataPdf(const char *name, const char *title,
           unsigned n_bin,
           double n_data_ev,
           RooAbsReal& lumi,
	         RooAbsReal& kappa,
	         const RooArgList& pretag_effs,
	         const RooArgList& xsecs,
	         const RooArgList& tag_effs,
	         const RooArgList& mean_b_jet_muls,
	         const RooArgList& mean_c_jet_muls,
	         const RooArgList& mean_l_jet_muls,
           RooAbsReal& c_jet_tag_eff,
           RooAbsReal& l_jet_tag_eff );
  PtBinDataPdf(const PtBinDataPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new PtBinDataPdf(*this,newname); }
  inline virtual ~PtBinDataPdf() { }

  void set_norms(std::vector<Norm> norms) { norms_ =norms;}

protected:

  unsigned n_bin_; 
  double n_data_ev_;

  RooRealProxy lumi_ ;
  RooRealProxy kappa_ ;
  RooListProxy pretag_effs_ ;
  RooListProxy xsecs_ ;
  RooListProxy tag_effs_ ;
  std::vector<Norm> norms_;
  RooListProxy mean_b_jet_muls_ ;
  RooListProxy mean_c_jet_muls_ ;
  RooListProxy mean_l_jet_muls_ ;
  RooRealProxy c_jet_tag_eff_ ;
  RooRealProxy l_jet_tag_eff_ ;

  Double_t evaluate() const ;

private:

  ClassDef(PtBinDataPdf,1) 
};
 
