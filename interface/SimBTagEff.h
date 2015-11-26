
#pragma once

#include <vector>

#include "RooAbsReal.h" 
#include "RooListProxy.h"


class SimBTagEff : public RooAbsReal {

  public:
    SimBTagEff() {};
    SimBTagEff( const char * name, const char * title,
                const RooArgList& mc_norms,
                const std::vector<double> tag_jets,
                const std::vector<double> good_jets);
    SimBTagEff(const SimBTagEff& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const { return new SimBTagEff(*this,newname); } 
    inline virtual ~SimBTagEff() { }

  protected:

    RooListProxy mc_norms_;
    std::vector<double> tag_jets_;  
    std::vector<double> good_jets_;  

    virtual Double_t evaluate() const;

  private:

    ClassDef(SimBTagEff, 1)
};
