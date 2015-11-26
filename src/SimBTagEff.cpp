#include "Riostream.h" 
#include "../interface/SimBTagEff.h" 

ClassImp(SimBTagEff)

SimBTagEff::SimBTagEff(const char * name, const char * title,
                       const RooArgList& mc_norms,
                       const std::vector<double> tag_jets,
                       const std::vector<double> good_jets) :
  RooAbsReal(name, title),
  mc_norms_("mc_norms","mc_norms", this),
  tag_jets_(tag_jets),
  good_jets_(good_jets)
{
  mc_norms_.add(mc_norms);
}

SimBTagEff::SimBTagEff(const SimBTagEff& other, const char* name) :
  RooAbsReal(other, name),
  mc_norms_("mc_norms", this, other.mc_norms_),
  tag_jets_(other.tag_jets_),
  good_jets_(other.good_jets_) 
{
}


  
  

double SimBTagEff::evaluate() const {

   std::size_t n_sam = mc_norms_.getSize();
   double tag_jet = 0.0;
   double good_jet = 0.0;
   for (std::size_t s_i=0; s_i < n_sam ; s_i++) { // for each sample 
      double norm_s = dynamic_cast<RooAbsReal&>(mc_norms_[s_i]).getVal();
      tag_jet += norm_s*tag_jets_.at(s_i);
      good_jet += norm_s*good_jets_.at(s_i);
   }

   return (tag_jet/good_jet);
}


