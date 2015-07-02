
#include "../interface/Model.h"

void Model::add_mc_component(std::string filename, double nEventGen,
                             double xSec, Norm n) {
  mc_comps_.emplace_back(filename, nEventGen, xSec, n);
  if (mc_comps_.size() == 1) { //init tag_eff_list
    for (std::size_t i=0; i< mc_comps_.back().get_n_cat(); i++) {
      std::string n_tag_eff = "tag_eff_"+std::to_string(i);
      tag_effs_.addOwned(*new RooRealVar(n_tag_eff.c_str(),
                                         n_tag_eff.c_str(), 0.0, 1.0));
      if (i==0 || i ==1) { // light and c jets eff from mc
        dynamic_cast<RooRealVar &>(tag_effs_[i]).setConstant();
      }
    }
  }
  const Component & c = mc_comps_.back();
  std::string n_pretag_eff = "pre_tag_eff_"+c.get_name();
  pretag_effs_.addOwned(*new RooRealVar(n_pretag_eff.c_str(),
                                        n_pretag_eff.c_str(),
                                        c.get_pretag_eff()[0]));
  std::string n_xsec = "xsec_"+c.get_name();
  if ( n == BKG ) { 
    xsecs_.addOwned(*new RooRealVar(n_xsec.c_str(), n_xsec.c_str(), xSec));
  } else {
    xsecs_.addOwned(*new RooRealVar(n_xsec.c_str(), n_xsec.c_str(),
                                    xSec, xSec*0.5, xSec*1.5));   
  }
  mc_norms_.emplace_back(n);
                                     
}

void Model::add_data_component(std::string filename) {
  data_comps_.emplace_back(filename, 1.0, 1.0, DATA);
}

void Model::set_category_mapping( std::vector<std::vector<int>> cat_mapping) {
  tag_effs_.removeAll();
  for (auto & mc_comp : mc_comps_) mc_comp.set_category_mapping(cat_mapping);
  for (auto & data_comp : data_comps_) data_comp.set_category_mapping(cat_mapping);
  for (std::size_t i=0; i< cat_mapping.size(); i++) {
    std::string n_tag_eff = "tag_eff_"+std::to_string(i);
    tag_effs_.addOwned(*new RooRealVar(n_tag_eff.c_str(),
                                       n_tag_eff.c_str(), 0.0, 1.0));
  }
}

void Model::set_tag_wp(std::string tag, double wp) {
  for (auto & mc_comp : mc_comps_) mc_comp.set_tag_wp(tag, wp);
  for (auto & data_comp : data_comps_) data_comp.set_tag_wp(tag, wp);

  std::vector<double> mc_effs = get_mc_tag_effs();
  for ( std::size_t i_c = 0; i_c < mc_effs.size(); i_c++) {
    dynamic_cast<RooRealVar &>(tag_effs_[i_c]).setVal(mc_effs[i_c]);
  } 

}

void Model::set_pdfs(int max_n_tag) {

  for (std::size_t n_t=0; n_t<std::size_t(max_n_tag); n_t++) {
    pdf_norms_.addOwned(*get_n_tag_pdf_ptr(n_t));
    std::string n_ext_pdf = "ext" + std::string(pdf_norms_[n_t].GetName());
    RooExtendPdf * ext_pdf = new RooExtendPdf(n_ext_pdf.c_str(), n_ext_pdf.c_str(),
                                 uni_, dynamic_cast<RooAbsPdf &>(pdf_norms_[n_t]));
    ext_pdfs_.addOwned(*ext_pdf);
    std::string n_tag_cat  = std::to_string(n_t) + "_tag_jets";
    mul_tag_.defineType(n_tag_cat.c_str());
    sim_pdf_.addPdf(dynamic_cast<RooAbsPdf &>(ext_pdfs_[n_t]), n_tag_cat.c_str());
  }
}

std::vector<double> Model::get_mc_tag_effs() const {

  // init vectors to zero
  std::size_t n_cat =  mc_comps_.at(0).get_n_cat();
  std::vector<double> good_cat_jets(n_cat, 0.0);
  std::vector<double> tag_cat_jets(n_cat, 0.0);
  std::vector<double> mc_tag_eff(n_cat, 0.0);

  // sum all good and tagged jets
  for (const auto & mc_comp : mc_comps_) {
   std::vector<double> c_good_cat_jets = mc_comp.get_good_cat_jets(); 
   std::vector<double> c_tag_cat_jets = mc_comp.get_tag_cat_jets(); 
   for (std::size_t i_cat = 0; i_cat < n_cat; i_cat++) {
     good_cat_jets.at(i_cat) += c_good_cat_jets.at(i_cat); 
     tag_cat_jets.at(i_cat) += c_tag_cat_jets.at(i_cat); 
   } 
  }

  // compute mc efficiencies
   for (std::size_t i_cat = 0; i_cat < n_cat; i_cat++) {
     mc_tag_eff.at(i_cat) = tag_cat_jets.at(i_cat) / good_cat_jets.at(i_cat);
   }
  
  return mc_tag_eff; 
}

std::vector<double> Model::get_data_tag_multiplicity() const {

  std::vector<double> tag_multiplicity = data_comps_.at(0).get_tag_multiplicity();
  
  for (std::size_t n_s = 1; n_s < data_comps_.size(); n_s++) {
    std::vector<double> c_tag_mul = data_comps_.at(n_s).get_tag_multiplicity();
    for (std::size_t i_j = 0; i_j < tag_multiplicity.size(); i_j++) {
     tag_multiplicity.at(i_j) += c_tag_mul.at(i_j);
   }
  }

  return tag_multiplicity;
}

ModelPdf Model::get_n_tag_pdf(unsigned n_tag) {

  // create pdf
  std::string pdf_name = std::to_string(n_tag) + "tag_pdf";
  ModelPdf m_pdf(pdf_name.c_str(), pdf_name.c_str(),
                 lumi_, kappa_,
                 pretag_effs_, xsecs_, tag_effs_);
  m_pdf.set_n_tag(n_tag);
  m_pdf.set_norms(mc_norms_);

  std::vector<std::vector<std::string>> cat;
  std::vector<std::vector<double>> frac;

  for (const auto & mc_comp : mc_comps_) {
    cat.emplace_back();
    frac.emplace_back();
    for (const auto & pair : mc_comp.get_cat_fractions()) {
      cat.back().emplace_back(pair.first);
      frac.back().emplace_back(pair.second[0]);
    }
  }

  m_pdf.set_cat_frac(cat, frac);

  return m_pdf;
}

ModelPdf * Model::get_n_tag_pdf_ptr(unsigned n_tag) {
  
  // create pdf pointer
  std::string pdf_name = std::to_string(n_tag) + "tag_pdf";
  ModelPdf * m_pdf = new ModelPdf(pdf_name.c_str(), pdf_name.c_str(),
                                  lumi_, kappa_,
                                  pretag_effs_, xsecs_, tag_effs_);
  m_pdf->set_n_tag(n_tag);
  m_pdf->set_norms(mc_norms_);

  std::vector<std::vector<std::string>> cat;
  std::vector<std::vector<double>> frac;

  for (const auto & mc_comp : mc_comps_) {
    cat.emplace_back();
    frac.emplace_back();
    for (const auto & pair : mc_comp.get_cat_fractions()) {
      cat.back().emplace_back(pair.first);
      frac.back().emplace_back(pair.second[0]);
    }
  }

  m_pdf->set_cat_frac(cat, frac);

  return m_pdf;
}

RooDataHist Model::get_data_hist(int max_n_tag) {
  std::vector<double> data_tag_mul = get_data_tag_multiplicity();
  RooDataHist data_hist("data","data", RooArgSet(mul_tag_));
  for (std::size_t n_t=0; n_t<std::size_t(max_n_tag); n_t++) {
    mul_tag_.setIndex(n_t);
    data_hist.add(RooArgSet(mul_tag_), data_tag_mul[n_t]);
  }
  return data_hist;
}




