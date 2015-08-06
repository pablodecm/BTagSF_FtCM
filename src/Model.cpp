
#include "../interface/Model.h"

void Model::add_mc_component(std::string filename, double nEventGen,
                             double xSec, Norm n) {

  mc_comps_.emplace_back(filename, nEventGen, xSec, n);
 
  // init tag_eff_list
  if (mc_comps_.size() == 1) { 
    for (std::size_t i=0; i< mc_comps_.back().get_n_cat(); i++) {
      std::string n_tag_eff = "tag_eff_"+std::to_string(i);
      tag_effs_.addOwned(*new RooRealVar(n_tag_eff.c_str(),
                                         n_tag_eff.c_str(), 0.0, 1.0));
      if (i==0 || i ==1) { // light and c jets eff from mc
        dynamic_cast<RooRealVar &>(tag_effs_[i]).setConstant();
      }
    }
  }

  // pretag_eff and cross section
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
  // BKG or SIGNAL
  mc_norms_.emplace_back(n);

  // save mean multiplicites for all flavours  to arg lists
  std::vector<double> mean_b_jet_mul = c.get_mean_b_jet_mul();
  for (std::size_t i_c=0; i_c < mean_b_jet_mul.size() ; i_c++) {
    std::string n_mean_b_jet_mul= "mean_b_jet_mul_bin"+std::to_string(i_c);
    if (mc_comps_.size() == 1) {
      mean_b_jet_muls_.emplace_back(new RooArgList(n_mean_b_jet_mul.c_str()));
    }
    mean_b_jet_muls_.at(i_c)->addOwned( 
        *new RooRealVar( (n_mean_b_jet_mul+"_"+ c.get_name()).c_str(),
                         (n_mean_b_jet_mul+"_"+ c.get_name()).c_str(),
                         mean_b_jet_mul.at(i_c)));
  }


  std::vector<double> mean_c_jet_mul = c.get_mean_c_jet_mul();
  for (std::size_t i_c=0; i_c < mean_c_jet_mul.size() ; i_c++) {
    std::string n_mean_c_jet_mul= "mean_c_jet_mul_bin"+std::to_string(i_c);
    if (mc_comps_.size() == 1) {
      mean_c_jet_muls_.emplace_back(new RooArgList(n_mean_c_jet_mul.c_str()));
      std::string n_c_jet_tag_eff = "c_jet_tag_eff_bin"+std::to_string(i_c);
      c_jet_tag_effs_.addOwned(*new RooRealVar(n_c_jet_tag_eff.c_str(),
                                               n_c_jet_tag_eff.c_str(), 
                                               0.0, 1.0));
      dynamic_cast<RooRealVar &>(c_jet_tag_effs_[i_c]).setConstant();

    }
    mean_c_jet_muls_.at(i_c)->addOwned( 
        *new RooRealVar( (n_mean_c_jet_mul+"_"+ c.get_name()).c_str(),
                         (n_mean_c_jet_mul+"_"+ c.get_name()).c_str(),
                         mean_c_jet_mul.at(i_c)));
  }

  std::vector<double> mean_l_jet_mul = c.get_mean_l_jet_mul();
  for (std::size_t i_c=0; i_c < mean_l_jet_mul.size() ; i_c++) {
    std::string n_mean_l_jet_mul= "mean_l_jet_mul_bin"+std::to_string(i_c);
    if (mc_comps_.size() == 1) {
      mean_l_jet_muls_.emplace_back(new RooArgList(n_mean_l_jet_mul.c_str()));
      std::string n_l_jet_tag_eff = "l_jet_tag_eff_bin"+std::to_string(i_c);
      l_jet_tag_effs_.addOwned(*new RooRealVar(n_l_jet_tag_eff.c_str(),
                                               n_l_jet_tag_eff.c_str(), 
                                               0.0, 1.0));
      dynamic_cast<RooRealVar &>(l_jet_tag_effs_[i_c]).setConstant();
    }
    mean_l_jet_muls_.at(i_c)->addOwned( 
        *new RooRealVar( (n_mean_l_jet_mul+"_"+ c.get_name()).c_str(),
                         (n_mean_l_jet_mul+"_"+ c.get_name()).c_str(),
                         mean_l_jet_mul.at(i_c)));
  }

  std::cout << "Component added" << std::endl;


   

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
  std::vector<double> c_jet_tag_effs = get_c_jet_tag_effs();
  for ( std::size_t i_c = 0; i_c < c_jet_tag_effs.size(); i_c++) {
    dynamic_cast<RooRealVar &>(c_jet_tag_effs_[i_c]).setVal(c_jet_tag_effs.at(i_c));
  }
  std::vector<double> l_jet_tag_effs = get_l_jet_tag_effs();
  for ( std::size_t i_c = 0; i_c < l_jet_tag_effs.size(); i_c++) {
    dynamic_cast<RooRealVar &>(l_jet_tag_effs_[i_c]).setVal(l_jet_tag_effs.at(i_c));
  }
}

void Model::set_pdfs(int min_n_tag, int max_n_tag) {

  // set tag multiplicity pdf
  for (std::size_t n_t=0; n_t<=std::size_t(max_n_tag); n_t++) {
    pdf_norms_.addOwned(*get_n_tag_pdf_ptr(n_t));
    std::string n_ext_pdf = "ext_" + std::string(pdf_norms_[n_t].GetName());
    RooExtendPdf * ext_pdf = new RooExtendPdf(n_ext_pdf.c_str(), n_ext_pdf.c_str(),
                                 uni_, dynamic_cast<RooAbsPdf &>(pdf_norms_[n_t]));
    ext_pdfs_.addOwned(*ext_pdf);
    if (n_t >= std::size_t(min_n_tag)) {
      std::string n_tag_cat  = std::to_string(n_t) + "_tag_jets";
      mul_tag_.defineType(n_tag_cat.c_str());
      sim_pdf_.addPdf(dynamic_cast<RooAbsPdf &>(ext_pdfs_[n_t]), n_tag_cat.c_str());
    }
  }

  // set kinematic bin pdf
  for (std::size_t n_b=0; n_b<mean_b_jet_muls_.size(); n_b++) {
    kin_pdf_norms_.addOwned(*get_pt_bin_pdf_ptr(n_b));
    std::string n_ext_pdf = "ext_" + std::string(kin_pdf_norms_[n_b].GetName());
    RooExtendPdf * ext_pdf = new RooExtendPdf(n_ext_pdf.c_str(), n_ext_pdf.c_str(),
                                 uni_, dynamic_cast<RooAbsPdf &>(kin_pdf_norms_[n_b]));
    ext_kin_pdfs_.addOwned(*ext_pdf);
    std::string kin_bin_cat  = std::to_string(n_b) + "_kin_bin";
    std::cout << kin_bin_cat << std::endl;
    kin_cat_.defineType(kin_bin_cat.c_str());
    sim_kin_pdf_.addPdf(dynamic_cast<RooAbsPdf &>(ext_kin_pdfs_[n_b]), kin_bin_cat.c_str());
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

std::vector<double> Model::get_b_jet_tag_effs() const {

  // init vectors to zero
  std::size_t n_cat =  mc_comps_.at(0).get_good_b_jets().size();
  std::vector<double> good_b_jets(n_cat, 0.0);
  std::vector<double> tag_b_jets(n_cat, 0.0);
  std::vector<double> b_jet_tag_eff(n_cat, 0.0);

  // sum all good and tagged jets
  for (const auto & mc_comp : mc_comps_) {
   std::vector<double> c_good_b_jets = mc_comp.get_good_b_jets(); 
   std::vector<double> c_tag_b_jets = mc_comp.get_tag_b_jets(); 
   for (std::size_t i_cat = 0; i_cat < n_cat; i_cat++) {
     good_b_jets.at(i_cat) += c_good_b_jets.at(i_cat); 
     tag_b_jets.at(i_cat) += c_tag_b_jets.at(i_cat); 
   } 
  }

  // compute mc efficiencies
   for (std::size_t i_cat = 0; i_cat < n_cat; i_cat++) {
     b_jet_tag_eff.at(i_cat) = tag_b_jets.at(i_cat) / good_b_jets.at(i_cat);
   }
  
  return b_jet_tag_eff; 
}

std::vector<double> Model::get_c_jet_tag_effs() const {

  // init vectors to zero
  std::size_t n_cat =  mc_comps_.at(0).get_good_c_jets().size();
  std::vector<double> good_c_jets(n_cat, 0.0);
  std::vector<double> tag_c_jets(n_cat, 0.0);
  std::vector<double> c_jet_tag_eff(n_cat, 0.0);

  // sum all good and tagged jets
  for (const auto & mc_comp : mc_comps_) {
   std::vector<double> c_good_c_jets = mc_comp.get_good_c_jets(); 
   std::vector<double> c_tag_c_jets = mc_comp.get_tag_c_jets(); 
   for (std::size_t i_cat = 0; i_cat < n_cat; i_cat++) {
     good_c_jets.at(i_cat) += c_good_c_jets.at(i_cat); 
     tag_c_jets.at(i_cat) += c_tag_c_jets.at(i_cat); 
   } 
  }

  // compute mc efficiencies
   for (std::size_t i_cat = 0; i_cat < n_cat; i_cat++) {
     c_jet_tag_eff.at(i_cat) = tag_c_jets.at(i_cat) / good_c_jets.at(i_cat);
   }
  
  return c_jet_tag_eff; 
}

std::vector<double> Model::get_l_jet_tag_effs() const {

  // init vectors to zero
  std::size_t n_cat =  mc_comps_.at(0).get_good_l_jets().size();
  std::vector<double> good_l_jets(n_cat, 0.0);
  std::vector<double> tag_l_jets(n_cat, 0.0);
  std::vector<double> l_jet_tag_eff(n_cat, 0.0);

  // sum all good and tagged jets
  for (const auto & mc_comp : mc_comps_) {
   std::vector<double> c_good_l_jets = mc_comp.get_good_l_jets(); 
   std::vector<double> c_tag_l_jets = mc_comp.get_tag_l_jets(); 
   for (std::size_t i_cat = 0; i_cat < n_cat; i_cat++) {
     good_l_jets.at(i_cat) += c_good_l_jets.at(i_cat); 
     tag_l_jets.at(i_cat) += c_tag_l_jets.at(i_cat); 
   } 
  }

  // compute mc efficiencies
   for (std::size_t i_cat = 0; i_cat < n_cat; i_cat++) {
     l_jet_tag_eff.at(i_cat) = tag_l_jets.at(i_cat) / good_l_jets.at(i_cat);
   }
  
  return l_jet_tag_eff; 
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

std::vector<double> Model::get_data_kin_categories() const {

  std::vector<double> kin_cat_counts = data_comps_.at(0).get_tag_x_jets();
  
  for (std::size_t n_s = 1; n_s < data_comps_.size(); n_s++) {
    std::vector<double> kin_cat_counts_add = data_comps_.at(n_s).get_tag_x_jets();
    for (std::size_t i_j = 0; i_j < kin_cat_counts.size(); i_j++) {
     kin_cat_counts.at(i_j) += kin_cat_counts_add.at(i_j);
   }
  }

  return kin_cat_counts;
}

std::vector<double> Model::get_mc_tag_multiplicity() const {

  std::vector<double> tag_multiplicity;
  
  double lumi = lumi_.getVal();
  double kappa = kappa_.getVal();

  for (std::size_t n_s = 0; n_s < mc_comps_.size(); n_s++) {
    double nEventGen =  mc_comps_.at(n_s).nEventGen_;
    double xsec =  dynamic_cast<RooAbsReal&>(xsecs_[n_s]).getVal();
    std::vector<double> tag_mul_mc = mc_comps_.at(n_s).get_tag_multiplicity();
    if ( n_s == 0 ) tag_multiplicity = std::vector<double> (tag_mul_mc.size(),0.0);
    for (std::size_t i_j = 0; i_j < tag_multiplicity.size(); i_j++) {
     if (mc_norms_.at(n_s) == SIGNAL ) { 
      tag_multiplicity.at(i_j) += tag_mul_mc.at(i_j)*lumi*xsec/nEventGen;
     } else if (mc_norms_.at(n_s) == BKG) {
      tag_multiplicity.at(i_j) += tag_mul_mc.at(i_j)*kappa*lumi*xsec/nEventGen;
     }
   }
  }

  return tag_multiplicity;
}

std::vector<double> Model::get_mc_kin_categories() const {

  std::vector<double> kin_cat_counts;
  
  double lumi = lumi_.getVal();
  double kappa = kappa_.getVal();

  for (std::size_t n_s = 0; n_s < mc_comps_.size(); n_s++) {
    double nEventGen =  mc_comps_.at(n_s).nEventGen_;
    double xsec =  dynamic_cast<RooAbsReal&>(xsecs_[n_s]).getVal();
    std::vector<double> kin_cat_mc = mc_comps_.at(n_s).get_good_jets();
    if ( n_s == 0 ) kin_cat_counts = std::vector<double> (kin_cat_mc.size(),0.0);
    for (std::size_t i_j = 0; i_j < kin_cat_counts.size(); i_j++) {
     if (mc_norms_.at(n_s) == SIGNAL ) { 
      kin_cat_counts.at(i_j) += kin_cat_mc.at(i_j)*lumi*xsec/nEventGen;
     } else if (mc_norms_.at(n_s) == BKG) {
      kin_cat_counts.at(i_j) += kin_cat_mc.at(i_j)*kappa*lumi*xsec/nEventGen;
     }
   }
  }

    return kin_cat_counts;
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

PtBinPdf * Model::get_pt_bin_pdf_ptr(unsigned n_bin) {

  // create pdf pointer
  std::string pdf_name = "pt_bin_" + std::to_string(n_bin) + "_pdf";
  PtBinPdf * p_pdf = new PtBinPdf(pdf_name.c_str(), pdf_name.c_str(),
                                  n_bin, lumi_, kappa_, pretag_effs_, 
                                  xsecs_, tag_effs_,
                                  *mean_b_jet_muls_.at(n_bin),
                                  *mean_c_jet_muls_.at(n_bin),
                                  *mean_l_jet_muls_.at(n_bin),
                                  dynamic_cast<RooAbsReal &>(c_jet_tag_effs_[n_bin]),
                                  dynamic_cast<RooAbsReal &>(l_jet_tag_effs_[n_bin]));
  p_pdf->set_norms(mc_norms_);
  return p_pdf;
}


RooDataHist Model::get_data_hist(int min_n_tag, int max_n_tag) {
  std::vector<double> data_tag_mul = get_data_tag_multiplicity();
  RooDataHist data_hist("data","data", RooArgSet(mul_tag_));
  for (std::size_t n_t=std::size_t(min_n_tag); n_t<=std::size_t(max_n_tag); n_t++) {
    std::string n_tag_cat  = std::to_string(n_t) + "_tag_jets";
    mul_tag_.setLabel(n_tag_cat.c_str());
    data_hist.add(RooArgSet(mul_tag_), data_tag_mul.at(n_t));
  }
  return data_hist;
}

RooDataHist Model::get_mc_hist(int min_n_tag, int max_n_tag) {
  std::vector<double> mc_tag_mul = get_mc_tag_multiplicity();
  RooDataHist mc_hist("mc_hist","mc_hist", RooArgSet(mul_tag_));
  for (std::size_t n_t=std::size_t(min_n_tag); n_t<=std::size_t(max_n_tag); n_t++) {
    std::string n_tag_cat  = std::to_string(n_t) + "_tag_jets";
    mul_tag_.setLabel(n_tag_cat.c_str());
    mc_hist.add(RooArgSet(mul_tag_), mc_tag_mul.at(n_t));
  }
  return mc_hist;
}

RooDataHist Model::get_data_kin_hist() {
  std::vector<double> data_kin_cat = get_data_kin_categories();
  RooDataHist data_kin_hist("data_kin_hist","data_kin_hist", RooArgSet(kin_cat_));
  for (std::size_t n_c=0; n_c<data_kin_cat.size(); n_c++) {
    kin_cat_.setIndex(n_c);
    data_kin_hist.add(RooArgSet(kin_cat_), data_kin_cat.at(n_c));
  }
  return data_kin_hist;
}

RooDataHist Model::get_mc_kin_hist() {
  std::vector<double> mc_kin_cat = get_mc_kin_categories();
  RooDataHist mc_kin_hist("mc_kin_hist","mc_kin_hist", RooArgSet(kin_cat_));
  for (std::size_t n_c=0; n_c<mc_kin_cat.size(); n_c++) {
    kin_cat_.setIndex(n_c);
    mc_kin_hist.add(RooArgSet(kin_cat_), mc_kin_cat.at(n_c));
  }
  return mc_kin_hist;
}


