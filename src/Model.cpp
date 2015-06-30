
#include "../interface/Model.h"

void Model::add_mc_component(std::string filename, double nEventGen,
                             double xSec, Normalization n) {
  mc_comps_.emplace_back(filename, nEventGen, xSec, n);
}

void Model::add_data_component(std::string filename) {
  data_comps_.emplace_back(filename, 1.0, 1.0, DATA);
}

void Model::set_category_mapping( std::vector<std::vector<int>> cat_mapping) {
  for (auto & mc_comp : mc_comps_) mc_comp.set_category_mapping(cat_mapping);
  for (auto & data_comp : data_comps_) data_comp.set_category_mapping(cat_mapping);
}

void Model::set_tag_wp(std::string tag, double wp) {
  for (auto & mc_comp : mc_comps_) mc_comp.set_tag_wp(tag, wp);
  for (auto & data_comp : data_comps_) data_comp.set_tag_wp(tag, wp);
}

std::vector<double> Model::get_mc_tag_eff() const {

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



