
#include "../interface/Model.h"

typedef std::string JetCategory;

Component::Component(std::string filename, double nEventGen,
                     double xSec, Norm norm = FIXED) :
                     c_name_(filename.substr(0,filename.find_last_of("."))),
                     nEventGen_(nEventGen),
                     xSec_(xSec),
                     norm_(norm)
{

  // load json file
  json j;
  std::ifstream f(filename);
  j << f;

  // load atributtes from JSON
  taggers_ = j["taggers"].get<std::vector<std::string>>();
  workPoints_ = j["workPoints"].get<std::vector<std::vector<double>>>();
  ptBins_ = j["ptBins"].get<std::vector<double>>();
  etaBins_ = j["etaBins"].get<std::vector<double>>();
  nEventPass_ = j["nEventPass"].get<std::vector<double>>();
  good_cat_jets_ = j["good_cat_jets"].get<std::vector<double>>();
  tag_cat_jets_ = j["tag_cat_jets"].get<std::vector<std::vector<std::vector<double>>>>();
  tagMultiplicity_ = j["tagMultiplicity"].get<std::vector<std::vector<std::vector<double>>>>();
  // Because get directly does not work
  for (auto it = j["cat_counts"].begin(); it != j["cat_counts"].end(); ++it) {
    cat_counts_[it.key()] = it.value().get<std::vector<double>>();
  }

  // default mapping
  for (std::size_t i = 0; i < good_cat_jets_.size(); i++) {
    cat_mapping_.emplace_back(1,i);
  }

} 

void Component::set_tag_wp(std::string tag, double wp) {

  // set tagger index
  auto it_tag = std::find(taggers_.begin(), taggers_.end(), tag);
  i_tag_ = std::distance(taggers_.begin(), it_tag);
  // set wp index
  auto it_wp = std::find(workPoints_[i_tag_].begin(), workPoints_[i_tag_].end(), wp);
  i_wp_ = std::distance(workPoints_[i_tag_].begin(), it_wp);

}


std::vector<double> Component::get_pretag_eff() const {
  std::vector<double> eff;
  // efficency
  eff.emplace_back(nEventPass_[1]/nEventGen_);
  // error
  double m = nEventPass_[0];
  double N = nEventGen_;
  eff.emplace_back(std::sqrt(m*(1-m/N))/N);
  return eff;
}


std::vector<double> Component::get_good_cat_jets() const {

  std::vector<double> good_cat_jets;
  for ( const auto & mapping : cat_mapping_ ) {
    good_cat_jets.emplace_back(0.0);
    for (const auto & cat : mapping ) {
      good_cat_jets.back()+=good_cat_jets_.at(cat);
    } 
  }
  return good_cat_jets;
}

std::vector<double> Component::get_tag_cat_jets() const {

  std::vector<double> tag_cat_jets;
  for ( const auto & mapping : cat_mapping_ ) {
    tag_cat_jets.emplace_back(0.0);
    for (const auto & cat : mapping ) {
      tag_cat_jets.back()+=tag_cat_jets_.at(i_tag_).at(i_wp_).at(cat);
    } 
  }
  return tag_cat_jets;
}

std::vector<double> Component::get_tag_multiplicity() const {
  return tagMultiplicity_.at(i_tag_).at(i_wp_);
}
    
std::map<std::string, std::vector<double>> Component::get_cat_fractions() const {

  std::map<std::string, std::vector<double>> cat_fractions; 

  for ( const auto & cat : cat_counts_) {
    // convert to new category mapping
    std::string key = "";
    std::vector<double> value = {0.0, 0.0};
    for ( const auto & mapping : cat_mapping_ ) {
      key.push_back('0');
      for (const auto & i_cat : mapping ) {
        key.back()+=(cat.first.at(i_cat)-'0');
        value.at(0) += cat.second.at(0);
        value.at(1) += cat.second.at(1);
      } 
    }
    // write to map (divide between total number)
    if (cat_fractions.count(key) > 0) {
      cat_fractions.at(key).at(0) += value.at(0);
      cat_fractions.at(key).at(1) += value.at(1);
    } else {
      cat_fractions[key] = value;
    }
  }

  // convert to fractions
  for ( auto & cat : cat_fractions ) {
    cat.second.at(0) /= nEventPass_.at(1) ;
    // error TODO
  }

  return cat_fractions;
}; 

void Model::add_mc_component(std::string filename, double nEventGen,
                             double xSec, Norm norm) {
  mc_comps_.emplace_back(filename, nEventGen, xSec, norm);
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



