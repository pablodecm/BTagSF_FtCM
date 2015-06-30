
#include "../interface/Component.h"

Component::Component(std::string filename, double nEventGen,
                     double xSec, Normalization n = FIXED) :
                     c_name_(filename.substr(0,filename.find_last_of("."))),
                     nEventGen_(nEventGen),
                     xSec_(xSec),
                     n_(n)
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

