
#include "../interface/Component.h"

Component::Component(std::string filename, std::string tagger,
                     double workPoint, double nEventGen) :
                     nEventGen_(nEventGen),
                     tagger_(tagger),
                     workPoint_(workPoint)
{

  if (nEventGen < 0) isData_ = true;

  // load json file
  json j;
  std::ifstream f(filename);
  j << f;

  // load atributtes from JSON
  ptBins_ = j["ptBins"].get<std::vector<double>>();
  etaBins_ = j["etaBins"].get<std::vector<double>>();
  nEventPass_ = j["nEventPass"].get<std::vector<double>>();

   // Because get directly does not work
   for (auto it = j["cat_counts"].begin(); it != j["cat_counts"].end(); ++it) {
     pretag_cat_counts_[it.key()] = std::map<std::string,std::vector<double>>();
     for (auto itt = it->begin(); itt != it->end(); ++itt) {
       pretag_cat_counts_[it.key()][itt.key()] = itt.value().get<std::vector<double>>();
     }
   }
 
   // Because get directly does not work
   std::string one_tag_name = tagger_+":"+std::to_string(workPoint_);
   for (auto it = j["tag_cat_counts-"+one_tag_name].begin();
        it != j["tag_cat_counts-"+one_tag_name].end(); ++it) {
     tag_cat_counts_[it.key()] = std::map<std::string,std::vector<double>>();
     for (auto itt = it->begin(); itt != it->end(); ++itt) {
       tag_cat_counts_[it.key()][itt.key()] = itt.value().get<std::vector<double>>();
     }
   }
   for (auto it = j["pretag_jet_counts-"+one_tag_name].begin(); 
        it != j["pretag_jet_counts-"+one_tag_name].end(); ++it) {
     pretag_jet_counts_[it.key()] = std::map<std::string,std::vector<double>>();
     for (auto itt = it->begin(); itt != it->end(); ++itt) {
       pretag_jet_counts_[it.key()][itt.key()] = itt.value().get<std::vector<double>>();
     }
   }
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

std::vector<double> Component::get_good_jets(
    const std::vector<int> & type,
    const std::set<std::pair<std::string,std::string>> & cat_set ) const
{
  std::vector<double> good_jets((ptBins_.size()-1)*(etaBins_.size()-1), 0.0);
  for (const auto & cat : cat_set) {
    if ( pretag_jet_counts_.count(cat.first) > 0) {
      if (pretag_jet_counts_.at(cat.first).count(cat.second) > 0) {
        const std::vector <double> & pretag_jet_count = pretag_jet_counts_.at(cat.first).at(cat.second);
        for (std::size_t b = 0; b < (ptBins_.size()-1)*(etaBins_.size()-1); b++) {
          for (std::size_t t = 0; t < type.size(); t++) {
            good_jets.at(b) += pretag_jet_count.at(b*4 +type.at(t))/nEventGen_;
          }
        }
      }
    }
  } 
  return good_jets;
}

std::vector<double> Component::get_tag_jets(
    const std::vector<int> & type,
    const std::set<std::pair<std::string,std::string>> & cat_set ) const
{
  std::vector<double> tag_jets((ptBins_.size()-1)*(etaBins_.size()-1), 0.0);
  for (const auto & cat : cat_set) {
      if ( tag_cat_counts_.count(cat.first) > 0) {
      const auto & tag_cat_count = tag_cat_counts_.at(cat.first);
      for ( const auto & tag_cat : tag_cat_count) {
        std::string s_tag_kin_cat((ptBins_.size()-1)*(etaBins_.size()-1),'0');
        for (std::size_t b = 0; b < (ptBins_.size()-1)*(etaBins_.size()-1); b++) {
          int sum = 0;
          for (std::size_t t=0; t < 4; t++) {
            sum += int(tag_cat.first.at(4*b+t)-'0');
          }
          s_tag_kin_cat.at(b) =  char(sum)+'0';
        }
        
        if (cat.second == s_tag_kin_cat) {
          for (std::size_t b = 0; b < (ptBins_.size()-1)*(etaBins_.size()-1); b++) {
            for (std::size_t t = 0; t < type.size(); t++) {
              tag_jets.at(b) += double(int(tag_cat.first.at(4*b+type.at(t))-'0'))*tag_cat.second[0]/nEventGen_;
            }
          }
        }
      }
    } 
  }

  return tag_jets;
}


double Component::get_tag_counts( const std::string & pretag_cat,
    const std::string & tag_cat) const {
  double counts = 0.0;
  if ( tag_cat_counts_.count(pretag_cat) > 0) {
  const std::map<std::string,std::vector<double>> & sub_map = tag_cat_counts_.at(pretag_cat);
  for (const auto & kv : sub_map) {
    const std::string & extended_cat = kv.first;
    std::string short_cat((ptBins_.size()-1)*(etaBins_.size()-1),'0');
    for (std::size_t b = 0; b < (ptBins_.size()-1)*(etaBins_.size()-1); b++) {
      int jet_sum = 0;
      for (std::size_t t = 0; t < 4; t++) {
        jet_sum += int(extended_cat.at(4*b+t)-'0');
      }
      short_cat.at(b) = char(jet_sum) + '0';
    }
    // check if short cat matches and add to counts
    if (short_cat == tag_cat) counts  += kv.second.at(0);
  }
  }
  return counts;
}

double Component::get_pretag_counts( const std::string & pretag_cat) const {

  double counts = 0.0;
  if ( tag_cat_counts_.count(pretag_cat) > 0) {
  const std::map<std::string,std::vector<double>> & sub_map = tag_cat_counts_.at(pretag_cat);
    for (const auto & kv : sub_map) {
      counts  += kv.second.at(0);
    }
  }
  return counts;
}

std::map<std::string, std::vector<double>> Component::get_flav_frac(
    const std::string & pretag_cat,
    std::vector<std::vector<int>> cat_mapping ) const {

  std::map<std::string, std::vector<double>> flav_frac; 
  if ( pretag_cat_counts_.count(pretag_cat) > 0) {
  const std::map<std::string,std::vector<double>> & sub_map = pretag_cat_counts_.at(pretag_cat);
  for ( const auto & cat : sub_map) {
    // convert to new category mapping
    std::string key = "";
    std::vector<double> value = {0.0, 0.0};
    for (std::size_t b = 0; b < (ptBins_.size()-1)*(etaBins_.size()-1); b++) {
      for ( const auto & mapping : cat_mapping ) {
        key.push_back('0');
        for (const auto & i_cat : mapping ) {
          key.back()+=(cat.first.at(b*4 + i_cat)-'0');
        }
      }
    }
    value.at(0) += cat.second.at(0);
    value.at(1) += cat.second.at(1);
    // write to map (divide between total number)
    if (flav_frac.count(key) > 0) {
      flav_frac.at(key).at(0) += value.at(0);
      flav_frac.at(key).at(1) += value.at(1);
    } else {
      flav_frac[key] = value;
    }
  }
  }

  // convert to efficiencies 
  for ( auto & cat : flav_frac ) {
    cat.second.at(0) /= nEventGen_ ;
    // error TODO
  }
  return flav_frac;
} 

