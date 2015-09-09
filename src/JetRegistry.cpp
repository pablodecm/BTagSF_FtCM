
#include "../interface/JetRegistry.h"

JetRegistry::JetRegistry(const std::vector<std::string> & taggers,
                         const std::vector<std::vector<double>> & workPoints,
                         const std::vector<double> & ptBins,
                         const std::vector<double> & etaBins) :
                         taggers_(taggers),
                         workPoints_(workPoints),
                         ptBins_(ptBins),
                         etaBins_(etaBins),
                         jetMultiplicity_(10, 0.0),
                         good_cat_jets_(3+(ptBins_.size()-1)*(etaBins_.size()-1), 0.0),
                         good_jets_("good_jets","", ptBins.size()-1, ptBins.data(), etaBins.size()-1, etaBins.data()),
                         good_b_jets_((ptBins_.size()-1)*(etaBins_.size()-1), 0.0),
                         good_c_jets_((ptBins_.size()-1)*(etaBins_.size()-1), 0.0),
                         good_l_jets_((ptBins_.size()-1)*(etaBins_.size()-1), 0.0),
                         good_x_jets_((ptBins_.size()-1)*(etaBins_.size()-1), 0.0)
{

  // init vector of histograms
  for (std::size_t t=0; t < taggers_.size(); t++) {
    tag_jets_.emplace_back();
    tag_b_jets_.emplace_back();
    tag_c_jets_.emplace_back();
    tag_l_jets_.emplace_back();
    tag_x_jets_.emplace_back();
    tagMultiplicity_.emplace_back();
    tag_cat_jets_.emplace_back();
    tag_cat_counts_.emplace_back();
    pretag_jet_counts_.emplace_back();
    for (std::size_t i=0; i < workPoints_[t].size(); i++) {
      std::string tag_jets_name = "tag_jets_"+taggers_[t]+"_"+std::to_string(i);
      tag_jets_.back().emplace_back(tag_jets_name.c_str(),"",
                                    ptBins_.size()-1, ptBins_.data(),
                                    etaBins_.size()-1, etaBins_.data());
      tag_b_jets_.back().emplace_back((ptBins_.size()-1)*(etaBins_.size()-1), 0.0);
      tag_c_jets_.back().emplace_back((ptBins_.size()-1)*(etaBins_.size()-1), 0.0);
      tag_l_jets_.back().emplace_back((ptBins_.size()-1)*(etaBins_.size()-1), 0.0);
      tag_x_jets_.back().emplace_back((ptBins_.size()-1)*(etaBins_.size()-1), 0.0);
      // init multiplicity counting (max 10 jet multiplicity)
      tagMultiplicity_.back().emplace_back(10, 0.0);
      // init tag cat jet counts
      tag_cat_jets_.back().emplace_back(3+(ptBins_.size()-1)*(etaBins_.size()-1), 0.0);
      tag_cat_counts_.back().emplace_back();
      pretag_jet_counts_.back().emplace_back();
    }
  }
}

void JetRegistry::registerJetMultiplicity(const int & nGoodJets, double eWeight) {
  jetMultiplicity_[nGoodJets] += eWeight;
} 

std::string JetRegistry::registerJet( const mut::Jet & jet,
                              TagKinematicCategory & tag_kin_cat,
                              double eWeight) {

  // get TH2D global bins and axis bins (1..)
  int global_bin = good_jets_.Fill( jet.pt() , jet.eta());
  int bin_pt = global_bin%(ptBins_.size()+1);
  int bin_eta = ((global_bin-bin_pt)/(ptBins_.size()+1))%(etaBins_.size()+1);
  // index (start at 0)
  int cat_index = (bin_pt-1)+(bin_eta-1)*(etaBins_.size()-1);

  // global index
  int glob_index = cat_index;

  // init jet type as unknown and char 47 (error) 
  std::string jet_type = "x/";
  jet_type.at(1) = char(cat_index); 

  int jet_flavour = jet.getPartonFlavour();

  if ( jet_flavour == 5) { // b jets
    good_b_jets_[cat_index] += eWeight;
    jet_type.at(0) = char(0);
  } else if ( jet_flavour == 4 ) { // c jets
    good_c_jets_[cat_index] += eWeight;
    jet_type.at(0) = char(1);
  } else if ( jet_flavour == 1 || jet_flavour == 2 ||
              jet_flavour == 3 || jet_flavour == 21 ) { // light jets
    good_l_jets_[cat_index] += eWeight;
    jet_type.at(0) = char(2);
  } else { //unknown jets
    good_x_jets_[cat_index] += eWeight;
    jet_type.at(0) = char(3);
  }

  glob_index += int(jet_type.at(0));

  // add to good jets count
  good_cat_jets_[glob_index] += eWeight;


  for (std::size_t t = 0; t < taggers_.size(); t++) {
    for (std::size_t i = 0; i < workPoints_[t].size(); i++) {
      // check if tagged
      bool isTagged = jet.getDiscriminator(taggers_[t]) > workPoints_[t][i];
      if (isTagged) {
//tag_jets_[t][i].Fill( jet.pt(), jet.eta());
        tag_cat_jets_.at(t).at(i).at(glob_index) += eWeight;
        if ( jet_flavour == 5) {
          tag_b_jets_.at(t).at(i).at(cat_index) += eWeight;
          tag_kin_cat.at(t).at(i).at(cat_index).at(0)++;
        } else  if ( jet_flavour == 4) {
          tag_c_jets_.at(t).at(i).at(cat_index) += eWeight;
          tag_kin_cat.at(t).at(i).at(cat_index).at(1)++;
        } else if ( jet_flavour == 1 || jet_flavour == 2 ||
                    jet_flavour == 3 || jet_flavour == 21 ) {
          tag_l_jets_.at(t).at(i).at(cat_index) += eWeight;
          tag_kin_cat.at(t).at(i).at(cat_index).at(2)++;
        } else { //unknown jets
          tag_x_jets_.at(t).at(i).at(cat_index) += eWeight;
          tag_kin_cat.at(t).at(i).at(cat_index).at(3)++;
        }
      }
    }
  }

  return jet_type;
}

bool JetRegistry::registerEvent( const KinematicCategory & kin_cat,
                                 const FlavourCategory & flav_cat,
                                 const TagKinematicCategory & tag_kin_cat,
                                 double weight)
{
  nEventPass_[0]++;
  nEventPass_[1]+=weight;
  nEventPass_[2]+=weight*weight;

  // concatenate flavour category (so it can be a JSON key)
  std::string j_flav_cat = "";
  for ( const auto & flav_count : flav_cat) j_flav_cat += flav_count;

  // update category map
  bool key_existed = false;
  if (cat_counts_.count(kin_cat) > 0) {
    if (cat_counts_[kin_cat].count(j_flav_cat) > 0) {
      key_existed = true;
      cat_counts_[kin_cat][j_flav_cat][0] += weight;
      cat_counts_[kin_cat][j_flav_cat][1] += weight*weight;
    } else {
      cat_counts_[kin_cat][j_flav_cat] = {weight,weight*weight};
    }
  } else {
    cat_counts_[kin_cat] = std::map<std::string,std::vector<double>>();
    cat_counts_[kin_cat][j_flav_cat] = {weight,weight*weight}; 
  }

    for (std::size_t t = 0; t < taggers_.size(); t++) {
    for (std::size_t i = 0; i < workPoints_[t].size(); i++) {
      
      // update tag category map
      CategoryCounts & one_tag_cat_counts = tag_cat_counts_.at(t).at(i); 
      // concatenate flavour category (so it can be a JSON key)
      std::string j_flav_tag_cat = "";
      for ( const auto & flav_count : tag_kin_cat.at(t).at(i)) j_flav_tag_cat += flav_count;
      int n_tags = 0;
      for (const auto & ch : j_flav_tag_cat) n_tags+= int(ch-'0'); 
      tagMultiplicity_.at(t).at(i).at(n_tags) += weight;

      if (one_tag_cat_counts.count(kin_cat) > 0) {
        if (one_tag_cat_counts[kin_cat].count(j_flav_tag_cat) > 0) {
          one_tag_cat_counts[kin_cat][j_flav_tag_cat][0] += weight;
          one_tag_cat_counts[kin_cat][j_flav_tag_cat][1] += weight*weight;
        } else {
          one_tag_cat_counts[kin_cat][j_flav_tag_cat] = {weight,weight*weight};
        }
      } else {
        one_tag_cat_counts[kin_cat] = std::map<std::string,std::vector<double>>();
        one_tag_cat_counts[kin_cat][j_flav_tag_cat] = {weight,weight*weight}; 
      }

      std::vector<double> ev_jet_count(4*(ptBins_.size()-1)*(etaBins_.size()-1), 0.0);
      for (std::size_t n_cat = 0; n_cat < flav_cat.size(); n_cat++) {
        for (std::size_t j_t = 0; j_t < flav_cat.at(n_cat).size(); j_t++) {
          ev_jet_count.at(4*n_cat+j_t) = double(int(flav_cat.at(n_cat).at(j_t)-'0'))*weight;
        }
      }
      std::string s_tag_kin_cat((ptBins_.size()-1)*(etaBins_.size()-1),'0');
      for (std::size_t n_cat = 0; n_cat < tag_kin_cat.size(); n_cat++) {
        int sum = 0;
        for ( const auto & cat : tag_kin_cat.at(t).at(i).at(n_cat)) {
          sum += int(cat-'0');
        }
        s_tag_kin_cat.at(n_cat) =  char(sum)+'0';
      }

      CategoryCounts & one_pretag_jet_counts = pretag_jet_counts_.at(t).at(i); 
      if (one_pretag_jet_counts.count(kin_cat) > 0) {
        if (one_pretag_jet_counts[kin_cat].count(s_tag_kin_cat) > 0) {
          for (std::size_t n=0; n<ev_jet_count.size(); n++) 
            one_pretag_jet_counts[kin_cat][s_tag_kin_cat].at(n) += ev_jet_count.at(n);
        } else {
          one_pretag_jet_counts[kin_cat][s_tag_kin_cat] = ev_jet_count;
        }
      } else {
        one_pretag_jet_counts[kin_cat] = std::map<std::string,std::vector<double>>();
        one_pretag_jet_counts[kin_cat][s_tag_kin_cat] = ev_jet_count; 
      }
    }
  }

  return key_existed;
}

void JetRegistry::serialize( std::string filename ) {

  std::ofstream output(filename, std::ofstream::out);
  serialize(output);
  output.flush();
  output.close();
}

void JetRegistry::serialize( std::ostream & os) {

  // set precision to maximum double precison
  os << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  // use modern JSON C++ library
  json j;

  j["taggers"] = taggers_;
  j["workPoints"] = workPoints_;
  j["ptBins"] = ptBins_;
  j["etaBins"] = etaBins_;
  j["nEventPass"] = nEventPass_;
  j["good_cat_jets"] = good_cat_jets_;
  j["good_b_jets"] = good_b_jets_;
  j["good_c_jets"] = good_c_jets_;
  j["good_l_jets"] = good_l_jets_;
  j["good_x_jets"] = good_x_jets_;
  j["tag_cat_jets"] = tag_cat_jets_;
  j["tag_b_jets"] = tag_b_jets_;
  j["tag_c_jets"] = tag_c_jets_;
  j["tag_l_jets"] = tag_l_jets_;
  j["tag_x_jets"] = tag_x_jets_;
  j["tagMultiplicity"] = tagMultiplicity_;
  j["jetMultiplicity"] = jetMultiplicity_;
  j["cat_counts"] = cat_counts_;

  for (std::size_t t = 0; t < taggers_.size(); t++) {
    for (std::size_t i = 0; i < workPoints_[t].size(); i++) {
      CategoryCounts tag_counts = tag_cat_counts_.at(t).at(i);
      j["tag_cat_counts-"+taggers_[t]+":"+std::to_string(workPoints_[t][i])] = tag_counts; 
      CategoryCounts pretag_counts = pretag_jet_counts_.at(t).at(i);
      j["pretag_jet_counts-"+taggers_[t]+":"+std::to_string(workPoints_[t][i])] = pretag_counts; 
    }
  }

  os << std::setw(4) << j << std::endl;

  }


