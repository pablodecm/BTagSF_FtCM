
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
    }
  }
}

void JetRegistry::registerJetMultiplicity(const int & nGoodJets, double eWeight) {
  jetMultiplicity_[nGoodJets] += eWeight;
} 

int JetRegistry::registerJet( const mut::Jet & jet,
                              TagNumber & tagNumber,
                              double eWeight) {

  // get TH2D global bins and axis bins (1..)
  int global_bin = good_jets_.Fill( jet.pt() , jet.eta());
  int bin_pt = global_bin%(ptBins_.size()+1);
  int bin_eta = ((global_bin-bin_pt)/(ptBins_.size()+1))%(etaBins_.size()+1);
  // index
  int cat_index = (bin_pt-1)+(bin_eta-1)*(etaBins_.size()-1);
  // global index
  int glob_index = cat_index;

  int jet_flavour = jet.getPartonFlavour();

  if ( jet_flavour == 5) { // b jets
    good_b_jets_[cat_index] += eWeight;
    glob_index += 3; // O is x, 1 is l and 2 is c
  } else if ( jet_flavour == 4 ) { // c jets
    good_c_jets_[cat_index] += eWeight;
    glob_index = 2; // O is x, 1 is l and 2 is c
  } else if ( jet_flavour == 1 || jet_flavour == 2 ||
              jet_flavour == 3 || jet_flavour == 21 ) { // light jets
    good_l_jets_[cat_index] += eWeight;
    glob_index = 1; // O is x, 1 is l and 2 is c
  } else { //unknown jets
    good_x_jets_[cat_index] += eWeight;
    glob_index = 0; // O is x, 1 is l and 2 is c
  }

  // add to good jets count
  good_cat_jets_[glob_index] += eWeight;


  for (std::size_t t = 0; t < taggers_.size(); t++) {
    for (std::size_t i = 0; i < workPoints_[t].size(); i++) {
      // check if tagged
      bool isTagged = jet.getDiscriminator(taggers_[t]) > workPoints_[t][i];
      if (isTagged) {
        tagNumber[t][i]++;
        tag_jets_[t][i].Fill( jet.pt(), jet.eta());
        tag_cat_jets_[t][i][glob_index] += eWeight;
        if ( jet_flavour == 5) {
          tag_b_jets_[t][i][cat_index] += eWeight;
        } else  if ( jet_flavour == 4) {
          tag_c_jets_[t][i][cat_index] += eWeight;
        } else if ( jet_flavour == 1 || jet_flavour == 2 ||
                    jet_flavour == 3 || jet_flavour == 21 ) {
          tag_l_jets_[t][i][cat_index] += eWeight;
        } else { //unknown jets
          tag_x_jets_[t][i][cat_index] += eWeight;
        }
      }
    }
  }

  return glob_index;
}

bool JetRegistry::registerEvent( const JetCategory & cat,
                                 const TagNumber & tagNumber,
                                 double weight)
{
  nEventPass_[0]++;
  nEventPass_[1]+=weight;
  nEventPass_[2]+=weight*weight;

  // update category map
  bool key_existed = false;
  if (cat_counts_.count(cat) > 0) {
    key_existed = true;
    cat_counts_[cat][0] += weight;
    cat_counts_[cat][1] += weight*weight;
  } else {
    cat_counts_[cat] = {weight,weight*weight};
  }

  for (std::size_t t = 0; t < taggers_.size(); t++) {
    for (std::size_t i = 0; i < workPoints_[t].size(); i++) {
      tagMultiplicity_[t][i][tagNumber[t][i]] += weight;
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

  os << std::setw(4) << j << std::endl;

  }


