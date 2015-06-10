
#include "../interface/JetRegistry.h"

JetRegistry::JetRegistry(const std::vector<std::string> & taggers,
                         const std::vector<std::vector<double>> & workPoints,
                         const std::vector<double> & ptBins,
                         const std::vector<double> & etaBins) :
                         taggers_(taggers),
                         workPoints_(workPoints),
                         ptBins_(ptBins),
                         etaBins_(etaBins), 
                         nCat_((ptBins.size()-1)*(etaBins.size()-1)+3),
                         good_jets_("good_jets","", ptBins.size()-1, ptBins.data(), etaBins.size()-1, etaBins.data()), 
                         good_b_jets_("good_b_jets","", ptBins.size()-1, ptBins.data(), etaBins.size()-1, etaBins.data()), 
                         good_c_jets_("good_c_jets","", ptBins.size()-1, ptBins.data(), etaBins.size()-1, etaBins.data()), 
                         good_l_jets_("good_l_jets","", ptBins.size()-1, ptBins.data(), etaBins.size()-1, etaBins.data()), 
                         good_x_jets_("good_x_jets","", ptBins.size()-1, ptBins.data(), etaBins.size()-1, etaBins.data()) 
{

  // init vector of histograms
  for (std::size_t t=0; t < taggers_.size(); t++) {
    tag_jets_.emplace_back();
    tag_b_jets_.emplace_back();
    tag_c_jets_.emplace_back();
    tag_l_jets_.emplace_back();
    tag_x_jets_.emplace_back();
    tagMultiplicity_.emplace_back();
    for (std::size_t i=0; i < workPoints_[t].size(); i++) {
      std::string tag_jets_name = "tag_jets_"+taggers_[t]+"_"+std::to_string(i);
      tag_jets_.back().emplace_back(tag_jets_name.c_str(),"",
                                    ptBins_.size()-1, ptBins_.data(),
                                    etaBins_.size()-1, etaBins_.data());
      std::string tag_b_jets_name = "tag_b_jets_"+taggers_[t]+"_"+std::to_string(i);
      tag_b_jets_.back().emplace_back(tag_b_jets_name.c_str(),"",
                                      ptBins_.size()-1, ptBins_.data(),
                                      etaBins_.size()-1, etaBins_.data());
      std::string tag_c_jets_name = "tag_c_jets_"+taggers_[t]+"_"+std::to_string(i);
      tag_c_jets_.back().emplace_back(tag_c_jets_name.c_str(),"",
                                      ptBins_.size()-1, ptBins_.data(),
                                      etaBins_.size()-1, etaBins_.data());
      std::string tag_l_jets_name = "tag_l_jets_"+taggers_[t]+"_"+std::to_string(i);
      tag_l_jets_.back().emplace_back(tag_l_jets_name.c_str(),"",
                                      ptBins_.size()-1, ptBins_.data(),
                                      etaBins_.size()-1, etaBins_.data());

      std::string tag_x_jets_name = "tag_x_jets_"+taggers_[t]+"_"+std::to_string(i);
      tag_x_jets_.back().emplace_back(tag_x_jets_name.c_str(),"",
                                      ptBins_.size()-1, ptBins_.data(),
                                      etaBins_.size()-1, etaBins_.data());
      // init multiplicity counting (max 10 b jet multiplicity)
      tagMultiplicity_.back().emplace_back(10, 0.0);
      
    }
  }
}

int JetRegistry::registerJet( const mut::Jet & jet,
                              TagNumber & tagNumber) {

  // get TH2D global bins and axis bins (1..)
  int global_bin = good_jets_.Fill( jet.pt() , jet.eta()); 
  int bin_pt = global_bin%(ptBins_.size()+1);  
  int bin_eta = ((global_bin-bin_pt)/(ptBins_.size()+1))%(etaBins_.size()+1);  
  // index
  int cat_index = (bin_pt-1)+(bin_eta-1)*(etaBins_.size()-1);

  
  int jet_flavour = jet.getPartonFlavour();
 
  if ( jet_flavour == 5) { // b jets
    good_b_jets_.Fill( jet.pt(), jet.eta());
    cat_index += 3; // O is x, 1 is l and 2 is c
  } else if ( jet_flavour == 4 ) { // c jets
    good_c_jets_.Fill( jet.pt(), jet.eta());
    cat_index = 2; // O is x, 1 is l and 2 is c
  } else if ( jet_flavour == 1 || jet_flavour == 2 || 
              jet_flavour == 3 || jet_flavour == 21 ) { // light jets
    good_l_jets_.Fill( jet.pt(), jet.eta());
    cat_index = 1; // O is x, 1 is l and 2 is c
  } else { //unknown jets
    good_x_jets_.Fill( jet.pt(), jet.eta());
    cat_index = 0; // O is x, 1 is l and 2 is c
  }
 
  for (std::size_t t = 0; t < taggers_.size(); t++) {
    for (std::size_t i = 0; i < workPoints_[t].size(); i++) {
      // check if tagged
      bool isTagged = jet.getDiscriminator(taggers_[t]) > workPoints_[t][i];  
      if (isTagged) {
        tagNumber[t][i]++;
        tag_jets_[t][i].Fill( jet.pt(), jet.eta());
        if ( jet_flavour == 5) {
          tag_b_jets_[t][i].Fill( jet.pt() , jet.eta()); 
        } else  if ( jet_flavour == 4) {
          tag_c_jets_[t][i].Fill( jet.pt() , jet.eta()); 
        } else if ( jet_flavour == 1 || jet_flavour == 2 || 
                    jet_flavour == 3 || jet_flavour == 21 ) { 
          tag_l_jets_[t][i].Fill( jet.pt(), jet.eta());
        } else { //unknown jets
          tag_x_jets_[t][i].Fill( jet.pt(), jet.eta());
        }
      }
    }
  }

  return cat_index;
}

bool JetRegistry::registerEvent( const ShortIntVector & cat, 
                                 const TagNumber & tagNumber,
                                 double weight)
{
  nEventPass_++;

  // update categoty map
  bool key_existed = false;
  if (cat_counts_.count(cat) > 0) {
    key_existed = true;
    cat_counts_[cat][0] += weight;
  } else {
    cat_counts_[cat] = {weight,0.0};
  }

  for (std::size_t t = 0; t < taggers_.size(); t++) {
    for (std::size_t i = 0; i < workPoints_[t].size(); i++) {
      tagMultiplicity_[t][i][tagNumber[t][i]] += weight;
    }
  }


  return key_existed;
}

std::ostream & extractCatCounts( std::ostream & output, const JetRegistry & jetRegistry) {

  for (const auto& catPair : jetRegistry.cat_counts_) {
    for (std::size_t i=0; i < catPair.first.size(); i++)  {
      if ( i == 0 )
      { 
        output << "[" << +catPair.first[i];  
      } else if ( i == catPair.first.size() - 1) {  
        output << "," << +catPair.first[i] << "]";  
      } else {
        output << "," << +catPair.first[i];  
      }
    }
    output << " = " << std::setw(15) << std::right;
    output << catPair.second[0] << "\u00B1";
    output << std::setw(15) << std::right << catPair.second[1] << std::endl;
  }

  return output;
}



