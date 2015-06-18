
#include "../interface/JetRegistry.h"

JetRegistry::JetRegistry(const std::vector<std::string> & taggers,
                         const std::vector<std::vector<double>> & workPoints,
                         const std::vector<double> & ptBins,
                         const std::vector<double> & etaBins) :
                         taggers_(taggers),
                         workPoints_(workPoints),
                         ptBins_(ptBins),
                         etaBins_(etaBins), 
                         good_cat_jets_(3+(ptBins_.size()-1)*(etaBins_.size()-1), 0.0),
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
    tag_cat_jets_.emplace_back(); 
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
      // init tag cat jet counts
      tag_cat_jets_.back().emplace_back(3+(ptBins_.size()-1)*(etaBins_.size()-1), 0.0); 
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
 
  // add to good jets count
  good_cat_jets_[cat_index]++;


  for (std::size_t t = 0; t < taggers_.size(); t++) {
    for (std::size_t i = 0; i < workPoints_[t].size(); i++) {
      // check if tagged
      bool isTagged = jet.getDiscriminator(taggers_[t]) > workPoints_[t][i];  
      if (isTagged) {
        tagNumber[t][i]++;
        tag_jets_[t][i].Fill( jet.pt(), jet.eta());
        tag_cat_jets_[t][i][cat_index]++; 
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
  nEventPass_[0]++;

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

void JetRegistry::serialize( std::string filename ) {
  
  std::ofstream output(filename, std::ofstream::out); 
  serialize(output);
  output.flush();
  output.close();  
}

void JetRegistry::serialize( std::ostream & os) {

  os << "{\n"; // first json brace

  // serialize tagger and workpoints
  os << "\t\"taggers\": [\n";
  for (std::size_t i = 0; i < taggers_.size() - 1; i ++) 
    os << "\t\t\"" << taggers_.at(i) << "\",\n";
  os << "\t\t\"" << taggers_.back() << "\"\n\t],\n";
  os << "\t\"workPoints\": [\n";
  for (std::size_t i = 0; i < workPoints_.size() - 1; i++) {
    os << "\t\t[";
    for (std::size_t ii = 0; ii < workPoints_.at(i).size() - 1; ii++) 
      os << workPoints_.at(i).at(ii) << ", ";
    os << workPoints_.at(i).back() << "],\n";
  } 
  os << "\t\t[";
  for (std::size_t ii = 0; ii < workPoints_.back().size() - 1; ii++) 
      os << workPoints_.back().at(ii) << ", ";
  os << workPoints_.back().back() << "]\n";
  os << "\t],\n"; // end of workpoints serialization

  // serialize pt and eta bin categories
  os << "\t\"ptBins\": [";
  for (std::size_t i = 0; i < ptBins_.size() - 1; i ++) 
    os << ptBins_.at(i) << ", ";
  os << ptBins_.back() << "],\n";
  os << "\t\"etaBins\": [";
  for (std::size_t i = 0; i < etaBins_.size() - 1; i ++) 
    os << etaBins_.at(i) << ", ";
  os << etaBins_.back() << "],\n";

  // number of events passed pre-tag cuts
  os << "\t\"nEventPass\": [" << nEventPass_[0] << ", "<< nEventPass_[1] << ", "
     << nEventPass_[2] << "],\n";

  // serialize good jets per category coounts
  os << "\t\"good_cat_jets\": [";
  for (std::size_t i = 0; i < good_cat_jets_.size() - 1; i ++) 
    os << good_cat_jets_.at(i) << ", ";
  os << good_cat_jets_.back() << "],\n";

  // serialize tagged jets per category counts
  os << "\t\"tag_cat_jets\": [\n";
  for (std::size_t i = 0; i < tag_cat_jets_.size() - 1; i++) {
    os << "\t\t[\n";
    for (std::size_t ii = 0; ii < tag_cat_jets_.at(i).size() - 1; ii++) {
      os << "\t\t\t[";
      for (std::size_t iii = 0; iii < tag_cat_jets_.at(i).at(ii).size() - 1; iii++) {
        os << tag_cat_jets_.at(i).at(ii).at(iii) << ", ";
      } 
      os << tag_cat_jets_.at(i).at(ii).back()<< "],\n";
    }
    os << "\t\t\t[";
    for (std::size_t iii = 0; iii < tag_cat_jets_.at(i).back().size() - 1; iii++) {
        os << tag_cat_jets_.at(i).back().at(iii) << ", ";
    }
    os << tag_cat_jets_.at(i).back().back()<< "]\n\t\t],\n";
  }
  os << "\t\t[\n";
  for (std::size_t ii = 0; ii < tag_cat_jets_.back().size() - 1; ii++) {
    os << "\t\t\t[";
    for (std::size_t iii = 0; iii < tag_cat_jets_.back().at(ii).size() - 1; iii++) {
      os << tag_cat_jets_.back().at(ii).at(iii) << ", ";
    } 
    os << tag_cat_jets_.back().at(ii).back()<< "],\n";
  }
  os << "\t\t\t[";
  for (std::size_t iii = 0; iii < tag_cat_jets_.back().back().size() - 1; iii++) {
    os << tag_cat_jets_.back().back().at(iii) << ", ";
  }
  os << tag_cat_jets_.back().back().back()<< "]\n\t\t]\n\t],\n";


  // serialize tagMultiplicity
  os << "\t\"tagMultiplicity\": [\n";
  for (std::size_t i = 0; i < tagMultiplicity_.size() - 1; i++) {
    os << "\t\t[\n";
    for (std::size_t ii = 0; ii < tagMultiplicity_.at(i).size() - 1; ii++) {
      os << "\t\t\t[";
      for (std::size_t iii = 0; iii < tagMultiplicity_.at(i).at(ii).size() - 1; iii++) {
        os << tagMultiplicity_.at(i).at(ii).at(iii) << ", ";
      } 
      os << tagMultiplicity_.at(i).at(ii).back()<< "],\n";
    }
    os << "\t\t\t[";
    for (std::size_t iii = 0; iii < tagMultiplicity_.at(i).back().size() - 1; iii++) {
        os << tagMultiplicity_.at(i).back().at(iii) << ", ";
    }
    os << tagMultiplicity_.at(i).back().back()<< "]\n\t\t],\n";
  }
  os << "\t\t[\n";
  for (std::size_t ii = 0; ii < tagMultiplicity_.back().size() - 1; ii++) {
    os << "\t\t\t[";
    for (std::size_t iii = 0; iii < tagMultiplicity_.back().at(ii).size() - 1; iii++) {
      os << tagMultiplicity_.back().at(ii).at(iii) << ", ";
    } 
    os << tagMultiplicity_.back().at(ii).back()<< "],\n";
  }
  os << "\t\t\t[";
  for (std::size_t iii = 0; iii < tagMultiplicity_.back().back().size() - 1; iii++) {
    os << tagMultiplicity_.back().back().at(iii) << ", ";
  }
  os << tagMultiplicity_.back().back().back()<< "]\n\t\t]\n\t],\n";

  // serialize category counts
  os << "\t\"cat_counts\": {\n";
  for (auto catPair= cat_counts_.begin(); catPair!= --cat_counts_.end(); catPair++) {
    // write category as skey string (e.g. "(0,0,1,3)") 
    os << "\t\t\"(" << +catPair->first.front();  
    for (std::size_t i=1; i < catPair->first.size() - 1; i++)
        os << "," << +catPair->first[i];  
    os << "," << +catPair->first.back() << ")\"";  
    // write [weighted_counts, error] as value
    os << ": [" << catPair->second[0] << ", " << catPair->second[1] << "],\n" ;
  }
  // last element is special case
  auto last =  --cat_counts_.end();
  os << "\t\t\"(" << +last->first.front();  
  for (std::size_t i=1; i < last->first.size() - 1; i++)
    os << "," << +last->first[i];  
  os << "," << +last->first.back() << ")\"";  
  os << ": [" << last->second[0] << ", " << last->second[1] << "]\n";
  os << "\t}"; // end of cat counts serialization

  os << "\n}\n"; // last json brace
}


