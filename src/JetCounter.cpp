#define JetCounter_cxx
#include "../interface/JetCounter.h"


// start of query (executed on client)
void JetCounter::Begin(TTree * /*tree*/)
{

   std::string option = GetOption();
   std::size_t i_isData = option.find("isData"); 
   if (i_isData != std::string::npos) {
     isData_ = true;
   } else {
     isData_ = false;
   }
}

// right after begin (executed on slave)
void JetCounter::SlaveBegin(TTree * /*tree*/)
{

   TString option = GetOption();

   all_good_jets = new TH2D("all_good_jets","good_jets",
                            ptBins_.size() - 1, ptBins_.data(),
                            etaBins_.size() - 1, etaBins_.data());

   jetRegistry_ = new JetRegistry(taggers_, workPoints_, ptBins_, etaBins_); 

}

// for each entry of the TTree
Bool_t JetCounter::Process(Long64_t entry)
{

  // set TTreeReader entry
  fReader.SetEntry(entry);

  // categories counts (init to zero)
  JetRegistry::KinematicCategory kin_cat((ptBins_.size()-1)*(etaBins_.size()-1),'0');
  JetRegistry::FlavourCategory flav_cat((ptBins_.size()-1)*(etaBins_.size()-1),"0000");
  // number of tagged jets (init at zero)
  JetRegistry::TagKinematicCategory tag_kin_cat;
  for (std::size_t t = 0; t < taggers_.size(); t++) {
      tag_kin_cat.emplace_back();
    for (std::size_t i = 0; i < workPoints_[t].size(); i++) {
      tag_kin_cat.back().emplace_back((ptBins_.size()-1)*(etaBins_.size()-1),"0000");
    }
  }

  // get event weight ( it will be 1 for data)
  double weight = getEventWeight();

  // fill jet multiplicity
  jetRegistry_->registerJetMultiplicity(pfjets->size(), weight);


  for ( std::size_t j = 0; j < pfjets->size(); j++) {
    // reference to current jet
    const auto & good_jet = pfjets->at(j);
    // fill will return -1 if not in pt/eta region
    if (all_good_jets->Fill( good_jet.pt() , good_jet.eta()) < 0) { 
      std::cout << "WARNING: some jets are outside of pt/eta ranges" << std::endl;
    }
    std::string jet_type = jetRegistry_->registerJet(good_jet, tag_kin_cat, weight);
    kin_cat[jet_type[1]]++; // increment kin category
    flav_cat[jet_type[1]][jet_type[0]]++; // increment flavour for kinematic cat
  }


  jetRegistry_->registerEvent(kin_cat, flav_cat, tag_kin_cat, weight);

  return true;

}

// all entries have been processed (executed in slave)
void JetCounter::SlaveTerminate()
{

}

// last function called (on client)
void JetCounter::Terminate()
{

}

void JetCounter::addTagger( std::string name, double min, double max, int num) {

  taggers_.emplace_back(name);
  workPoints_.emplace_back(num);

  // set working points for the tagger
  double step = (max - min)/double(num);
  for (int i=0; i< num; i++) workPoints_.back().at(i) = min + step*double(i);

}

void JetCounter::addTagger( std::string name, std::vector<double> workPoints) {

  taggers_.emplace_back(name);
  workPoints_.emplace_back(workPoints);

}

double JetCounter::getEventWeight() {
  
  double weight = 1.0;

  if (!isData_) {
    for ( auto eWeight : eWeights_ ) {
      weight*=eventInfo->getWeight(eWeight);
    }
  }   

  return weight;
}


