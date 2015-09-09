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


  // basic event selection
  std::vector<int> good_jets_index;
  bool pass_event_sel = false;
  if ( pfmet->Et() >= 20.0 &&
       muons->at(0).pt() > 35.0 &&
       std::abs(muons->at(0).eta()) < 2.1 &&
       muons->at(0).getLeptonIso("relIso") < 0.125 &&
       eventInfo->getFilter("goodLumi") ) 
  {
    // check number of jets in pt/eta region 
    for ( std::size_t j = 0; j < pfjets->size(); j++) {
      const auto & pfjet = pfjets->at(j);
      // fill will return -1 if not in pt/eta region
      if (all_good_jets->Fill( pfjet.pt() , pfjet.eta()) > 0) { 
        // check if each jet passes the corresponding pt criteria
        if (good_jets_index.size() < min_pt_jets_.size()) {
          if (pfjet.pt() >=  min_pt_jets_.at(good_jets_index.size())) {
            good_jets_index.emplace_back(j);
          }
        } else {
          if (pfjet.pt() >=  min_pt_jets_.back()) {
            good_jets_index.emplace_back(j);
          }
        } 
      }
    } 

    // pass event selection if 4 or more good jets
    if (good_jets_index.size() >= 4) {
      pass_event_sel = true;
    } else {
      return false; // no enough good jets
    }
  } else {
    return false; // no passes event selection
  }

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
  jetRegistry_->registerJetMultiplicity(good_jets_index.size(), weight);

  // for each good jet index
  for (auto j : good_jets_index) {
    // reference to current jet
    const auto & good_jet = pfjets->at(j);
    std::string jet_type = jetRegistry_->registerJet(good_jet, tag_kin_cat, weight);
    kin_cat[jet_type[1]]++; // increment kin category
    flav_cat[jet_type[1]][jet_type[0]]++; // increment flavour for kinematic cat
  }


  jetRegistry_->registerEvent(kin_cat, flav_cat, tag_kin_cat, weight);

  return pass_event_sel;

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
    if (useOldMuonSF_) { // apply muon SF as a function of eta
      double muon_eta = muons->at(0).eta();
      if (std::abs(muon_eta) <= 0.9) {weight *= 0.994*0.993*0.976;}
      else if (std::abs(muon_eta) > 0.9 && std::abs(muon_eta) <= 1.2 ) 
        {weight *= 0.992*0.998*0.961;}
      else if (std::abs(muon_eta) > 1.2 && std::abs(muon_eta) <= 2.1 )
        {weight *= 0.998*1.002*0.983; }
    }
  }   

  return weight;
}


