#define JetCounter_cxx
#include "../interface/JetCounter.h"


// start of query (executed on client)
void JetCounter::Begin(TTree * /*tree*/)
{

   TString option = GetOption();

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
      if (all_good_jets->Fill( pfjet.pt() , pfjet.eta()) > 0) 
        good_jets_index.emplace_back(j);
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

  // vector categories counts (init to zero)
  JetRegistry::ShortIntVector cat((ptBins_.size()-1)*(etaBins_.size()-1)+3,0);
  // number of tagged jets (init at zero)
  JetRegistry::TagNumber tagNumber;
  for (std::size_t t = 0; t < taggers_.size(); t++) {
      tagNumber.emplace_back();
    for (std::size_t i = 0; i < workPoints_[t].size(); i++) {
      tagNumber.back().emplace_back(0);
    }
  }


  // for each good jet index
  for (auto j : good_jets_index) {
    // reference to current jet
    const auto & good_jet = pfjets->at(j);
    int cat_index = jetRegistry_->registerJet(good_jet, tagNumber);
    cat[cat_index]++;
  }

  jetRegistry_->registerEvent(cat, tagNumber);

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

