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

   all_jets = new TH2D("all_jets","all_jets",
                       ptBins_.size() - 1, ptBins_.data(),
                       etaBins_.size() - 1, etaBins_.data());

   

}

// for each entry of the TTree
Bool_t JetCounter::Process(Long64_t entry)
{

  // set TTreeReader entry
  fReader.SetEntry(entry);

  // fill all_jets histogram
  for ( auto pfjet : *pfjets ) all_jets->Fill( pfjet.pt() , pfjet.eta());

  bool pass_event_sel = false;
  // basic event selection
  if ( pfmet->Et() >= 20.0 &&
       muons->at(0).pt() > 35.0 &&
       muons->at(0).eta() < 2.1 &&
       muons->at(0).getLeptonIso("relIso") < 0.125 &&
       eventInfo->getFilter("goodLumi") ) 
  {
    pass_event_sel = true;
  };

  std::cout << pass_event_sel << std::endl;


  // for each jet loop
  for (std::size_t j = 0; j < pfjets->size(); j++) {

    
    

  }

  return kTRUE;
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

