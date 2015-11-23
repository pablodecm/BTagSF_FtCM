#define SingleLeptonSelector_cxx
#include "../interface/SingleLeptonSelector.h"


// start of query (executed on client)
void SingleLeptonSelector::Begin(TTree * /*tree*/)
{

   std::string option = GetOption();

}

// right after begin (executed on slave)
void SingleLeptonSelector::SlaveBegin(TTree * /*tree*/)
{
   std::string  option = GetOption();

   std::size_t i_ofile = option.find("ofile="); 
   if (i_ofile != std::string::npos) {
     std::size_t length = (option.find(";", i_ofile) -  option.find("=", i_ofile) - 1);
     o_filename = option.substr(option.find("=", i_ofile)+1 , length );
   } else {
     o_filename = "output.root";
   }

   std::cout << "Output filename: " << o_filename << std::endl;
   n_entries = 0;

   o_file = new TFile(o_filename.c_str(), "RECREATE");

  ttree = new TTree("tree","Physics Object based TTree");

  ttree->Branch("eventInfo","mut::EventInfo", &eventInfo_ptr, 64000,1);
  ttree->Branch("pfmet","mut::MET", &pfmet_ptr, 64000,1);
  ttree->Branch("pfjets","std::vector<mut::Jet>", &pfjets_ptr, 64000,1);
  ttree->Branch("muons","std::vector<mut::Lepton>", &muons_ptr, 64000,1);
  ttree->Branch("elecs","std::vector<mut::Lepton>", &elecs_ptr, 64000,1);

  ttree->AutoSave();



}

// for each entry of the TTree
Bool_t SingleLeptonSelector::Process(Long64_t entry)
{

  n_entries++;
  if ((n_entries%10000) == 0) std::cout << "processing " << n_entries << " entry" << std::endl; 


  // set TTreeReader entry
  fReader.SetEntry(entry);


  // basic event selection (only mu+jets channel)
  std::vector<int> good_jets_index;
  if (muons->size() != 1) return false; 
  // check good lumi filter (Run I)
  if (eventInfo->hasFilter("goodLumi")) {
    if (!eventInfo->getFilter("goodLumi")) {
      return false; // it has but no passes filter
    }
  }
  // check if MC pass any muon trigger 
  if (eventInfo->hasFilter("pass_any_muon_trig") && !eventInfo->isRealData()) {
    if (!eventInfo->getFilter("pass_any_muon_trig")) {
      return false; // it has but no passes filter
    }
  }

  bool pass_event_sel = false;
  if ( pfmet->Et() > 20.0 &&
       muons->at(0).pt() > 35.0 &&
       std::abs(muons->at(0).eta()) < 2.1 &&
       muons->at(0).getLeptonIso("relIso") < 0.125 ) 
  {
    // check number of jets in pt/eta region 
    for ( std::size_t j = 0; j < pfjets->size(); j++) {
      const auto & pfjet = pfjets->at(j);
      // fill will return -1 if not in pt/eta region
      if (std::abs(pfjet.eta()) < jet_max_eta_ ) { 
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

  // if event passes event selection save objects
  eventInfo_ptr = new mut::EventInfo(*eventInfo);
  pfmet_ptr = new mut::MET(*pfmet);
  muons_ptr = new std::vector<mut::Lepton>(*muons);
  // only save selected jets
  pfjets_ptr = new std::vector<mut::Jet>();

  // for each good jet index
  for (auto j : good_jets_index) {
    // reference to current jet
    const auto & good_jet = pfjets->at(j);
    pfjets_ptr->emplace_back(good_jet);
  }


  ttree->Fill();

  delete eventInfo_ptr;
  delete pfjets_ptr;
  delete pfmet_ptr;
  delete muons_ptr;

  return pass_event_sel;

}

// all entries have been processed (executed in slave)
void SingleLeptonSelector::SlaveTerminate()
{
  ttree->AutoSave();
  delete ttree;
  delete o_file;

}

// last function called (on client)
void SingleLeptonSelector::Terminate()
{

}

