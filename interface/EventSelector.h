
#ifndef EventSelector_h
#define EventSelector_h

#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <memory>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

// mut_dataformats includes
#include "mut_framework/mut_dataformats/interface/EventInfo.h"
#include "mut_framework/mut_dataformats/interface/Lepton.h"
#include "mut_framework/mut_dataformats/interface/Jet.h"
#include "mut_framework/mut_dataformats/interface/MET.h"



class EventSelector : public TSelector {
public :

   long n_entries;
  
   // asociated with a TTree 
   TTreeReader fReader;

   // asociated with a TBranch
   TTreeReaderValue<mut::EventInfo> eventInfo;
   TTreeReaderValue<std::vector<mut::Lepton>> muons;
   TTreeReaderValue<std::vector<mut::Jet>> pfjets;
   TTreeReaderValue<mut::MET> pfmet;

    // mut objects to save to TTree
    mut::EventInfo * eventInfo_ptr = nullptr;
    std::vector<mut::Jet> * pfjets_ptr = nullptr;
    std::vector<mut::Lepton> * muons_ptr = nullptr;
    std::vector<mut::Lepton> * elecs_ptr = nullptr;
    mut::MET * pfmet_ptr = nullptr;

    // output filename
    std::string o_filename;
    // output TTree pointer
    TTree * ttree;
    // output TFile poinyer
    TFile * o_file;

   // mimimum pt for each good jet
   std::vector<double> min_pt_jets_ = { 70.0, 50.0, 30.0};
   double jet_max_eta_ = 2.4;

   // default constructor
   EventSelector(TTree * /*tree*/ =0) :
      eventInfo(fReader, "eventInfo"),
      muons(fReader, "muons"),
      pfjets(fReader, "pfjets"),
      pfmet(fReader, "pfmet") { }
   // destructor
   virtual ~EventSelector() { }

   // TSelector functions
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   void set_min_pt_jets( std::vector<double> min_pt_jets) { min_pt_jets_ = min_pt_jets; }

};

#endif

#ifdef EventSelector_cxx

// each new tree is opened
void EventSelector::Init(TTree *tree)
{
  fReader.SetTree(tree);
}

// each new file is opened
Bool_t EventSelector::Notify()
{

   return kTRUE;
}

#endif // #ifdef EventSelector_cxx
