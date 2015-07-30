
#ifndef JetCounter_h
#define JetCounter_h

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
#include <TH2.h>

// mut_dataformats includes
#include "mut_framework/mut_dataformats/interface/EventInfo.h"
#include "mut_framework/mut_dataformats/interface/Lepton.h"
#include "mut_framework/mut_dataformats/interface/Jet.h"
#include "mut_framework/mut_dataformats/interface/MET.h"

// FtCM JetRegistry
#include "mut_framework/BTagSF_FtCM/interface/JetRegistry.h"


class JetCounter : public TSelector {
public :
  
   // asociated with a TTree 
   TTreeReader fReader;

   // is actual data
   bool isData_;
   // use 2012 simplified Muon SFs
   bool useOldMuonSF_;

   // asociated with a TBranch
   TTreeReaderValue<mut::EventInfo> eventInfo;
   TTreeReaderValue<std::vector<mut::Lepton>> muons;
   TTreeReaderValue<std::vector<mut::Jet>> pfjets;
   TTreeReaderValue<mut::MET> pfmet;

   // vector of taggers to be used
   std::vector<std::string> taggers_;
   // vector of working points for each tagger
   std::vector<std::vector<double>> workPoints_;

   // pt and eta bins
   std::vector<double> ptBins_;
   std::vector<double> etaBins_;

   // pt-eta histogram of good jets
   TH2D * all_good_jets;

   // event weights to take into account (set defaults)
   std::vector<std::string> eWeights_ = 
    { "weight_PU_nom" };

   // class to manage and save results 
   JetRegistry * jetRegistry_;
    
   // default constructor
   JetCounter(TTree * /*tree*/ =0) :
      eventInfo(fReader, "eventInfo"),
      muons(fReader, "muons"),
      pfjets(fReader, "pfjets"),
      pfmet(fReader, "pfmet") { }
   // destructor
   virtual ~JetCounter() { }

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

   void setTaggers( std::vector<std::string> taggers ) { taggers_ = taggers; }
   void addTagger( std::string name, double min, double max, int num);
   void addTagger( std::string name, std::vector<double> workPoints);

   void setPtBins( std::vector<double> ptBins ) { ptBins_ = ptBins; }
   void setEtaBins( std::vector<double> etaBins ) { etaBins_ = etaBins; }

   void setEventWeights( std::vector<std::string> eWeights ) { eWeights_ = eWeights; }
   double getEventWeight();
   void useOldMuonSF( bool value) { useOldMuonSF_ = value; } 

   void serialize(std::string filename) { jetRegistry_->serialize(filename); } 
   void resetJetRegistry() { 
     delete jetRegistry_;
     delete all_good_jets; 
   }

};

#endif

#ifdef JetCounter_cxx

// each new tree is opened
void JetCounter::Init(TTree *tree)
{
  fReader.SetTree(tree);
}

// each new file is opened
Bool_t JetCounter::Notify()
{

   return kTRUE;
}

#endif // #ifdef JetCounter_cxx
