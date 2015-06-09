
#pragma once

#include <vector>
#include <string>
#include <map>
#include <iostream>

#include <TH2.h>

// mut_dataformats includes
#include "mut_framework/mut_dataformats/interface/Jet.h"

// class to manage the results of the FtCM JetCounter selector
// this will be written to a root file
class JetRegistry {
 
   public:

    typedef std::vector<std::vector<TH2D>> TagTH2D;
    typedef unsigned char ShortInt;  
    typedef std::vector<ShortInt> ShortIntVector;  

    
    unsigned long int nEventPass_ = 0;   

     // vector of taggers to be used
    std::vector<std::string> taggers_;
    // vector of working points for each tagger
    std::vector<std::vector<double>> workPoints_;
    // vector with ptBins
    std::vector<double> ptBins_;
    // vector with etaBins
    std::vector<double> etaBins_;
    // number of caterories
    unsigned int nCat_;
 
    // histograms for the jets that pass selection (no btagging) 
    TH2D good_jets_; 
    TH2D good_b_jets_; 
    TH2D good_c_jets_; 
    TH2D good_l_jets_; 
    TH2D good_x_jets_; 
    // vector of vector of histograms for tagged jets for each tagger and WP
    TagTH2D tag_jets_;
    TagTH2D tag_b_jets_;
    TagTH2D tag_c_jets_;
    TagTH2D tag_l_jets_;
    TagTH2D tag_x_jets_;

    // event counts for each event category pair<[x,l,c,b1,...,bn],[counts, error]>
    std::map<ShortIntVector, std::vector<double>> cat_counts_;   

    // constructor
    JetRegistry(const std::vector<std::string> & taggers,
                const std::vector<std::vector<double>> & workPoints,
                const std::vector<double> & ptBins,
                const std::vector<double> & etaBins);

    // destructor
    ~JetRegistry() {}

    // add a jet to all the histograms and returns the category of the jet
    int registerJet( const mut::Jet & jet);  
    // count event in the corresponding category (return true if created)
    bool registerEvent( const ShortIntVector & cat, double weight = 1.);  
    
    friend std::ostream& extractCatCounts( std::ostream & out, const JetRegistry & jetRegistry); 
     
};
