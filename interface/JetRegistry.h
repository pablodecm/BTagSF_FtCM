
#pragma once

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream> 
#include <limits>

#include <TH2.h>

// mut_dataformats includes
#include "mut_framework/mut_dataformats/interface/Jet.h"
#include "mut_framework/mut_utils/interface/json.h"

// class to manage the results of the FtCM JetCounter selector
// this will be written to a root or JSON file
class JetRegistry {
 
   public:

    typedef std::vector<std::vector<TH2D>> TagTH2D;
    typedef std::vector<std::vector<std::vector<double>>> TagVector;
    typedef std::vector<std::vector<int>> TagNumber;
    typedef std::string KinematicCategory;
    typedef std::vector<std::string> FlavourCategory;

    // [unweighted_events, weighthed_events, sumw2]
    std::vector<double> nEventPass_ = {0.0, 0.0, 0.0};   

     // vector of taggers to be used
    std::vector<std::string> taggers_;
    // vector of working points for each tagger
    std::vector<std::vector<double>> workPoints_;
    // vector with ptBins
    std::vector<double> ptBins_;
    // vector with etaBins
    std::vector<double> etaBins_;

    // this is to count the tagged jet multiplicity
    std::vector<std::vector<std::vector<double>>> tagMultiplicity_;
    // also pre tag jet multiplicity
    std::vector<double> jetMultiplicity_;
    
    // vector to count good jets and tagged jets of each category
    std::vector<double> good_cat_jets_;
    std::vector<std::vector<std::vector<double>>> tag_cat_jets_;
 
    // histograms for the jets that pass selection (no btagging) 
    TH2D good_jets_; 
    std::vector<double> good_b_jets_;
    std::vector<double> good_c_jets_;
    std::vector<double> good_l_jets_;
    std::vector<double> good_x_jets_;
    // vector of vector of histograms for tagged jets for each tagger and WP
    TagTH2D tag_jets_;
    TagVector tag_b_jets_;
    TagVector tag_c_jets_;
    TagVector tag_l_jets_;
    TagVector tag_x_jets_;

    // event counts for each event category pair <[x,l,c,b1,...,bn],[counts, sumw2]>
    std::map<KinematicCategory, std::map<std::string,std::vector<double>>> cat_counts_;   

    // constructor
    JetRegistry(const std::vector<std::string> & taggers,
                const std::vector<std::vector<double>> & workPoints,
                const std::vector<double> & ptBins,
                const std::vector<double> & etaBins);

    // default constructor (required by ROOT)
    JetRegistry() {}

    // destructor
    ~JetRegistry() {}

    // save jet multiplicity (pretag)
    void registerJetMultiplicity(const int & nGoodJets, double eWeight);

    // add a jet to all the histograms and returns the category of the jet
    std::string registerJet( const mut::Jet & jet,
                     TagNumber & tagNumber,
                     double eWeight = 1.);  
    // count event in the corresponding category (return true if created)
    bool registerEvent( const KinematicCategory & kin_cat,
                        const FlavourCategory & flav_cat,
                        const TagNumber & tagNumber,
                        double weight = 1.);  

    void serialize(std::ostream & os);
    void serialize(std::string filename);
    
};


