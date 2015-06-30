
#pragma once

#include <vector>
#include <map>
#include <fstream>
#include <string>

#include "mut_framework/mut_utils/interface/json.h"

enum Normalization {
  SIGNAL,
  FLOATING,
  FIXED,
  BKG,
  DATA
};


class Component {
  
  public:

    std::string c_name_;

    // [unweighted_events, weighthed_events, sumw2]
    std::vector<double> nEventPass_;   
    double nEventGen_;
    double xSec_;
    Normalization n_;
    

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
    
    // vector to count good jets and tagged jets of each category
    std::vector<double> good_cat_jets_;
    std::vector<std::vector<std::vector<double>>> tag_cat_jets_;

    // event counts for each event category pair <[x,l,c,b1,...,bn],[counts, sumw2]>
    std::map<std::string, std::vector<double>> cat_counts_;   

    // category mapping (default is [[1],[2],...] )
    std::vector<std::vector<int>> cat_mapping_; 

    // tag and wqp indexes
    std::size_t i_tag_;
    std::size_t i_wp_;


    Component() {}
    Component(std::string filename, double nEventGen, double xSec, Normalization n );
    ~Component() {}

    void set_category_mapping( std::vector<std::vector<int>> cat_mapping) 
      { cat_mapping_ = cat_mapping; }
    void set_tag_wp( std::string tag, double wp); 

    std::size_t get_n_cat() const { return cat_mapping_.size(); }
    std::vector<double> get_pretag_eff() const;
    std::vector<double> get_good_cat_jets() const;
    std::vector<double> get_tag_cat_jets() const;
    std::vector<double> get_tag_multiplicity() const;
    std::map<std::string, std::vector<double>> get_cat_fractions() const;

};


