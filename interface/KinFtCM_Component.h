
#pragma once

#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <set>

#include "mut_framework/mut_utils/interface/json.h"


namespace KinFtCM {

  enum Norm {
    SIGNAL,
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
      Norm n_;
      

      // tagger to use
      std::string tagger_;
      // working point to use
      double workPoint_;
      // vector with ptBins
      std::vector<double> ptBins_;
      // vector with etaBins
      std::vector<double> etaBins_;

      // flavour counts for each pretag category 
      std::map<std::string, std::map<std::string,std::vector<double>>> pretag_cat_counts_;   
      // counts for each pretag-tag category (with flavour info for MC)
      std::map<std::string, std::map<std::string,std::vector<double>>> tag_cat_counts_;
      // pretag jet counts for each pretag-tag category 
      std::map<std::string, std::map<std::string,std::vector<double>>> pretag_jet_counts_;

      Component() {}
      Component(std::string filename, std::string tagger,
                double workPoint, double nEventGen, double xSec, Norm n );
      ~Component() {}

      std::string get_name() const { return c_name_; }
      std::size_t get_n_cat() const { return ((ptBins_.size()-1)*(etaBins_.size()-1)); }

      std::vector<double> get_pretag_eff() const;
      std::vector<double> get_good_jets(const std::vector<int> & type,
                                        const std::set<std::pair<std::string,std::string>> & cat_set ) const;
      std::vector<double> get_tag_jets(const std::vector<int> & type,
                                       const std::set<std::pair<std::string,std::string>> & cat_set ) const;

      double get_tag_counts( const std::string & pretag_cat, const std::string & tag_cat) const;

      std::map<std::string, std::vector<double>> get_flav_frac(const std::string & pretag_cat,
          std::vector<std::vector<int>> cat_mapping= {{0}, {1}, {2,3}}) const; 

  };

}

