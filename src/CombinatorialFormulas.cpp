
#include "../interface/CombinatorialFormulas.h"

namespace FtCM {

  std::vector<std::string> multiset(int n, int k) {
    if (k == 0) return std::vector<std::string>(1, std::string(n, '0'));
    if (n == 0) return std::vector<std::string>();
    if (n == 1) return std::vector<std::string>(1, std::string(1, char(k)+'0'));
  std::vector<std::string> first;
  for (auto comb : multiset(n-1,k))  first.push_back("0"+comb);
  std::vector<std::string> second;
  for (auto comb : multiset(n, k-1)) {
    std::string element(comb);
    element[0]++; 
    second.push_back(element);
  }
  first.insert( first.end(), second.begin(), second.end());
  return first; 
  }

  std::vector<std::string> submultiset(std::string m, int k) {
    std::vector<std::string> mset = multiset(m.size(), k);
    std::vector<std::string> smset; 
    for (auto comb : mset) {
      bool in_subset = true;
      for (std::size_t i_c=0; i_c < m.size(); i_c++) {
        if ( comb[i_c] < '0' || comb[i_c] > m[i_c] ) {
          in_subset = false;
        } 
      } 
      if (in_subset) smset.push_back(comb);
    }  
    return smset;
  }


  std::vector<std::vector<std::string>> cart_product( const std::vector<std::vector<std::string>> & v) {
    std::vector<std::vector<std::string>> s = {{}};
    for (auto& u : v) {
        std::vector<std::vector<std::string>> r;
        for (auto& x : s) {
            for (auto y : u) {
                r.push_back(x);
                r.back().push_back(y);
            }
        }
        s.swap(r);
    }
    return s;
  }


  std::vector<std::string> all_tag_cats(const std::string & pretag_cat) {

    std::vector<std::vector<std::string>> list_comb;
    
    for (const char & c: pretag_cat) {
      list_comb.emplace_back();
      char cc = c;  
      while ( cc != '0') { 
        std::string s;  
        s.push_back(cc);
        list_comb.back().emplace_back(s);
        cc--;  
      }
      std::string s;  
      s.push_back(cc);
      list_comb.back().emplace_back(s);
    }

    std::vector<std::vector<std::string>> all_comb = cart_product(list_comb);
    std::vector<std::string> tag_cats;
    for (const auto & v_tag_cats : all_comb) {
      std::string tag_cat;
      for (const auto & n : v_tag_cats) {
        tag_cat += n;
      }
      tag_cats.emplace_back(tag_cat);
    }

    return tag_cats;
  }
}
  
