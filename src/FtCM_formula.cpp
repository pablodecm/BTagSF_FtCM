
#include <vector>
#include <string>


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
    for (auto comp : mset) {
      bool in_subset = true;
      for (std::size_t i_c=0; i_c < m.size(); i_c++) {
        if ( comp[i_c] < '0' || comp[i_c] > m[i_c] ) {
          in_subset = false;
        } 
      } 
      if (in_subset) smset.push_back(comp);
    }  
    return smset;
  }

}
  
