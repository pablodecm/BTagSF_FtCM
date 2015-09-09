
#pragma once

#include <vector>
#include <string>

namespace FtCM {


  std::vector<std::string> multiset(int n, int k); 

  std::vector<std::string> submultiset(std::string m, int k); 

  std::vector<std::vector<std::string>> cart_product( const std::vector<std::vector<std::string>> & v); 
}
