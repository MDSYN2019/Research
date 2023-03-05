#include "pteros/pteros.h"
#include <string>
#include <iostream>
#include <vector>


pteros::System load_system(std::string input_data) {
  /*
    sample system loading function
  */ 
  pteros::System sample_system;
  sample_system(input_data);
  return sample_system;
}



int main() {
  pteros::System sys1;
  sys1.load("/home/sang/Desktop/CG_NP_ANALYSIS/CG/NP/HYDROPHOBIC/1_NP_sim/1_NP_sim_NPT.gro");
  pteros::Selection sel0 = sys1("all");
  std::cout << sel0 << std::endl;
  return 0;
}
