#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// boost libraries
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

// custom library
#include "analysis.h"
#include "analysis_pteros.h"

int main(void) {
  // edit forcefield files 
  //analysis initial("/home/sang/Desktop/git/Research/C12E2_work/parameters/parm.cg_cmm_fake");
  //initial.read_file();
  //initial.generate_parameter();

  pteros_traj_analysis pteros_data(
      "/home/sang/Desktop/git/Research/C12E2_work/trajectories/dump.gro",
      "/home/sang/Desktop/git/Research/C12E2_work/trajectories/dump.xtc");

  //pteros_data.compute_order_parameter();
  //pteros_data.compute_composition_around_np();
  //pteros_data.compute_grid(100);
  //pteros_data.compute_order_parameter();
  pteros_data.compute_budding_rate();
  pteros_data.collect_orderphobic_parameter_vector();
  pteros_data.print_structs(50);
  
  return 0;
}
