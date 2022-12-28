#ifndef __analysis_pteros__
#define __analysis_pteros__

#include <pteros/pteros.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <complex> // library for complex numbers


typedef struct { // Used to input the center of masses for each lipid
  double dist;
  int top_phi_C12E2_count;
  int top_phi_C12E2M_count;
  int bot_phi_C12E2_count;
  int bot_phi_C12E2M_count;
} phi_struct;

typedef struct { // Used to identify the group and distance to compute the
                 // orderphobic effect
  float resid_self;
  float resid_other;
  float dist;
  float self_x_coord;
  float self_y_coord;
  float self_z_coord;
  float x_coord;
  float y_coord;
  float z_coord;
  float orderphobic_value_self = 0.0;
} oph_struct;


// functions to work with 
bool compare_by_distance(const oph_struct&, const oph_struct&);

// classes to work with 
class pteros_traj_analysis {
public:
  pteros_traj_analysis(std::string, std::string); // constructor 
  ~pteros_traj_analysis(); // destructor
  // void methodologies
  void compute_order_parameter();
  void compute_budding_rate();
  //void compute_composition_around_np();
  void compute_grid(int);
  void compute_resid(int); // new 

  // orderphobic parameter collection and computation
  void collect_orderphobic_parameter_vector();
  void compute_orderphobic_parameter();
  // protected methods are accesible in the class that defines them and in classes
  // that inherit from that class.
  double cartesian_distance(float&, float&, float&, float&,
			    float&, float&);

  // print functions
  void print_structs(int);
  
  
private:  // private members are only accessible within the class defining them 
  
  std::vector<std::vector<oph_struct> > orderphobic_storage;
  std::set<int> resid;
  std::set<int> surfactant_ids;
  std::set<int>::iterator set_itr;
  // unordered map - the equivalent of a dictionary
  std::unordered_map<int, std::vector<int> > resid_indices;
  pteros::System sample_system;
};

#endif
