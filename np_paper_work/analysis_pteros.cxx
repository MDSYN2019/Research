
#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>
// dictionary equivalent
#include <complex> // library for complex numbers
#include <unordered_map>

// pteros library and custom class import
#include "analysis_pteros.h"
#include <pteros/core/grid.h>
#include <pteros/pteros.h>

// filling in functions

// random new comment for the sake of using magit


bool compare_by_distance(const oph_struct &a, const oph_struct &b) {
  /*
   */
  return a.dist < b.dist;
}

// filling in classes

pteros_traj_analysis::pteros_traj_analysis(std::string input_data_gro,
                                           std::string input_data_xtc) {
  /*
    read the gro and xtc files
  */

  sample_system.load(input_data_gro);
  sample_system.load(input_data_xtc);
  std::cout << sample_system.num_frames() << std::endl;
};

void pteros_traj_analysis::compute_order_parameter() {
  // identify the nanoparticle
  std::string id = "index ";
  std::string id2 = "index ";

  pteros::Selection sel0;
  pteros::Selection sel1;
  pteros::Selection nanoparticle =
      sample_system("name 14"); // nanoparticle selection

  // loop over each residue, compute the normal and compute
  // the general representative vector and compute the angle between
  for (auto &it : resid_indices) {

    sel0 = sample_system(id.append(std::to_string(it.second[0])));
    sel1 = sample_system(id2.append(std::to_string(it.second[6])));

    // for (auto it_new = begin(it.second); it_new != end (it.second); ++it_new)
    // {
    //   std::cout << it.first << " " << *it_new << std::endl;
    // }

    //    for (auto& it : ) {

    std::cout << it.first << " " << it.second[1] << " " << sel0[0].name()
              << " " // sel0.center << " "
              << sel1[0].name() << " "
              << std::endl; // sel1.center()  << std::endl;

    id = "index ";
    id2 = "index ";
  }

  // resid_indices;
};

void pteros_traj_analysis::compute_budding_rate() {
  /*
     compute the budding rate by comparing the z coordinate absolute
     distance between the nanoparticle center from the central z coordinate
     of the 'lower' bilayer
  */
  std::string resname = "resid ";
  pteros::Selection sel0;
  pteros::Selection nanoparticle =
      sample_system("name 14"); // nanoparticle selection
  pteros::Selection all_surfactants(
      sample_system, "name 7 6 3 4 13 12 10 9"); // surfactant selection

  // as surfactant_ids is a set, it will only store the unique
  // residue ids. Any duplicated id will not be stored at all and ignored

  for (pteros::Selection::iterator it = all_surfactants.begin();
       it != all_surfactants.end(); it++) {
    surfactant_ids.insert(it->resid());
  }

  // explain!!
  for (set_itr = surfactant_ids.begin(); set_itr != surfactant_ids.end();
       set_itr++) {

    resname.append(std::to_string(*set_itr));
    sel0 = sample_system(resname);
    resid_indices[*set_itr] = std::vector<int>();

    for (pteros::Selection::iterator it = sel0.begin(); it != sel0.end();
         it++) {
      // create vector entry for dictionary
      resid_indices[*set_itr].push_back(it->index());
    }
    resname = "resid "; // reset resname string
  }
}

// void pteros_traj_analysis::compute_composition_around_np() {
//   // identify the residues around the nanoparticle
//
//   std::string resid = "resid ";
//   std::string not_resid = "not resid ";
//
//   pteros::Selection sel0;
//   pteros::Selection sel1;
//
//   //pteros::Selection composition_np(
//   //    sample_system, "within 3.0 of (name 14) and name 7 6 3 4 13 12 10
//   9");
//   // iterator-based loop:
//
//   for (set_itr = resid.begin(); set_itr != resid.end(); set_itr++) {
//     not_resid.append(std::to_string(*set_itr));
//     sel0 = sample_system(not_resid);
//
//     //for (pteros::Selection::iterator it = sel0.begin();
//     //	 it !=sel0.end(); it++) {
//
//     //}
//   }
//
// };

void pteros_traj_analysis::compute_grid(int grid_unit) {
  /*
    create a grid of grid_unit * grid_unit * grid_unit
  */
  pteros::Grid g(grid_unit, grid_unit, grid_unit);
  g.populate_periodic(sample_system.select_all());

  for (int i = 0; i < grid_unit; ++i)
    for (int j = 0; j < grid_unit; ++j)
      for (int k = 0; k < grid_unit; ++k)
        std::cout << g.cell(i, j, k).size() << std::endl;
};

void pteros_traj_analysis::compute_orderphobic_parameter() {
  /*
    take data from orderphobic_storage and compute the
    orderphobic parameter for that particular surfactant

    phi_l = |(1/6) \sum_{allpairs} exp(6 * i * theta_{ij})|^2

    where theta is the angle between an arbitrary axis and a vector
    connecting tail-end particle i to tail-end particle j

  */
  
  std::complex<double> inital_complex(0.0, 0.0);
  // initalize a complex fraction 1/6 with a 0 in the complex plane
  std::complex<double> inital_sixth((double)1 / 6, 0.0);
  
  for (int i = 0; i <= orderphobic_storage.size(); i++) {
    for (std::vector<oph_struct>::iterator it =
	   orderphobic_storage.at(i).begin();
         it != orderphobic_storage.at(i).end(); ++it) {

      
    }
  }

};

void pteros_traj_analysis::collect_orderphobic_parameter_vector() {
  /*
    compute the orderphobic parameter based on an arbitrary
    vector pointing out and
  */
  double u, v;
  double dist;
  // assign arbitrary vector to compare with
  double reference_vector[2];
  // initialize a complex number with 0s
  // in both the real and complex plane
  std::complex<double> inital_complex(0.0, 0.0);
  // initalize a complex fraction 1/6 with a 0 in the complex plane
  std::complex<double> inital_sixth((double)1 / 6, 0.0);
  // make sure to exclude waters

  pteros::Selection new_system;
  pteros::Selection new_system_II;
  pteros::Selection adhoc_system;

  reference_vector[0] = 10.0;
  reference_vector[1] = 10.0;
  u = pow((pow(reference_vector[0], 2) + pow(reference_vector[1], 2)), 0.5);

  //     sample_system("name 14"); // nanoparticle selection
  std::string adhoc_resname = "not name 1 and resid ";
  oph_struct orderphobic_data;
  std::vector<oph_struct> single_storage;

  std::string resid_1 = "resid ";
  std::string resid_2 = "resid ";

  // where to store orderphobic data
  int count = 0;
  for (auto &it : resid_indices)

    if (count == 50) {
      break;
    } else {
      resid_1.append(std::to_string(it.first)); // get residue id
      new_system = sample_system(resid_1);

      for (auto &it_2 : resid_indices) {

        if (it.first != it_2.first) {

          resid_2.append(std::to_string(it_2.first)); // get residue id
          new_system_II = sample_system(resid_2);

          dist = this->cartesian_distance(
              new_system_II[4].x(), new_system_II[4].y(), new_system_II[4].z(),
              new_system[4].x(), new_system[4].y(), new_system[4].z());

          // std::cout << resid_1 << " " << resid_2 << " " <<
          // new_system_II[4].name() <<  " " << new_system[4].name() << " " <<
          // dist << std::endl;

          // --- storing the orderphobic data ---
          orderphobic_data.resid_self = it.first;
          orderphobic_data.resid_other = it_2.first;
          orderphobic_data.dist = dist;
          // store cartesian coordinates of the self
          orderphobic_data.self_x_coord = new_system_II[4].x();
          orderphobic_data.self_y_coord = new_system_II[4].y();
          orderphobic_data.self_z_coord = new_system_II[4].z();
          // store cartesian coordinates of the cooordinates away
          orderphobic_data.x_coord = new_system[4].x();
          orderphobic_data.y_coord = new_system[4].y();
          orderphobic_data.z_coord = new_system[4].z();
          orderphobic_data.orderphobic_value_self =
              0.0; // placeholder for the orderphobic value
          // append entry
          single_storage.push_back(orderphobic_data);
          // reset strings
          // std::cout << newcount << std::endl;
          resid_2 = "resid ";
        }
      }

      resid_1 = "resid ";
      // sort entry and then push_back
      std::sort(single_storage.begin(), single_storage.end(),
                compare_by_distance);

      orderphobic_storage.push_back(single_storage);
      // clear all entries for next iteration
      single_storage.clear();
      count += 1;
    }
};

void pteros_traj_analysis::print_structs(int limit = 5) {
  /*
    print out vector of vector entries for the orderphobic struct inputs
   */
  for (int i = 0; i <= orderphobic_storage.size(); i++) {
    for (std::vector<oph_struct>::iterator it =
             orderphobic_storage.at(i).begin();
         it != orderphobic_storage.at(i).end(); ++it) {
      std::cout << i << " " << it->resid_self << " " << it->resid_other << " "
                << it->dist << std::endl;
      // break if conditions
      if (it - orderphobic_storage.at(i).begin() == limit) {
        break;
      }
    }
  }
};

double pteros_traj_analysis::cartesian_distance(float &com_1x, float &com_1y,
                                                float &com_1z, float &com_2x,
                                                float &com_2y, float &com_2z) {
  /*
    compute absolute distance between 3d cartesian points
  */

  double cdist = pow((pow(com_1x - com_2x, 2) + pow(com_1y - com_2y, 2) +
                      pow(com_1z - com_2z, 2)),
                     0.5);
  return cdist;
};

pteros_traj_analysis::~pteros_traj_analysis() {
  /*
    class destructor
  */
  std::cout << "destructor called" << std::endl;
};