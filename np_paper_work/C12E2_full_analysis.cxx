/*

 ------------------------------------------------------------------------
 |                    Analysis of C12E2 Mixed Bilayers                  |
 ------------------------------------------------------------------------

 Author = "Sang Young Noh"

 Version = "0.0.1"

 Updated = 04/04/2023

 The papers referenced for this work is:

 1. "F. C. MacKintosh and S. A. Safran, Phase Separation and curvature of
 bilayer membranes, Physical Review E, 47, 2, 1993"

 2. "R. Lipowsky, Domain-induced budding of fluid membranes, Biophysical
 Society, 64, 1993, 1133 - 1138"

 3. "J. Wolff and S. Komura and D. Andelman, Budding of domains in mixed bilayer
 domains, 91, Physical Review E, 91, 012708, 2015"

 4. S. Katira and K. K. Madadapu and S. Vaikuntanathan and B. Smit and D.
 Chandler, Pre-transition effects mediate forces of assembly between
 transmembrane proteins, elife, 2016, 5, e13150

 5. Physical Considerations of the Organization of Inclusions in Lipid Bilayer
 Systems, Spring 2015, Shachi Katira

*/

#include <algorithm>
#include <bits/stdc++.h>
#include <boost/progress.hpp>
#include <cerrno>
#include <cmath>
#include <complex> // Library for complex numbers
#include <cstdlib>
#include <functional> // std::multiplies
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <map>
#include <numeric> // std::adjacent_differenc#include <cmath>
#include <tuple>
#include <utility>
#include <vector>
#include <chrono> 

#define _USE_MATH_DEFINES

template <typename T> class basic_stopwatch : T {
  typedef typename T BaseTimer;
public:
  // create, optionally start timing an activity
  explicit basic_stopwatch(bool start);
  explicit basic_stopwatch(char const* activity = "Stopwatch",
			   bool start=true);
  basic_stopwatch(std::ostream& log,
		  char const* activity="Stopwatch",
		  bool start=true);
  // stop and destroy a stopwatch
  ~basic_stopwatch();
  // get last lap time (time of last stop)
  unsigned LapGet() const;
  // predicate: return true if the stopwatch is running
  bool IsStarted() const;
  // show accumulated time, keep running, set/return lap
  unsigned Show(char const* event="show");
  // (re)start a stopwatch, set/return lap time
  unsigned Start(char const* event_name="start");
  // stop a running stopwatch, set/return lap time
  unsigned Stop(char const* event_name="stop");
private:
  // members
  char const* m_activity; // "activity" string
  unsigned m_lap;
  // lap time (time of last stop)
  std::ostream& m_log;
  // stream on which to log events
};


// Constants for use later with analysis
const int number_of_polymers =
    998; // The number of polymers of each type - C12E2 or mimic
const int number_of_atoms = 71313; // Total number of beads in the simulation
const int index_cg = 7;
const int number_of_snapshots = 100; /*The number of screenshots in the dump file*/
const int box_dim = 3;

typedef struct {
  int index[7];
} C12E2_skeleton;

typedef struct {
  int a;
  int b;
  double x;
  double y;
  double z;
} inputCoord;

typedef struct { // Used to input the center of masses for each lipid
  int index;
  double x;
  double y;
  double z;
} COMstruct;

typedef struct { 
  double dist;
  int topPhiC12E2Count;
  int topPhiC12E2MCount;
  int botPhiC12E2Count;
  int botPhiC12E2MCount;
} phiStruct;

typedef struct { 
  double phip;
  double phim;
  std::vector<double> phipVec;
  std::vector<double> phimVec;
} phipm;

typedef struct { 
  double X, Y, Z;
  double orderphobicVal;
  std::vector<int> surrounding_residue_indices;
} OPVal;

typedef struct { // Used to identify the group and distance to compute the
                 // orderphobic effect
  int index;
  double dist;
  double selfXcoord;
  double selfYcoord;
  double selfZcoord;
  double Xcoord;
  double Ycoord;
  double Zcoord;

} OPHstruct;

double calcAngle(std::vector<OPHstruct> *OPHInput) {

  double refVector[2];
  double dot;
  double newx;
  double newy;
  double u, v;

  std::complex<double> angle;

  double phi;

  refVector[0] = 10.0;

  refVector[1] = 10.0;

  std::complex<double> mycomplex(0.0, 0.0);
  std::complex<double> sixth((double)1 / 6, 0.0);

  double absval;

  for (unsigned int i = 0; i < OPHInput->size(); ++i) {

    newx = OPHInput->at(i).Xcoord;

    newy = OPHInput->at(i).Ycoord;

    u = pow((pow(refVector[0], 2) + pow(refVector[1], 2)), 0.5);

    v = pow((pow(newx, 2) + pow(newy, 2)), 0.5);

    dot = (newx * refVector[0]) + (newy * refVector[1]);

    std::complex<double> placeholder(
        0, (6 * 180 * M_PI * acos((dot / (u * v))))); // TODO
    
    angle = exp(placeholder);
    mycomplex += angle;

    //  std::cout << i << " " << absval << std::endl;
  }
  mycomplex = sixth * mycomplex;
  absval = pow(abs(mycomplex), 2);
  // std::cout  << absval << std::endl;
  return absval;
}

double trueDist(double *COM1x, double *COM1y, double *COM1z, double *COM2x,
                double *COM2y, double *COM2z) {
  double dist = pow((pow(*COM1x - *COM2x, 2) + pow(*COM1y - *COM2y, 2) +
                     pow(*COM1z - *COM2z, 2)),
                    0.5);
  return dist;
}

// Function used for sorting vector of structs
bool compareByIndex(const inputCoord &a, const inputCoord &b) {
  return a.a < b.a;
}

bool compareByIndexOPh(const OPHstruct &a, const OPHstruct &b) {
  return a.dist < b.dist;
}

void CenterOfMass(std::vector<int> *vec1, std::vector<int> *vec2,
                  std::vector<std::vector<inputCoord>> *inputVec,
                  std::vector<COMstruct> *COM1, std::vector<COMstruct> *COM2,
                  std::vector<std::vector<COMstruct>> *C12E2final,
                  std::vector<std::vector<COMstruct>> *C12E2Mfinal) {
  //    CenterOfMass(&C12E2IndexVector, &C12E2MIndexVector, &inputTotal,
  //    &C12E2COM, &C12E2MCOM, &C12E2TotalCOMArray, &C12E2MTotalCOMArray);

  double C12E2comX, C12E2McomX; // X coordinate COM
  double C12E2comY, C12E2McomY; // Y coordinate COM
  double C12E2comZ, C12E2McomZ; // Z coordinate cCOM

  COMstruct C12E2input;  // struct to store COM coordinates for C12E2
  COMstruct C12E2Minput; // struct to store cOM coordinates for C12E2M

  for (unsigned int i = 0; i <= inputVec->size() - 1; ++i) {

    for (unsigned int index = 0; index <= vec1->size() - 1; ++index) {

      C12E2comX = (inputVec->at(i).at(vec1->at(index)).x +
                   inputVec->at(i).at((vec1->at(index)) + 1).x +
                   inputVec->at(i).at((vec1->at(index)) + 2).x +
                   inputVec->at(i).at((vec1->at(index)) + 3).x +
                   inputVec->at(i).at((vec1->at(index)) + 4).x +
                   inputVec->at(i).at((vec1->at(index)) + 5).x +
                   inputVec->at(i).at((vec1->at(index)) + 6).x) /
                  7.0;

      C12E2comY = (inputVec->at(i).at(vec1->at(index)).y +
                   inputVec->at(i).at((vec1->at(index)) + 1).x +
                   inputVec->at(i).at((vec1->at(index)) + 2).y +
                   inputVec->at(i).at((vec1->at(index)) + 3).y +
                   inputVec->at(i).at((vec1->at(index)) + 4).y +
                   inputVec->at(i).at((vec1->at(index)) + 5).y +
                   inputVec->at(i).at((vec1->at(index)) + 6).y) /
                  7.0;

      C12E2comZ = (inputVec->at(i).at(vec1->at(index)).z +
                   inputVec->at(i).at((vec1->at(index)) + 1).z +
                   inputVec->at(i).at((vec1->at(index)) + 2).z +
                   inputVec->at(i).at((vec1->at(index)) + 3).z +
                   inputVec->at(i).at((vec1->at(index)) + 4).z +
                   inputVec->at(i).at((vec1->at(index)) + 5).z +
                   inputVec->at(i).at((vec1->at(index)) + 6).z) /
                  7.0;

      C12E2input.index = vec1->at(index);
      C12E2input.x = C12E2comX;
      C12E2input.y = C12E2comY;
      C12E2input.z = C12E2comZ;
      COM1->push_back(C12E2input);
    }

    for (unsigned int index = 0; index <= vec2->size() - 1; ++index) {

      C12E2McomX = (inputVec->at(i).at(vec2->at(index)).x +
                    inputVec->at(i).at((vec2->at(index)) + 1).x +
                    inputVec->at(i).at((vec2->at(index)) + 2).x +
                    inputVec->at(i).at((vec2->at(index)) + 3).x +
                    inputVec->at(i).at((vec2->at(index)) + 4).x +
                    inputVec->at(i).at((vec2->at(index)) + 5).x +
                    inputVec->at(i).at((vec2->at(index)) + 6).x) /
                   7.0;

      C12E2McomY = (inputVec->at(i).at(vec2->at(index)).y +
                    inputVec->at(i).at((vec2->at(index)) + 1).x +
                    inputVec->at(i).at((vec2->at(index)) + 2).y +
                    inputVec->at(i).at((vec2->at(index)) + 3).y +
                    inputVec->at(i).at((vec2->at(index)) + 4).y +
                    inputVec->at(i).at((vec2->at(index)) + 5).y +
                    inputVec->at(i).at((vec2->at(index)) + 6).y) /
                   7.0;

      C12E2McomZ = (inputVec->at(i).at(vec2->at(index)).z +
                    inputVec->at(i).at((vec2->at(index)) + 1).z +
                    inputVec->at(i).at((vec2->at(index)) + 2).z +
                    inputVec->at(i).at((vec2->at(index)) + 3).z +
                    inputVec->at(i).at((vec2->at(index)) + 4).z +
                    inputVec->at(i).at((vec2->at(index)) + 5).z +
                    inputVec->at(i).at((vec2->at(index)) + 6).z) /
                   7.0;

      C12E2Minput.index = vec2->at(index);
      C12E2Minput.x = C12E2comX;
      C12E2Minput.y = C12E2comY;
      C12E2Minput.z = C12E2comZ;
      // TODO
      COM1->push_back(C12E2Minput);
    }

    C12E2final->push_back(*COM1);
    //  C12E2final->push_back(*COM2);
    // TODO
    C12E2Mfinal->push_back(*COM2);
    COM1->clear();
    COM2->clear();
  }
}

bool sortbysec_int(const std::pair<int, int> &a, const std::pair<int, int> &b) {
  return (a.second < b.second);
}

bool sortbysec_double(const std::pair<int, double> &a,
                      const std::pair<int, double> &b) {
  return (a.second < b.second);
}

class compute {
public:
  compute(){};                 // Default constructor
  compute(std::string file){}; // constructor reading file

  void storeFile() {

    boost::progress_display show_progress(number_of_snapshots);

    // open file for reading
    ipf = fopen("dump.mixed2", "r"); // Needs correction
    // check for error opening file */
    if (ipf == NULL) {
      std::cout << "Error opening file\n";
      exit(1);
    }
    // loop over the values
    for (int SSno = 0; SSno < number_of_snapshots; SSno++) {

      l = 0;
      n = 0;
      index = 0;
      // Make sure to clear all vectors before going on to the next snapshot
      inputVector.clear();
      for (int k = 0; k < nlines; k++) {
        // get a line from the file
        // fgets() returns NULL when it reaches an error or end of file

        fgets(line, sizeof(line), ipf);
        // while (fgets(line,sizeof(line),ipf) != NULL) {

        if (l < 5) {
          /* We are doing nothing */
          //	printf("%s l= %d", line,l);
          l++;
        } else if ((l > 4 && l < 8)) {
          /*We are scanning the bit with just the box parameters*/
          sscanf(line, "%lf %lf", &box1, &box2);
          //	printf("%lf %lf \n",box1,box2 );
          boxlength[l - 5] = box2 - box1;
          l++;
        } else if (l == 8) {
          // printf(" **** l= %d \n",l);
          /* We are doing nothing */
          l++;
        } else {
          //	printf(" ***** l = %d \n",l );
          /* convert the text to numbers */

          sscanf(line, "%d %d %lf %lf %lf", &index, &atomtype, &x, &y, &z);
          // std::cout << index << " " << atomtype << " " << x << " " << y << "
          // " << z; a.push_back(std::make_pair(index, index)); // Push back
          // indices b.push_back(std::make_pair(index, atomtype)); // Push back
          // atomtypes xco.push_back(std::make_pair(index, x*boxlength[0])); //
          // Push back boxlengths - x coordinates
          // yco.push_back(std::make_pair(index, y*boxlength[1])); // Push back
          // boxlengths - y coordinates zco.push_back(std::make_pair(index,
          // z*boxlength[2])); // Push back boxlengths - z coordinates
          inputline.a = index;
          inputline.b = atomtype;
          inputline.x = x * boxlength[0];
          inputline.y = y * boxlength[1];
          inputline.z = z * boxlength[2];
          inputVector.push_back(inputline);
          n++;
          l++;
          // printf("%d %d %lf %lf %lf\n",index,atomtype,x,y,z);
        }
      }
      inputTotal.push_back(inputVector);
      ++show_progress;
    }
  }

  void sortVectors() { // Sort vector values based on indices of the dump files
    for (unsigned int i = 0; i < inputTotal.size(); ++i) {
      sort(inputTotal[i].begin(), inputTotal[i].end(), compareByIndex);
    }
  }

  void printVectorElements() {
    for (unsigned int i = 0; i < inputTotal.size(); ++i) {
      for (unsigned int j = 0; j < inputTotal[1].size(); ++j) {
        // std::cout << i << " " << j << " " << " " << inputTotal[i][j].a << " "
        // << inputTotal[i][j].b << " " <<  inputTotal[i][j].x << " " <<
        // inputTotal[i][j].y << " " << inputTotal[i][j].z << std::endl;
      }
    }
  }

  /*
  void CompositionProfile() {
    //Calculating the composition profile
    for (int i = 0; i <= sizeof(C12E2_struct)/sizeof(C12E2_struct[1])-1; i++) {
  // TODO

      //double truedist(double COM1x, double COM1y, double COM1z, double COM2x,
  double COM2y, double COM2z) {
      // See if the distance between the nanopparticle is smaller than 15
      //std::cout << zco[71312] << std::endl;

      if (trueDist(xco[71312], yco[71312], zco[71312], headGroupC12_E2xCOM[i],
  headGroupC12_E2yCOM[i], headGroupC12_E2zCOM[i]) <= 25.000) {
        // Check if the c12e2 molecule is on the top layer or the bottom layer
        if (headGroupC12_E2zCOM[i] > zco[71312]) { // If on the top layer
          tophead++;
        }
        else if (headGroupC12_E2zCOM[i] < zco[71312]) { // If on the bottom
  layer downhead++;
        }
      }

      // This time, working with mimics
       if(trueDist(xco[71312], yco[71312], zco[71312], headGroupMIMCxCOM[i],
  headGroupMIMICyCOM[i], headGroupMIMICzCOM[i] ) <= 25.000) { if
  (headGroupMIMICzCOM[i] > zco[71312]) { mimictophead++;
         }
         else if (headGroupMIMICzCOM[i] < zco[71312]) {
           mimicdownhead++;
           }
         }
       }


    //std::cout << SSno << " " << tophead/(tophead + mimictophead) << " " <<
  mimictophead/(tophead + mimictophead) << " " << downhead/(downhead +
  mimicdownhead) << " " << mimicdownhead/(downhead + mimicdownhead) << " " <<
  (tophead-mimictophead)/(tophead + mimictophead) << " " <<
  (downhead-mimicdownhead)/(downhead + mimicdownhead) << " " <<
  ((tophead-mimictophead)/(tophead + mimictophead) +
  (downhead-mimicdownhead)/(downhead + mimicdownhead))/2- << " " <<
  ((tophead-mimictophead)/(tophead + mimictophead) -
  (downhead-mimicdownhead)/(downhead + mimicdownhead))/2 <<   std::endl;
    }
  */

  void check() {
    for (unsigned int i = 0; i < inputTotal.size(); ++i) {
      for (unsigned int j = 0; j <= inputTotal[1].size(); j++) {
        if (inputTotal[i][j].b == 7) {
          std::cout << inputTotal[i][j].a << " " << inputTotal[i][j].b
                    << std::endl;
          std::cout << inputTotal[i][j + 1].a << " " << inputTotal[i][j + 1].b
                    << std::endl;
          std::cout << inputTotal[i][j + 2].a << " " << inputTotal[i][j + 2].b
                    << std::endl;
          std::cout << inputTotal[i][j + 3].a << " " << inputTotal[i][j + 3].b
                    << std::endl;
          std::cout << inputTotal[i][j + 4].a << " " << inputTotal[i][j + 4].b
                    << std::endl;
          std::cout << inputTotal[i][j + 5].a << " " << inputTotal[i][j + 5].b
                    << std::endl;
          std::cout << inputTotal[i][j + 6].a << " " << inputTotal[i][j + 6].b
                    << std::endl;
        } else if (inputTotal[i][j].b == 13) {
          std::cout << inputTotal[i][j].a << " " << inputTotal[i][j].b
                    << std::endl;
          std::cout << inputTotal[i][j + 1].a << " " << inputTotal[i][j + 1].b
                    << std::endl;
          std::cout << inputTotal[i][j + 2].a << " " << inputTotal[i][j + 2].b
                    << std::endl;
          std::cout << inputTotal[i][j + 3].a << " " << inputTotal[i][j + 3].b
                    << std::endl;
          std::cout << inputTotal[i][j + 4].a << " " << inputTotal[i][j + 4].b
                    << std::endl;
          std::cout << inputTotal[i][j + 5].a << " " << inputTotal[i][j + 5].b
                    << std::endl;
          std::cout << inputTotal[i][j + 6].a << " " << inputTotal[i][j + 6].b
                    << std::endl;
        }
      }
    }
  }

  void headGroupVectorFormation() {
    for (unsigned int j = 0; j <= inputTotal[1].size(); j++) {
      if (inputTotal[1][j].b == 7) {
        C12E2IndexVector.push_back(
            inputTotal[1][j].a); // push back C12E2 bead 7 indices (headgroups)
      } else if (inputTotal[1][j].b == 13) {
        C12E2MIndexVector.push_back(
            inputTotal[1][j].a); // push back C12E2 bead 7 indices (headgroups)
      }
    }

    for (unsigned int index = 0; index <= C12E2IndexVector.size(); ++index) {
      std::cout << C12E2IndexVector[index] << std::endl;
    }

    for (unsigned int index = 0; index <= C12E2MIndexVector.size(); ++index) {
      std::cout << C12E2MIndexVector[index] << std::endl;
    }
  }
  
  void ComputePhi() {
    
    // std::cout << C12E2IndexVector.size() << " " <<  C12E2MIndexVector.size()
    // << " " << inputTotal.size() << std::endl;

    /*
      Following Python's comment system """ """
      
      Computes the phi, or the mismatch between the bilayer leaflets around the NP

      We are explcitily ignoring flip-flops in the case of this study
    */

    for (unsigned int index = 0; index < C12E2IndexVector.size(); ++index) {
      if (inputTotal[1][C12E2IndexVector[index]].z > 120.0) {

	topC12E2Index.push_back(C12E2IndexVector[index]);
      } else if (inputTotal[1][C12E2IndexVector[index]].z < 120.0) {
        botC12E2Index.push_back(C12E2IndexVector[index]);
      }
    }

    for (unsigned int index = 0; index < C12E2MIndexVector.size(); ++index) {
      if (inputTotal[1][C12E2MIndexVector[index]].z > 120.0) {
        topC12E2MIndex.push_back(C12E2MIndexVector[index]);
      } else if (inputTotal[1][C12E2MIndexVector[index]].z < 120.0) {
        botC12E2MIndex.push_back(C12E2MIndexVector[index]);
      }
    }

    std::cout << "check" << std::endl;

    for (std::vector<int>::iterator it = topC12E2Index.begin();
         it != topC12E2Index.end(); it++) {
      // std::cout << *it << " " << std::endl;
      //	std::cout << topC12E2Index.size() << std::endl;
    }

    for (std::vector<int>::iterator it = botC12E2Index.begin();
         it != botC12E2Index.end(); it++) {
      // std::cout << *it << " " << std::endl;
      //	std::cout << botC12E2Index.size() << std::endl;
    }

    std::cout << topC12E2Index.size() << std::endl;
    std::cout << botC12E2Index.size() << std::endl;

    std::cout << topC12E2MIndex.size() << std::endl;
    std::cout << botC12E2MIndex.size() << std::endl;
    phiStruct phitemplate;
    std::cout << "check" << std::endl;

    for (unsigned int i = 0; i <= inputTotal.size() - 1; ++i) {

      for (unsigned phiIndex = 0; phiIndex <= 100; ++phiIndex) {
        phitemplate.dist = phiIndex;
        phitemplate.topPhiC12E2Count = 0;
        phitemplate.botPhiC12E2Count = 0;
        phitemplate.topPhiC12E2MCount = 0;
        phitemplate.botPhiC12E2MCount = 0;
        phiCount.push_back(phitemplate);
      }

      NPX = inputTotal[i][71312].x; // x coordinate of the NP
      NPY = inputTotal[i][71312].y; // y coordinate of the NP
      NPZ = inputTotal[i][71312].z; // z coordinate of the NP

      int topC12E2Counter, botC12E2Counter, topC12E2MCounter, botC12E2MCounter;

      for (unsigned int dist = 0; dist <= 100; ++dist) {

        for (unsigned int index = 0; index < topC12E2Index.size(); index++) {
          DistVec1 =
              trueDist(&NPX, &NPY, &NPZ, &inputTotal[i][topC12E2Index[index]].x,
                       &inputTotal[i][topC12E2Index[index]].y,
                       &inputTotal[i][topC12E2Index[index]].z);
          if (DistVec1 <= dist + 0.5 && DistVec1 >= dist - 0.5) {
            phiCount[dist].topPhiC12E2Count += 1;
          }
        }

        for (unsigned int index = 0; index < botC12E2Index.size(); index++) {
          DistVec2 =
              trueDist(&NPX, &NPY, &NPZ, &inputTotal[i][botC12E2Index[index]].x,
                       &inputTotal[i][botC12E2Index[index]].y,
                       &inputTotal[i][botC12E2Index[index]].z);
          if (DistVec2 <= dist + 0.5 && DistVec2 >= dist - 0.5) {
            phiCount[dist].botPhiC12E2Count += 1;
          }
        }

        for (unsigned int index = 0; index < topC12E2MIndex.size(); index++) {
          DistVec3 = trueDist(&NPX, &NPY, &NPZ,
                              &inputTotal[i][topC12E2MIndex[index]].x,
                              &inputTotal[i][topC12E2MIndex[index]].y,
                              &inputTotal[i][topC12E2MIndex[index]].z);
          if (DistVec3 <= dist + 0.5 && DistVec3 >= dist - 0.5) {
            phiCount[dist].topPhiC12E2MCount += 1;
          }
        }

        for (unsigned int index = 0; index < botC12E2MIndex.size(); index++) {
          DistVec4 = trueDist(&NPX, &NPY, &NPZ,
                              &inputTotal[i][botC12E2MIndex[index]].x,
                              &inputTotal[i][botC12E2MIndex[index]].y,
                              &inputTotal[i][botC12E2MIndex[index]].z);
          if (DistVec4 <= dist + 0.5 && DistVec4 >= dist - 0.5) {
            phiCount[dist].botPhiC12E2MCount += 1;
          }
        }
      }
      
      phiTotal.push_back(phiCount);
      phiCount.clear();
    }
  }

  void PhiPrint() {
    phipm A;
    for (unsigned phiIndex = 0; phiIndex <= 100; ++phiIndex) {
      A.phim = 0.0;
      A.phip = 0.0;
      A.phimVec.clear();
      A.phipVec.clear();
      NewNew.push_back(A);
    }

    double phi1, phi2;
    double phip, phim;

    for (unsigned int index = 0; index <= phiTotal.size() - 1; index++) {
      for (unsigned int index2 = 0; index2 <= phiTotal[0].size(); index2++) {

        phi1 = 0.0;
        phi2 = 0.0;
        phip = 0.0;
        phim = 0.0;

        // std::cout << index << " " << index2 << " " <<
        // phiTotal[index][index2].topPhiC12E2Count << " " <<
        // phiTotal[index][index2].botPhiC12E2Count << " " <<
        // phiTotal[index][index2].topPhiC12E2MCount << " " <<
        // phiTotal[index][index2].botPhiC12E2MCount << std::endl;
        if ((double)(phiTotal[index][index2].topPhiC12E2Count +
                     phiTotal[index][index2].topPhiC12E2MCount) == 0.0 ||
            (double)(phiTotal[index][index2].botPhiC12E2Count +
                     phiTotal[index][index2].botPhiC12E2MCount) == 0) {
          // Do nothing

        } else {
          phi1 = (double)(phiTotal[index][index2].topPhiC12E2Count -
                          phiTotal[index][index2].topPhiC12E2MCount) /
                 (double)(phiTotal[index][index2].topPhiC12E2Count +
                          phiTotal[index][index2].topPhiC12E2MCount);
          phi2 = (double)(phiTotal[index][index2].botPhiC12E2Count -
                          phiTotal[index][index2].botPhiC12E2MCount) /
                 (double)(phiTotal[index][index2].botPhiC12E2Count +
                          phiTotal[index][index2].botPhiC12E2MCount);
        }

        phip = phi2 + phi1 / 2.0;
        phim = phi2 - phi1 / 2.0;

        NewNew[index2].phim += phim;
        NewNew[index2].phip += phip;
        NewNew[index2].phimVec.push_back(phim);
        NewNew[index2].phipVec.push_back(phip);

        std::cout << index << " " << index2 << " " << phip << " " << phim
                  << std::endl;
      }
    }
  }

  void ComputePhiStandardDev() {
    double sum, mean, sq_sum, stdev;
    for (std::vector<phipm>::iterator it = NewNew.begin(); it != NewNew.end();
         it++) {
      // std::cout << " "  << it - NewNew.begin() << " " <<
      // (it->phim)/inputTotal.size() << " " << (it->phip)/inputTotal.size() <<
      // " " <<  std::endl; std::cout << (it->phimVec).size() << " " <<
      // (it->phipVec).size() << std::endl;
      for (unsigned int i = 0; i != NewNew[it - NewNew.begin()].phimVec.size();
           ++i) {

        sum = std::accumulate(NewNew[it - NewNew.begin()].phimVec.begin(),
                              NewNew[it - NewNew.begin()].phimVec.end(),
                              0.0); // Compute sum

        mean = sum / NewNew[it - NewNew.begin()].phimVec.size(); // Compute Mean

        sq_sum = std::inner_product(NewNew[it - NewNew.begin()].phimVec.begin(),
                                    NewNew[it - NewNew.begin()].phimVec.end(),
                                    NewNew[it - NewNew.begin()].phimVec.begin(),
                                    0.0); // Compute square sum

        stdev = std::sqrt(sq_sum / NewNew[it - NewNew.begin()].phimVec.size() -
                          mean * mean) /
                (pow(number_of_snapshots, 0.5));
      }

      std::cout << " " << it - NewNew.begin() << " "
                << (it->phip) / inputTotal.size() << " " << stdev << std::endl;
    }
  }

  void ComputeOrderphobic() { // Computes the phi, or the mismatch between the
                              // bilayer leaflets around the NP

    double DIST, DIST2;   // TODO
    double DISTM, DISTM2; // TODO

    for (unsigned int i = 0; i < inputTotal.size(); ++i) {

      for (unsigned int index = 0; index < C12E2IndexVector.size(); ++index) {

        for (unsigned int newindex = 0; newindex < C12E2IndexVector.size();
             ++newindex) {

          DIST = trueDist(&inputTotal[i][C12E2IndexVector[index] + 4].x,
                          &inputTotal[i][C12E2IndexVector[index] + 4].y,
                          &inputTotal[i][C12E2IndexVector[index] + 4].z,
                          &inputTotal[i][C12E2IndexVector[newindex] + 4].x,
                          &inputTotal[i][C12E2IndexVector[newindex] + 4].y,
                          &inputTotal[i][C12E2IndexVector[newindex] + 4].z);

	  C12E2sample.index = C12E2IndexVector[index];
          C12E2sample.dist = DIST;
          C12E2sample.selfXcoord = inputTotal[i][C12E2IndexVector[index] + 4].x;
          C12E2sample.selfYcoord = inputTotal[i][C12E2IndexVector[index] + 4].y;
          C12E2sample.selfZcoord = inputTotal[i][C12E2IndexVector[index] + 4].z;
          C12E2sample.Xcoord = inputTotal[i][C12E2IndexVector[newindex] + 4].x;
          C12E2sample.Ycoord = inputTotal[i][C12E2IndexVector[newindex] + 4].y;
          C12E2sample.Ycoord = inputTotal[i][C12E2IndexVector[newindex] + 4].z;
          C12E2orderphobic.push_back(C12E2sample);
        }

        for (unsigned int newindex = 0; newindex < C12E2IndexVector.size();
             ++newindex) {

          DISTM = trueDist(&inputTotal[i][C12E2IndexVector[index] + 4].x,
                           &inputTotal[i][C12E2IndexVector[index] + 4].y,
                           &inputTotal[i][C12E2IndexVector[index] + 4].z,
                           &inputTotal[i][C12E2MIndexVector[newindex] + 4].x,
                           &inputTotal[i][C12E2MIndexVector[newindex] + 4].y,
                           &inputTotal[i][C12E2MIndexVector[newindex] + 4].z);

          C12E2sample.index = C12E2IndexVector[index];
          C12E2sample.dist = DISTM;
          C12E2sample.selfXcoord = inputTotal[i][C12E2IndexVector[index] + 4].x;
          C12E2sample.selfYcoord = inputTotal[i][C12E2IndexVector[index] + 4].y;
          C12E2sample.selfZcoord = inputTotal[i][C12E2IndexVector[index] + 4].z;
          C12E2sample.Xcoord = inputTotal[i][C12E2MIndexVector[newindex] + 4].x;
          C12E2sample.Ycoord = inputTotal[i][C12E2MIndexVector[newindex] + 4].y;
          C12E2sample.Ycoord = inputTotal[i][C12E2MIndexVector[newindex] + 4].z;
          C12E2orderphobic.push_back(C12E2sample);
        }

        sort(C12E2orderphobic.begin(), C12E2orderphobic.end(),
             compareByIndexOPh);
        C12E2orderphobic.erase(C12E2orderphobic.begin());
        C12E2orderphobic.erase(C12E2orderphobic.begin() + 6,
                               C12E2orderphobic.end());
        //	C12E2orderphobic.resize(6);
        C12E2orderphobicVec.push_back(C12E2orderphobic);

        for (unsigned int newindex = 0; newindex < C12E2MIndexVector.size();
             ++newindex) {
          DISTM2 = trueDist(&inputTotal[i][C12E2MIndexVector[index] + 4].x,
                            &inputTotal[i][C12E2MIndexVector[index] + 4].y,
                            &inputTotal[i][C12E2MIndexVector[index] + 4].z,
                            &inputTotal[i][C12E2MIndexVector[newindex] + 4].x,
                            &inputTotal[i][C12E2MIndexVector[newindex] + 4].y,
                            &inputTotal[i][C12E2MIndexVector[newindex] + 4].z);

          C12E2Msample.index = C12E2MIndexVector[index];
          C12E2Msample.dist = DISTM2;

          C12E2Msample.selfXcoord =
              inputTotal[i][C12E2MIndexVector[index] + 4].x;
          C12E2Msample.selfYcoord =
              inputTotal[i][C12E2MIndexVector[index] + 4].y;
          C12E2Msample.selfZcoord =
              inputTotal[i][C12E2MIndexVector[index] + 4].z;

          C12E2Msample.Xcoord =
              inputTotal[i][C12E2MIndexVector[newindex] + 4].x;
          C12E2Msample.Ycoord =
              inputTotal[i][C12E2MIndexVector[newindex] + 4].y;
          C12E2Msample.Ycoord =
              inputTotal[i][C12E2MIndexVector[newindex] + 4].z;

          C12E2Morderphobic.push_back(C12E2Msample);
        }

        for (unsigned int newindex = 0; newindex < C12E2MIndexVector.size();
             ++newindex) {
          DIST2 = trueDist(&inputTotal[i][C12E2MIndexVector[index] + 4].x,
                           &inputTotal[i][C12E2MIndexVector[index] + 4].y,
                           &inputTotal[i][C12E2MIndexVector[index] + 4].z,
                           &inputTotal[i][C12E2IndexVector[newindex] + 4].x,
                           &inputTotal[i][C12E2IndexVector[newindex] + 4].y,
                           &inputTotal[i][C12E2IndexVector[newindex] + 4].z);

          C12E2Msample.index = C12E2MIndexVector[index];

          C12E2Msample.dist = DIST2;

          C12E2Msample.selfXcoord =
              inputTotal[i][C12E2MIndexVector[index] + 4].x;
          C12E2Msample.selfYcoord =
              inputTotal[i][C12E2MIndexVector[index] + 4].y;
          C12E2Msample.selfZcoord =
              inputTotal[i][C12E2MIndexVector[index] + 4].z;

          C12E2Msample.Xcoord = inputTotal[i][C12E2IndexVector[newindex] + 4].x;
          C12E2Msample.Ycoord = inputTotal[i][C12E2IndexVector[newindex] + 4].y;
          C12E2Msample.Ycoord = inputTotal[i][C12E2IndexVector[newindex] + 4].z;

          C12E2Morderphobic.push_back(C12E2Msample);
        }

        // std::sort(c12E2Morderphobic.begin(), c12E2Morderphobic.end(),
        // compareByIndexOPh); std::sort(C12E2Morderphobic.begin(),
        // C12E2Morderphobic.end(), std::greater<CustomStruct>());
        // C12E2Morderphobic.erase(C12E2Morderphobic.begin());
        // C12E2Morderphobic.resize(6);

        sort(C12E2Morderphobic.begin(), C12E2Morderphobic.end(),
             compareByIndexOPh);
        C12E2Morderphobic.erase(C12E2Morderphobic.begin());
        C12E2Morderphobic.erase(C12E2Morderphobic.begin() + 6,
                                C12E2Morderphobic.end());

        // C12E2Morderphobic.resize(6);

        C12E2MorderphobicVec.push_back(C12E2Morderphobic);

        C12E2orderphobic.clear();
        C12E2Morderphobic.clear();
      }

      // Final push_back
      orderphobicC12E2.push_back(C12E2orderphobicVec);
      orderphobicC12E2M.push_back(C12E2MorderphobicVec);

      C12E2MorderphobicVec.clear();
      C12E2orderphobicVec.clear();
    }
  }

  void OrderphobicSort() {
    // Computes the phi, or the mismatch between the
    // bilayer leaflets around the NP

    OPVal Val;
    std::vector<OPVal> tempVal;
    double output;

    for (unsigned int i = 0; i < orderphobicC12E2.size() - 1; ++i) {
      for (unsigned int index = 0; index < C12E2IndexVector.size(); ++index) {
        output = calcAngle(&orderphobicC12E2[i][index]);

        std::cout << output << " " << i << " " << index << " "
                  << orderphobicC12E2[i][index][0].selfXcoord << " "
                  << orderphobicC12E2[i][index][0].selfYcoord << " "
                  << orderphobicC12E2[i][index][0].selfZcoord << std::endl;

        Val.X = orderphobicC12E2[i][index][0].selfXcoord;
        Val.Y = orderphobicC12E2[i][index][0].selfYcoord;
        Val.Z = orderphobicC12E2[i][index][0].selfZcoord;
        Val.orderphobicVal = output;
        tempVal.push_back(Val);
      }

      for (unsigned int index = 0; index < C12E2MIndexVector.size(); ++index) {
        output = calcAngle(&orderphobicC12E2M[i][index]);
        std::cout << output << " " << i << " " << index << " "
                  << orderphobicC12E2M[i][index][0].selfXcoord << " "
                  << orderphobicC12E2M[i][index][0].selfYcoord << " "
                  << orderphobicC12E2M[i][index][0].selfZcoord << std::endl;
        Val.X = orderphobicC12E2M[i][index][0].selfXcoord;
        Val.Y = orderphobicC12E2M[i][index][0].selfYcoord;
        Val.Z = orderphobicC12E2M[i][index][0].selfZcoord;
        Val.orderphobicVal = output;
        tempVal.push_back(Val);
      }

      orderphobicVectorFinal.push_back(tempVal);
      tempVal.clear();
    }
  }

  void printop() {

    OPVal PlaceHolder;
    double distP;
    for (unsigned int dist = 0; dist <= 100; ++dist) {
      PlaceHolder.X = 0.0;
      PlaceHolder.Y = 0.0;
      PlaceHolder.Z = 0.0;
      ABC.push_back(PlaceHolder);
    }

    for (unsigned int i = 0; i != orderphobicVectorFinal.size(); i++) {
      for (unsigned int j = 0; j != orderphobicVectorFinal[i].size(); j++) {

        NPX = inputTotal[i][71312].x; // x coordinate of the NP
        NPY = inputTotal[i][71312].y; // y coordinate of the NP
        NPZ = inputTotal[i][71312].z; // z coordinate of the NP

        // std::cout <<  i << " " << j << " " << orderphobicVectorFinal[i][j].X
        // << " " << orderphobicVectorFinal[i][j].Y << " " <<
        // orderphobicVectorFinal[i][j].Z << " " <<
        // orderphobicVectorFinal[i][j].orderphobicVal << " "  << NPX << " " <<
        // NPY << " " << NPZ << std::endl;

        distP = trueDist(&NPX, &NPY, &NPZ, &orderphobicVectorFinal[i][j].X,
                         &orderphobicVectorFinal[i][j].Y,
                         &orderphobicVectorFinal[i][j].Z);

        for (unsigned int dist = 0; dist <= 100; ++dist) {
          if (distP >= dist - 0.5 && distP <= dist + 0.5) {
            ABC[dist].VV.push_back(orderphobicVectorFinal[i][j].orderphobicVal);
          }
        }
      }
    }
  }

  void LargePrint() {
    double sum, mean;
    for (unsigned int i = 0; i != ABC.size(); i++) {
      sum = std::accumulate(ABC[i].VV.begin(), ABC[i].VV.end(),
                            0.0);    // Compute sum
      mean = sum / ABC[i].VV.size(); // Compute Mean
      std::cout << i << " " << mean << " " << std::endl;
    }
  }

  void allocateCOM() {
    /*
     */
    CenterOfMass(&C12E2IndexVector, &C12E2MIndexVector, &inputTotal, &C12E2COM,
                 &C12E2MCOM, &C12E2TotalCOMArray, &C12E2MTotalCOMArray);
  }

private:
  std::vector<OPHstruct> inputVectorStruct; // push back all structs
  std::vector<std::vector<std::vector<OPHstruct>>>
      orderphobicC12E2; // This needs to be a simplified
  std::vector<std::vector<std::vector<OPHstruct>>>
      orderphobicC12E2M; // This needs to be simplfied also..
  // Vectors to store trajectory values
  std::vector<inputCoord> inputVector; // push back all structs
  std::vector<std::vector<inputCoord>>
      inputTotal; // push back vector of structs
  std::vector<int>
      C12E2IndexVector; // push back C12E2 bead 7 indices (headgroups)
  std::vector<int>
      C12E2MIndexVector; // push back all C12E2M bead 13 indices (headgroups)
  // Vectors to store the Centers of Mass
  std::vector<COMstruct> C12E2COM;
  std::vector<COMstruct> C12E2MCOM;
  std::vector<phiStruct> phiCount;
  std::vector<std::vector<phiStruct>> phiTotal;
  std::vector<std::vector<COMstruct>> C12E2TotalCOMArray;
  std::vector<std::vector<COMstruct>> C12E2MTotalCOMArray;
  // top head groups C12E2
  std::vector<int> topC12E2Index;
  // top head groups C12E2M
  std::vector<int> botC12E2Index;
  // bot head groups C12E2M
  std::vector<int> topC12E2MIndex;
  // bot head groups C12E2M
  std::vector<int> botC12E2MIndex;
  // TODO
  std::vector<phipm> NewNew; // TODO - need to rename

  FILE *ipf; /* input file */
  int atomtype;
  int index, l, n;
  int nlines = number_of_atoms + 9;
  double x, y, z; /*coordinates for the atoms in the box*/
  double box1;
  double box2;
  double NPX, NPY, NPZ;
  double boxlength[box_dim];
  char line[100];
  inputCoord inputline;
  double DistVec1, DistVec2, DistVec3, DistVec4;

  OPHstruct C12E2sample;
  OPHstruct C12E2Msample;

  std::vector<std::vector<OPHstruct>> C12E2orderphobicVec;  // TODO
  std::vector<std::vector<OPHstruct>> C12E2MorderphobicVec; // TODO
  std::vector<OPHstruct> C12E2orderphobic;
  std::vector<OPHstruct> C12E2Morderphobic;
  std::vector<std::vector<OPVal>> orderphobicVectorFinal; // TODO
  std::vector<OPVal> ABC;
};

class testClass {};

// Intitiate test class

compute C12E2PhiOrderphobic;

int main(int argc, char *argv[]) {
  C12E2PhiOrderphobic.storeFile(); // 
  C12E2PhiOrderphobic.sortVectors(); // 
  C12E2PhiOrderphobic.check(); // 
  C12E2PhiOrderphobic.headGroupVectorFormation(); // 
  C12E2PhiOrderphobic.allocateCOM(); // 
  C12E2PhiOrderphobic.ComputePhi(); // 
  C12E2PhiOrderphobic.PhiPrint(); // 
  C12E2PhiOrderphobic.ComputePhiStandardDev(); // 
  C12E2PhiOrderphobic.ComputeOrderphobic(); //
  C12E2PhiOrderphobic.OrderphobicSort(); // 
  C12E2PhiOrderphobic.printop(); // 
  C12E2PhiOrderphobic.LargePrint(); //

  return 0;
}
