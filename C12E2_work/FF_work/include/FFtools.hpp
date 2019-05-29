#ifndef _FFtools_hpp_
#define _FFtools_hpp_

class FFtools {
public:
  FFtools(); // default contructor
  FFtools(const FFtools& ref);  // Copy constructor
  FFtools& operator=(const FFtools ref); // Assignment operator overloaded
  virtual ~FFtools(); // Destructor
  
  double LJ(double dist_ab, double eps_ab, double sigma_ab);  
  double distance(double ax,double ay,double az, double bx, double by, double bz); 
  double minDist(double ax,double ay,double az, double bx, double by, double bz, double bxx, double byy, double bzz);
  void minimumImage(double *a,double *b,double *c,double bxx,double byy,double bzz); // implementing the minimum image conventon
  double centreOfMass(int Oa, int H1a, int H2a,int Ob, int H1b, int H2b,int Oc, int H1c, int H2c);
  void openDumpFile();
  
private:  
  /* We are now declaring the atom coordinates of each individual atom in the water triad */ 
  double Obcoordinatex; /*cartesian coordinates for the oxygen on water b */ 
  double Obcoordinatey; /*cartesian coordinates for the oxygen on water b */ 
  double Obcoordinatez; /*cartesian coordinates for the oxygen on water b */ 
  
  double Occoordinatex; /*cartesian coordinates for the oxygen on water c */ 
  double Occoordinatey; /*cartesian coordinates for the oxygen on water c */
  double Occoordinatez; /*cartesian coordinates for the oxygen on water c */
  
  double Odcoordinatex; /*cartersian coordinates for the oxygen on water d */
  double Odcoordinatey; /*cartersian coordinates for the oxygen on water d */
  double Odcoordinatez; /*cartersian coordinates for the oxygen on water d */

  /* LJ parameters */
  double Oepsilon = 0.15207;
  double Osigma = 3.1507;
  double Hepsilon = 0.000;
  double Hsigma = 0.0;
  double goldSigma = 2.629; 
  double goldEpsilon = 5.29;

  const int numberOfAtoms = 3271; 
  const int boxdim = 3; 
  const int nbins = 100;
};


#endif
