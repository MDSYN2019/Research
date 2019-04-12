#include <iostream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
#include <utility>
#include <cerrno>


#define minClust 5 // the minimum sizes cluster we want to count as a cluster

double CenterOfMass(double H7, double H6_1,double H6_2,double T3_1, double T3_2, double T3_3, double T4);
double truedist(double COM1x, double COM1y, double COM1z, double COM2x, double COM2y, double COM2z);

void die (const char *message) {
  if (errno) {
    perror(message);
  } else {
    std::cout << message << std::endl;
  }
  exit(1);
}

struct pointCluster {
  std::vector<int> clusterIndices;
};

pointCluster* regionQuery (int , double*, double*, double*, double*, double*, double*);


// The struct to construct the array 

struct C12E2_skeleton {
  int index[7];
};

// struct to add up the indices of the cluster around it, if the number of the surrounding is

std::vector<int> C12E2I;
std::vector<int> C12E2MI;
std::map<int, C12E2_skeleton> C12E2_M;
std::map<int, C12E2_skeleton> C12E2M_M;

// The regionQuery returns the number of points that are within the criteria of being a cluster - i.e. a distance costraint and a type
// constaint for the lipid. We already have a type constraint in terms of the struct arrays for the C12E2 and the C12E2-M structs so we
// dont need to worry about that, but we need to worry about the region and indices we want to return

int MimicCounter = 0;    
int PolymerCounter = 0;



/* -------------- 
   
   Calculate the clustering of the mixed lipid bilayer through a DBSCAN algorithm --- How does the DBSCAN algorithm work?
   
   -> Point P in a cluster is a 'core' point if there are a critical number of the same type of points within a distance E.
   -> Point Q is a point in the cluster if there is a clear vector towards that from any of the elements of the P vector 
   -> All points not reachable from any other points are outliers
   
   Two main algorithms are required, where there are 
   
   -------------- */

/* Function to find the center of mass of the C12E2 */

// --- notification of end of a function --- //





//double regionQuery (int P, int D);



// DBSCAN algorithm -  Pseudocode:
// -------------------------------//

/*
DBSCAN(D, eps, MinPts) {
C = 0
for each point P in dataset D {
  if P is visited 
     continue next point 
  mark P as visited
  NeighbourPts = regionQuery(P, eps)
  if Sizeof(NeighbourPts) < MinPts
     mark P as NOISE
  else {
     C = next cluster
     expandCluster(P, NeighbourPts, C, eps, MinPits)
  }
 }
}

 */

/*
void DBSCAN(double *, eps, Minpts);
void regionQuery(double Distcriteria, double (*COMref) (double)) {
  COMref = CenterOfMass;  
}  
*/
const int numberofatoms = 71313;  
int boxdim = 3;
int atomtype; 
double x,y,z; /*coordinates for the atoms in the box*/
double xco[numberofatoms],yco[numberofatoms],zco[numberofatoms];


int main () {  
  /* Declaration of variables */ 
  
  
  /* x,y,z array to store all the values we read in from the dump file in order */
  
 
  /* ditto for a and b, which represent index and atomtype respectively */

  int a[numberofatoms],b[numberofatoms];
 
  /*box dimensions*/

  double box1;
  double box2; 
  double boxlength[boxdim];
 
  /* Parameters to loop over an entire dump file */ 

  int numberofSS = 1000; /*The number of screenshots in the dump file*/
  int SSno = 0; /*the nth screenshot we are at*/

  /* MISC */

  char line[100];
  int n = 0; 
  int index;
  int i; 
  int j;
  int k = 0;
  int l = 0;

  // -------------------------------------------------//
  /* Defining C12E2 atoms (A7 -- A6_1 -- A6_2 -- A3_1 -- A3_2 -- A3_3 -- A4*/ 
  // -------------------------------------------------//

  // First Batch

  int A7; 
  int A6_1;
  int A6_2;
  int A3_1;
  int A3_2; 
  int A3_3; 
  int A4;
 
  // Second Batch
  
  int A7_2; 
  int A6_1_2;
  int A6_2_2;
  int A3_1_2;
  int A3_2_2; 
  int A3_3_2; 
  int A4_2;
  
  // Third Batch
  
  int A7_3; 
  int A6_1_3;
  int A6_2_3;
  int A3_1_3;
  int A3_2_3; 
  int A3_3_3; 
  int A4_3;
  
 // Fourth Batch
 
  int A7_4; 
  int A6_1_4;
  int A6_2_4;
  int A3_1_4;
  int A3_2_4; 
  int A3_3_4; 
  int A4_4;

  // -------------------------------------------------//

  
  // -------------------------------------------------//
  /* Defining mimics (A13 -- A12_1 -- A12_2 -- A9_1 -- A9_2 -- A9_3 -- A10) */
  // -------------------------------------------------//

  // First Batch

  int A13; 
  int A12_1;
  int A12_2;
  int A9_1;
  int A9_2; 
  int A9_3; 
  int A10;

  // Second Batch

  int A13_2; 
  int A12_1_2;
  int A12_2_2;
  int A9_1_2;
  int A9_2_2; 
  int A9_3_2; 
  int A10_2;

  // Third Batch

  int A13_3; 
  int A12_1_3;
  int A12_2_3;
  int A9_1_3;
  int A9_2_3; 
  int A9_3_3; 
  int A10_3;

  // Fourth Batch

  int A13_4; 
  int A12_1_4;
  int A12_2_4;
  int A9_1_4;
  int A9_2_4; 
  int A9_3_4; 
  int A10_4;

  // -------------------------------------------------//

  /* Array parameters to store in the COM of each molecule */ 
  /* Array to store in C12E2 only */ 

  /* This will have to change */
  
  int numberOfPolymers = 1000; // The number of polymers of each type - C12E2 or mimic 
  int indexCG = 7;
    
  /* Array for the C12E2 groups */

  double headGroupC12_E2xCOM[numberOfPolymers];
  double headGroupC12_E2yCOM[numberOfPolymers];
  double headGroupC12_E2zCOM[numberOfPolymers];
  
  /* Array to store in mimics only! */  

  double headGroupMIMICxCOM[numberOfPolymers];
  double headGroupMIMICyCOM[numberOfPolymers];
  double headGroupMIMICzCOM[numberOfPolymers];

  // int ww[numberofwaters];

  /* Start of code */ 
  /* Reading in the dump file */
  /* pointer for the dump file */

  FILE *ipf; /* input file */
  
  /* open file for reading */
  ipf = fopen("dump.out", "r");
  
  /* check for error opening file */
    
  if(ipf == NULL) {  
    printf("Error opening file\n");
    exit(1);
  }
  
  /* get a line from the file */
  /* fgets() returns NULL when it reaches an error or end of file */   
  
  int nlines = numberofatoms + 9;      
 
  // Ensure that the vectors are cleared first.

  C12E2I.clear();
  C12E2MI.clear();

  // Structs to keep the index of each polymer and it's individual CG bead indices
  
  struct C12E2_skeleton C12E2_struct[numberOfPolymers];
  struct C12E2_skeleton C12E2M_struct[numberOfPolymers];

  // Structs to keep the indices of the CG beads 
  
  for (int i = 0; i <= 249; i++) {
    
    // Allocate the first batch of indices    
    // --- C12E2 --- //
    // First batch
    
    A7 = 7*(i);
    A6_1= 7*(i) + 1;
    A6_2= 7*(i) + 2;
    A3_1 = 7*(i) + 3;
    A3_2 = 7*(i) + 4;
    A3_3 = 7*(i) + 5;
    A4 = 7*(i) + 6;
    
    // Second batch
    
    A7_2 = 7*(i) + 3501;
    A6_1_2 = 7*(i) + 3502;
    A6_2_2 = 7*(i) + 3503;
    A3_1_2 = 7*(i) + 3504;
    A3_2_2 = 7*(i) + 3505;
    A3_3_2 = 7*(i) + 3506;
    A4_2 = 7*(i) + 3507;
    
    // Third batch
      
    A7_3 = 7*(i) + 7001;
    A6_1_3 = 7*(i) + 7002;
    A6_2_3 = 7*(i) + 7003;
    A3_1_3 = 7*(i) + 7004;
    A3_2_3 = 7*(i) + 7005;
    A3_3_3 = 7*(i) + 7006;
    A4_3 = 7*(i) + 7007;
    
    // Fourth batch
    
    A7_4 = 7*(i)  + 10501;
    A6_1_4 = 7*(i) + 10502;
    A6_2_4 = 7*(i) + 10503;
    A3_1_4 = 7*(i) + 10504;
    A3_2_4 = 7*(i) + 10505;
    A3_3_4 = 7*(i) + 10506;
    A4_4 = 7*(i) + 10507;
    
    //  --- C12E2-M --- //
    
    
    A13 = 7*(i) + 1751;
    A12_1= 7*(i)+ 1752;
    A12_2= 7*(i)+ 1753;
    A9_1 = 7*(i)+ 1754;
    A9_2 = 7*(i)+ 1755;
    A9_3 = 7*(i)+ 1756;
    A10 = 7*(i)+ 1757;
    
    /* Mimic atoms - Second batch */
      
    A13_2 = 7*(i) + 5251;
    A12_1_2 = 7*(i)+ 5252;
    A12_2_2 = 7*(i)+ 5253;
    A9_1_2 = 7*(i)+ 5254;
    A9_2_2 = 7*(i)+ 5255;
    A9_3_2 = 7*(i)+ 5256;
    A10_2 = 7*(i)+ 5257;
    
    /* Mimic atoms - Third batch */
    
    A13_3 = 7*(i) + 8751;
    A12_1_3 = 7*(i)+ 8752;
    A12_2_3 = 7*(i)+ 8753;
    A9_1_3 = 7*(i)+ 8754;
    A9_2_3 = 7*(i)+ 8755;
    A9_3_3 = 7*(i)+ 8756;
    A10_3 = 7*(i)+ 8757;
    
    /* Mimic atoms - Fourth batch */
  
    A13_4 = 7*(i) + 12251;
    A12_1_4 = 7*(i)+ 12252;
    A12_2_4 = 7*(i)+ 12253;
    A9_1_4 = 7*(i)+ 12254;
    A9_2_4 = 7*(i)+ 12255;
    A9_3_4 = 7*(i)+ 12256;
    A10_4 = 7*(i)+ 12257;
        
    // ----------- C12E2 indicies ----------- //

    //std::map<int, int[7]> C12E2_M;
    //std::map<int, int[7]> C12E2M_M;


    // C12E2 assignment

    C12E2_struct[i].index[0] = A7;
    C12E2_struct[i].index[1] = A6_1;
    C12E2_struct[i].index[2] = A6_2;
    C12E2_struct[i].index[3] = A3_1;
    C12E2_struct[i].index[4] = A3_2;
    C12E2_struct[i].index[5] = A3_3;
    C12E2_struct[i].index[6] = A4;

    C12E2_struct[i+250].index[0] = A7_2;
    C12E2_struct[i+250].index[1] = A6_1_2;
    C12E2_struct[i+250].index[2] = A6_2_2;
    C12E2_struct[i+250].index[3] = A3_1_2;
    C12E2_struct[i+250].index[4] = A3_2_2;
    C12E2_struct[i+250].index[5] = A3_3_2;
    C12E2_struct[i+250].index[6] = A4_2;
    
    C12E2_struct[i+500].index[0] = A7_3;
    C12E2_struct[i+500].index[1] = A6_1_3;
    C12E2_struct[i+500].index[2] = A6_2_3;
    C12E2_struct[i+500].index[3] = A3_1_3;
    C12E2_struct[i+500].index[4] = A3_2_3;
    C12E2_struct[i+500].index[5] = A3_3_3;
    C12E2_struct[i+500].index[6] = A4_3;

    C12E2_struct[i+750].index[0] = A7_4;
    C12E2_struct[i+750].index[1] = A6_1_4;
    C12E2_struct[i+750].index[2] = A6_2_4;
    C12E2_struct[i+750].index[3] = A3_1_4;
    C12E2_struct[i+750].index[4] = A3_2_4;
    C12E2_struct[i+750].index[5] = A3_3_4;
    C12E2_struct[i+750].index[6] = A4_4;
    
    // Mimic assignment
    
    C12E2M_struct[i].index[0] = A13;
    C12E2M_struct[i].index[1] = A12_1;
    C12E2M_struct[i].index[2] = A12_2;
    C12E2M_struct[i].index[3] = A9_1;
    C12E2M_struct[i].index[4] = A9_2;
    C12E2M_struct[i].index[5] = A9_3;
    C12E2M_struct[i].index[6] = A10;

    C12E2M_struct[i+250].index[0] = A13_2;
    C12E2M_struct[i+250].index[1] = A12_1_2;
    C12E2M_struct[i+250].index[2] = A12_2_2;
    C12E2M_struct[i+250].index[3] = A9_1_2;
    C12E2M_struct[i+250].index[4] = A9_2_2;
    C12E2M_struct[i+250].index[5] = A9_3_2;
    C12E2M_struct[i+250].index[6] = A10_2;
    
    C12E2M_struct[i+500].index[0] = A13_3;
    C12E2M_struct[i+500].index[1] = A12_1_3;
    C12E2M_struct[i+500].index[2] = A12_2_3;
    C12E2M_struct[i+500].index[3] = A9_1_3;
    C12E2M_struct[i+500].index[4] = A9_2_3;
    C12E2M_struct[i+500].index[5] = A9_3_3;
    C12E2M_struct[i+500].index[6] = A10_3;

    C12E2M_struct[i+750].index[0] = A13_4;
    C12E2M_struct[i+750].index[1] = A12_1_4;
    C12E2M_struct[i+750].index[2] = A12_2_4;
    C12E2M_struct[i+750].index[3] = A9_1_4;
    C12E2M_struct[i+750].index[4] = A9_2_4;
    C12E2M_struct[i+750].index[5] = A9_3_4;
    C12E2M_struct[i+750].index[6] = A10_4;

    
    C12E2_skeleton arrM1 = {A13, A12_1, A12_2, A9_1, A9_2, A9_3, A10};  // First array batch
    C12E2_skeleton arrM2 = {A13_2, A12_1_2, A12_2_2, A9_1_2, A9_2_2, A9_3_2, A10_2}; // Second array batch
    C12E2_skeleton arrM3 = {A13_3, A12_1_3, A12_2_3, A9_1_3, A9_2_3, A9_3_3, A10_3}; // Third array batch
    C12E2_skeleton arrM4 = {A13_4, A12_1_4, A12_2_4, A9_1_4, A9_2_4, A9_3_4, A10_4}; // Fourth array batch

    // C12E2M_M.insert(std::pair<int,C12E2_skeleton>(i, arrM1));
    //C12E2M_M.insert(std::pair<int,C12E2_skeleton>(i+1, arrM2));
    // C12E2M_M.insert(std::pair<int,C12E2_skeleton>(i+2, arrM3));
    //C12E2M_M.insert(std::pair<int,C12E2_skeleton>(i+3, arrM4));

  }

  // loop over the values
  
  

  for (int i = 0; i <= 999; i++){

    // Loop to ensure that each struct has been constructed properly 
    // This has been checked
    // std::cout << C12E2M_struct[i].index[1] << " " << C12E2_struct[i].index[1] << " " << i  << std::endl;
  }
  
  for(SSno = 0; SSno < numberofSS; SSno++) {
    
     //  printf("This is the data for trajectory no %d \n", SSno); 
    l = 0;
    n = 0;
    
    for(k=0;k<nlines;k++) { 
   
      fgets(line,sizeof(line),ipf);
      
      //while (fgets(line,sizeof(line),ipf) != NULL) {
      if (l < 5) {	
	/* We are doing nothing */
	//	printf("%s l= %d", line,l);
	l++;
      }
      
      else if ((l > 4 && l < 8)) {
      
	/*We are scanning the bit with just the box parameters*/
	sscanf(line, "%lf %lf", &box1, &box2);
	//	printf("%lf %lf \n",box1,box2 );
	boxlength[l-5] = box2-box1; 
	l++;
      }
      else if (l == 8) {
      	// printf(" **** l= %d \n",l);
	/* We are doing nothing */
	l++;
      }
    
      else {

	//	printf(" ***** l = %d \n",l );
	/* convert the text to numbers */
	sscanf(line,"%d %d %lf %lf %lf",&index,&atomtype,&x,&y,&z);
	a[index-1] = index;
	b[index-1] = atomtype;
	xco[index-1] = x*boxlength[0]; 
	yco[index-1] = y*boxlength[1];
	zco[index-1] = z*boxlength[2];  
	n++;
	l++;
       	//printf("%d %d %lf %lf %lf\n",index,atomtype,x,y,z);
      }
    }
    
    // We can hold a value in the i array, then check it with the j loop
    // we want to put the r threshold (the zone where if in it, we want to consider it as phase separated)
    
    // We need to first catagorize polymer 7 6 3 4 as one type of polymer, and set 12 11 10 9 as the second. 
    // make sure that we do not read atoms of the same molecule 
    

    /* There are 1000 polymers, of C12E2 and mimic type, but because the  */
    
    // --------------------- Array of the center of masses --------------------- //
    
    // We are implementing a pseudo DBSCAN algorithm
    
    // To measure the distance between the points, we measure the distances between //
    // the centers of masses (COM). The DBSCAN algorithm originally requires to     //
    // calculate the `anchor' COMs, which are decided through a nearest neighbour   //
    // calculation. In the schematic below, 
    
    // O  O  O ---> X  -|
    //  \/ \/           | -- the O's show the 'anchor', and the X's show the increasing cluster
    //   O  O ---- X   -|
    
    // Calculating the COM of the C12E2
    // The COM components are divided into the x, y, and z coordinates

      //  std::cout << *it << std::endl;
      //  std::cout << C12E2I.size() << "\n" << std::endl;

    for (int i = 0; i <= sizeof(C12E2_struct)/sizeof(C12E2_struct[1])-1; i++) {
      
      // Calculating the centers of mass for the C12E2      
       headGroupC12_E2xCOM[i] = CenterOfMass(xco[C12E2_struct[i].index[0]],xco[C12E2_struct[i].index[1]],xco[C12E2_struct[i].index[2]],xco[C12E2_struct[i].index[3]],xco[C12E2_struct[i].index[4]],xco[C12E2_struct[i].index[5]],xco[C12E2_struct[i].index[6]]); 
       headGroupC12_E2yCOM[i] = CenterOfMass(yco[C12E2_struct[i].index[0]],yco[C12E2_struct[i].index[1]],yco[C12E2_struct[i].index[2]],yco[C12E2_struct[i].index[3]],yco[C12E2_struct[i].index[4]],yco[C12E2_struct[i].index[5]],yco[C12E2_struct[i].index[6]]); 
       headGroupC12_E2zCOM[i] = CenterOfMass(zco[C12E2_struct[i].index[0]],zco[C12E2_struct[i].index[1]],zco[C12E2_struct[i].index[2]],zco[C12E2_struct[i].index[3]],zco[C12E2_struct[i].index[4]],zco[C12E2_struct[i].index[5]],zco[C12E2_struct[i].index[6]]); 
     
      // Calculating the COM of the C12E2 mimics

       headGroupMIMICxCOM[i] = CenterOfMass(xco[C12E2M_struct[i].index[0]],xco[C12E2M_struct[i].index[1]],xco[C12E2M_struct[i].index[2]],xco[C12E2M_struct[i].index[3]],xco[C12E2M_struct[i].index[4]],xco[C12E2M_struct[i].index[5]],xco[C12E2M_struct[i].index[6]]); 
       headGroupMIMICyCOM[i] = CenterOfMass(yco[C12E2M_struct[i].index[0]],yco[C12E2M_struct[i].index[1]],yco[C12E2M_struct[i].index[2]],yco[C12E2M_struct[i].index[3]],yco[C12E2M_struct[i].index[4]],yco[C12E2M_struct[i].index[5]],yco[C12E2M_struct[i].index[6]]);
       headGroupMIMICzCOM[i] = CenterOfMass(zco[C12E2M_struct[i].index[0]],zco[C12E2M_struct[i].index[1]],zco[C12E2M_struct[i].index[2]],zco[C12E2M_struct[i].index[3]],zco[C12E2M_struct[i].index[4]],zco[C12E2M_struct[i].index[5]],zco[C12E2M_struct[i].index[6]]); 
    
       //std::cout << headGroupC12_E2xCOM[i] << " " << headGroupC12_E2yCOM[i] << " " << headGroupC12_E2zCOM[i]  << " " << i << " " << std::endl;
       // std::cout << headGroupMIMICxCOM[i] << " " << headGroupMIMICyCOM[i] << " " << headGroupMIMICzCOM[i]  << i << " " << std::endl;
       
    }

    /*      
      Calculating the composition profile 
    */

    double tophead = 0;
    double downhead = 0;

    double mimictophead = 0;
    double mimicdownhead = 0;
  
    
    for (int i = 0; i <= sizeof(C12E2_struct)/sizeof(C12E2_struct[1])-1; i++) {

      //double truedist(double COM1x, double COM1y, double COM1z, double COM2x, double COM2y, double COM2z) {

      // See if the distance between the nanopparticle is smaller than 15 
      //std::cout << zco[71312] << std::endl;

       if(truedist(xco[71312], yco[71312], zco[71312], headGroupC12_E2xCOM[i], headGroupC12_E2yCOM[i], headGroupC12_E2zCOM[i]) <= 35.000) { 

	 // Check if the c12e2 molecule is on the top layer or the bottom layer
	 
	 if (headGroupC12_E2zCOM[i] > zco[71312]) {
	   tophead++;
	 }
	 
	 else if (headGroupC12_E2zCOM[i] < zco[71312]) {
	   downhead++;
	 }
	
       }
     
      
       if(truedist(xco[71312], xco[71312], xco[71312], headGroupMIMICxCOM[i], headGroupMIMICyCOM[i], headGroupMIMICzCOM[i] ) <= 35.000) { 
	 
	 if (headGroupMIMICzCOM[i] > zco[71312]) {
	   mimictophead++;
	 }
	 
	 else if (headGroupMIMICzCOM[i] < zco[71312]) {
	   mimicdownhead++;
	 }
       }
             
    }
    
    std::cout << SSno << " " << tophead << " " << downhead << " " <<  mimictophead << " " <<  mimicdownhead <<  std::endl;
	   
    //C12E2M_struct[i+250].index[0] = A13_2;


  // point to the array of the COM
    PolymerCounter = 0;
    MimicCounter = 0;
    
    for (int i = 0; i <= sizeof(C12E2_struct)/sizeof(C12E2_struct[1])-1; i++) {
      for (int j = 0; j <= sizeof(C12E2_struct)/sizeof(C12E2_struct[1])-1; j++) {
	
	if(headGroupC12_E2xCOM[i] !=headGroupC12_E2xCOM[j] && headGroupC12_E2yCOM[i] !=headGroupC12_E2yCOM[j] && headGroupC12_E2zCOM[i] !=headGroupC12_E2zCOM[j]) {
	  
	  if(headGroupC12_E2xCOM[i] - headGroupC12_E2xCOM[j] !=0 && sqrt(pow(headGroupC12_E2xCOM[i] - headGroupC12_E2xCOM[j],2)) <= 7.000 ) {    
	  if(headGroupC12_E2yCOM[i] - headGroupC12_E2yCOM[j] !=0 && sqrt(pow(headGroupC12_E2yCOM[i] - headGroupC12_E2yCOM[j],2)) <= 7.000) {

	    if(headGroupC12_E2zCOM[i] - headGroupC12_E2zCOM[j] !=0 && sqrt(pow(headGroupC12_E2zCOM[i] - headGroupC12_E2zCOM[j],2)) <= 7.000) {
	      
	    // printf("%i %i %lf %lf %lf %i\n",b[i],b[j],headGroupC12_E2xCOM[i] - headGroupC12_E2xCOM[j],headGroupC12_E2yCOM[i] - headGroupC12_E2yCOM[j],headGroupC12_E2zCOM[i] - headGroupC12_E2zCOM[j],NewCounter);
	    
	    
		//NewCounter takes into account the total number of signifcant interactions between like particles 
	    
	    PolymerCounter++;
	    
	    }
	    
	  }
	  
	  }
	  
	}
      }
     
    }

    for (int i = 0; i <= sizeof(C12E2_struct)/sizeof(C12E2_struct[1])-1; i++) {
      for (int j = 0; j <= sizeof(C12E2_struct)/sizeof(C12E2_struct[1])-1; j++) {
	
	if(headGroupMIMICxCOM[i] !=headGroupMIMICxCOM[j] && headGroupMIMICyCOM[i] !=headGroupC12_E2yCOM[j] && headGroupMIMICzCOM[i] !=headGroupMIMICzCOM[j]) {
	  
	  if(headGroupMIMICxCOM[i] - headGroupMIMICxCOM[j] !=0 && sqrt(pow(headGroupMIMICxCOM[i] - headGroupMIMICxCOM[j],2)) <= 7.000 ) {    

	    if(headGroupMIMICyCOM[i] - headGroupMIMICyCOM[j] !=0 && sqrt(pow(headGroupMIMICyCOM[i] - headGroupMIMICyCOM[j],2)) <= 7.000) {

	    if(headGroupMIMICzCOM[i] - headGroupMIMICzCOM[j] !=0 && sqrt(pow(headGroupMIMICzCOM[i] - headGroupMIMICzCOM[j],2)) <= 7.000) {
	      
	    // printf("%i %i %lf %lf %lf %i\n",b[i],b[j],headGroupC12_E2xCOM[i] - headGroupC12_E2xCOM[j],headGroupC12_E2yCOM[i] - headGroupC12_E2yCOM[j],headGroupC12_E2zCOM[i] - headGroupC12_E2zCOM[j],NewCounter);
	    
	    
		//NewCounter takes into account the total number of signifcant interactions between like particles 
	    
	    MimicCounter++;
	    
	    }
	    
	  }
	  
	  }
	  
	}
      }
     
    }
	

	
    //	std::cout << SSno << " " << PolymerCounter << " " << MimicCounter <<  std::endl;

  }
    
  return 0;

    
}
    
/* Function to calculate the centre of mass of each molecule */

double CenterOfMass(double H7, double H6_1, double H6_2, double T3_1, double T3_2, double T3_3, double T4)
  
{  
  double H7coord = H7;
  double H6coord_1 = H6_1;
  double H6coord_2 = H6_2;
  double T3coord_1 = T3_1; 
  double T3coord_2 = T3_2; 
  double T3coord_3 = T3_3; 
  double T4coord = T4; 
  double COM;
  
  COM = (( H7coord ) + ( H6coord_1 ) + ( H6coord_2 ) + ( T3coord_1 )  + ( T3coord_2 )  +( T3coord_3 )  + ( T4coord ))/7; 

  /* With this COM definition we now know the COM in each cartesian coordinate */ 

  return COM; 
}  

double truedist(double COM1x, double COM1y, double COM1z, double COM2x, double COM2y, double COM2z) {

  double dist = pow((pow(COM1x-COM2x,2.0) + pow(COM1y-COM2y,2.0) + pow(COM1z-COM2z,2.0)),0.5);

  return dist;
}

pointCluster* regionQuery (int D, double *xarray1, double *yarray1, double *zarray1, double *xarray2, double *yarray2, double *zarray2) {

  pointCluster A[1000];

  for (int i = 0; i <= sizeof(xarray1)/sizeof(int)-1; i++) {

    for (int j = 0; i <= sizeof(xarray1)/sizeof(int)-1; j++) {
      
      if (i == j) continue;
      
      if (truedist(xarray1[i], xarray1[i], xarray1[i], xarray2[j], xarray2[j], xarray2[j]) <= D) {
	A[i].clusterIndices.push_back(j);
      }      

      return A;
    }
  }
}

//void DBSCAN() {
//}
/*Function to take into account the periodic boundary conditions*/ 


/* void minimum_image(double *xcoord7, double *xcoord6_1, double *xcoord6_2, double *xcoord3_1, double *xcoordin3_2, double xcoordin3_3, double xcoordin4,double boxlength)  */

/* { */
  
/*   double box = xboxlength;  */

/*   //with 7, we move the headgroup  */

/*   if (*xcoord7 - *xcoordx6_1 > boxlength/2 ) { */

/*     *xccord7 = *xccord7 - xboxlength;  */
 
/*   } */

/*   else if (*xcoord7 - *xcoordx6_1 > boxlength/2) { */




/*   } */
  


/* } */


