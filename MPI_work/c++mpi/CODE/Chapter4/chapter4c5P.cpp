
#include<iostream.h>
#include<math.h>
#include<mpi.h>
#include "SCchapter4.h"

double func(double x);

int main(int argc, char * argv[]){
  int mynode, totalnodes;
  const double global_a = -50.0;
  const double global_b =  50.0;
  const int levels = 10;
  double local_a,local_b,local_sum,answer;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  local_a = global_a + mynode    *(global_b-global_a)/totalnodes;
  local_b = global_a + (mynode+1)*(global_b-global_a)/totalnodes;

  local_sum = MidpointRule(levels, local_a, local_b, func);

  MPI_Reduce(&local_sum,&answer,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if(mynode == 0){
    cout << "The value of the integral is: " << answer << endl;
  }

  MPI_Finalize();
  
}

double func(double x){
  return(sin(x)/x);
}






