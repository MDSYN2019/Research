#ifdef __PP__
#define __PP__

class MPI_PP {
public:
  MPI_PP();
  MPI_PP(int, int);
  virtual ~MPI_PP();
  float serial_dot(float*, float*, int);
  
private:
  int world_rank;
  int world_size;
  int i;
  float sum = 0.0;

  
};

#endif
