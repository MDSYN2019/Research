#ifdef __PP__
#define __PP__

class MPI_PP {
public:
  MPI_PP();
  MPI_PP(int, int);
  virtual ~MPI_PP();

private:
  int world_rank;
  int world_size;
  
};

#endif
