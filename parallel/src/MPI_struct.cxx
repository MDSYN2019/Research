struct Partstruct 
{ 
  int    class;  /* particle class */ 
  double d[6];   /* particle coordinates */ 
  char   b[7];   /* some additional information */ 
}; 
 
struct Partstruct    particle[1000];  
int                  i, dest, rank; 
MPI_Comm     comm; 
 
 
/* build datatype describing structure */ 

MPI_Datatype Particletype; 
MPI_Datatype type[3] = {MPI_INT, MPI_DOUBLE, MPI_CHAR}; 
int          blocklen[3] = {1, 6, 7}; 
MPI_Aint     disp[3]; 
int          base; 
 
 
/* compute displacements of structure components */ 
 
MPI_Address( particle, disp); 
MPI_Address( particle[0].d, disp+1); 
MPI_Address( particle[0].b, disp+2); 
base = disp[0]; 
for (i=0; i <3; i++) disp[i] -= base; 

MPI_Type_struct( 3, blocklen, disp, type, &Particletype); 
 
/* If compiler does padding in mysterious ways, 
   the following may be safer */ 
 
MPI_Datatype type1[4] = {MPI_INT, MPI_DOUBLE, MPI_CHAR, MPI_UB}; 
int          blocklen1[4] = {1, 6, 7, 1}; 
MPI_Aint     disp1[4]; 
 
/* compute displacements of structure components */ 
 
MPI_Address( particle, disp1); 
MPI_Address( particle[0].d, disp1+1); 
MPI_Address( particle[0].b, disp1+2); 
MPI_Address( particle+1, disp1+3); 
base = disp1[0]; 
for (i=0; i <4; i++) disp1[i] -= base; 
 
/* build datatype describing structure */ 

MPI_Type_struct( 4, blocklen1, disp1, type1, &Particletype); 


/* 4.1: 
   send the entire array */ 

MPI_Type_commit( &Particletype); 
MPI_Send( particle, 1000, Particletype, dest, tag, comm); 
 
 
/* 4.2: 
   send only the entries of class zero particles, 
   preceded by the number of such entries */ 

MPI_Datatype Zparticles;   /* datatype describing all particles 
                              with class zero (needs to be recomputed 
                              if classes change) */ 
MPI_Datatype Ztype; 
 
MPI_Aint     zdisp[1000]; 
int zblock[1000], j, k; 
int zzblock[2] = {1,1}; 
MPI_Aint     zzdisp[2]; 
MPI_Datatype zztype[2]; 
 
/* compute displacements of class zero particles */ 
j = 0; 
for(i=0; i < 1000; i++) 
  if (particle[i].class==0) 
     { 
     zdisp[j] = i; 
     zblock[j] = 1; 
     j++; 
     } 
 
/* create datatype for class zero particles  */ 
MPI_Type_indexed( j, zblock, zdisp, Particletype, &Zparticles); 
 
/* prepend particle count */ 
MPI_Address(&j, zzdisp); 
MPI_Address(particle, zzdisp+1); 
zztype[0] = MPI_INT; 
zztype[1] = Zparticles; 
MPI_Type_struct(2, zzblock, zzdisp, zztype, &Ztype); 
 
MPI_Type_commit( &Ztype); 
MPI_Send( MPI_BOTTOM, 1, Ztype, dest, tag, comm); 
 
 
       /* A probably more efficient way of defining Zparticles */ 
 
/* consecutive particles with index zero are handled as one block */ 
j=0; 
for (i=0; i < 1000; i++) 
  if (particle[i].index==0) 
    { 
    for (k=i+1; (k < 1000)&&(particle[k].index == 0) ; k++); 
    zdisp[j] = i; 
    zblock[j] = k-i; 
    j++; 
    i = k; 
    } 
MPI_Type_indexed( j, zblock, zdisp, Particletype, &Zparticles); 
 
 
                /* 4.3: 
          send the first two coordinates of all entries */ 
 
MPI_Datatype Allpairs;      /* datatype for all pairs of coordinates */ 
 
MPI_Aint sizeofentry; 
 
MPI_Type_extent( Particletype, &sizeofentry); 
 
     /* sizeofentry can also be computed by subtracting the address 
        of particle[0] from the address of particle[1] */ 
 
MPI_Type_hvector( 1000, 2, sizeofentry, MPI_DOUBLE, &Allpairs); 
MPI_Type_commit( &Allpairs); 
MPI_Send( particle[0].d, 1, Allpairs, dest, tag, comm); 
 
      /* an alternative solution to 4.3 */ 
 
MPI_Datatype Onepair;   /* datatype for one pair of coordinates, with 
                          the extent of one particle entry */ 
MPI_Aint disp2[3]; 
MPI_Datatype type2[3] = {MPI_LB, MPI_DOUBLE, MPI_UB}; 
int blocklen2[3] = {1, 2, 1}; 
 
MPI_Address( particle, disp2); 
MPI_Address( particle[0].d, disp2+1); 
MPI_Address( particle+1, disp2+2); 
base = disp2[0]; 
for (i=0; i<2; i++) disp2[i] -= base; 
 
MPI_Type_struct( 3, blocklen2, disp2, type2, &Onepair); 
MPI_Type_commit( &Onepair); 
MPI_Send( particle[0].d, 1000, Onepair, dest, tag, comm); 
