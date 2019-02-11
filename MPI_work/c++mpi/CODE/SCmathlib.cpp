/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCmathlib.cpp (definition file)
* Header file:  SCmathlib.h   (header file)      
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
* Description:
*
* Available Classes Declarations/Definitions:
*    SCPoint
*    SCVector
*    SCMatrix
*
*************************************************************************/


#include "SCmathlib.h"

SCPoint::SCPoint(int dim){
  dimension = dim;
  data = new double[dimension];

  for(int i=0;i<dimension;i++)
    data[i] = 0.0;
}


SCPoint::SCPoint(const SCPoint &v){
  dimension = v.Dimension();
  data = new double[dimension];

  for(int i=0;i<dimension;i++)
    data[i] = v.data[i];
}


SCPoint::~SCPoint(){
  dimension = 0;
  delete[] data;
  data = NULL;
}


int SCPoint::Dimension() const{
  return(dimension);
}


double SCPoint::operator()(const int i) const{
  if(i>=0 && i<dimension)
    return data[i];

  cerr << "SCPoint::Invalid index " << i << " for SCPoint of dimension " << dimension << endl;
  return(0);
}



double& SCPoint::operator()(const int i){
  if(i>=0 && i<dimension)
    return data[i];

  cerr << "SCPoint::Invalid index " << i << " for SCPoint of dimension " << dimension << endl;
  return(data[0]);
}


SCPoint& SCPoint::operator=(const SCPoint &v) {
  dimension = v.Dimension();
  for(int i=0;i<dimension;i++)
    data[i] = v.data[i];
  return *this;
};

void SCPoint::Print() const{
  cout << endl;
  cout << "[ ";
  if(dimension>0)
    cout << data[0];
  for(int i=1;i<dimension;i++)
    cout << "; " << data[i];
  cout << " ]" << endl;
}

SCVector::SCVector(){
  dimension = 0;
  data = NULL;
}


SCVector::SCVector(int dim){
  dimension = dim;
  data = new double[dimension];

  for(int i=0;i<dimension;i++)
    data[i] = 0.0;
}


SCVector::SCVector(const SCVector &v){
  dimension = v.Dimension();
  data = new double[dimension];

  for(int i=0;i<dimension;i++)
    data[i] = v.data[i];
}


SCVector::SCVector(int col, const SCMatrix &A){
  dimension = A.Rows();

  data = new double[dimension];
  
  for(int i=0;i<A.Rows();i++)
    data[i] = A(i,col);

}


SCVector::~SCVector(){
  dimension = 0;
  delete[] data;
  data = NULL;
}


void SCVector::Initialize(int dim){
  if(dimension!=0)
    delete[] data;

  dimension = dim;
  data = new double[dimension];
  
  for(int i=0;i<dimension;i++)
    data[i] = 0.0;
}


int SCVector::Dimension() const{
  return(dimension);
}


double SCVector::operator()(const int i) const{
  if(i>=0 && i<dimension)
    return data[i];

  cerr << "SCVector::Invalid index " << i << " for SCVector of dimension " << dimension << endl;
  return(0);
}



double& SCVector::operator()(const int i){
  if(i>=0 && i<dimension)
    return data[i];

  cerr << "SCVector::Invalid index " << i << " for SCVector of dimension " << dimension << endl;
  return(data[0]);
}


SCVector& SCVector::operator=(const SCVector &v) {
  dimension = v.Dimension();
  for(int i=0;i<dimension;i++)
    data[i] = v.data[i];
  return *this;
};

void SCVector::Print() const{
  cout << endl;
  cout << "[ ";
  if(dimension>0)
    cout << data[0];
  for(int i=1;i<dimension;i++)
    cout << "; " << data[i];
  cout << " ]" << endl;
}


double SCVector::Norm_l1(){
  double sum = 0.0;
  for(int i=0;i<dimension;i++)
    sum += fabs(data[i]);
  return(sum);
}


double SCVector::Norm_l2(){
  double sum = 0.0;
  for(int i=0;i<dimension;i++)
    sum += data[i]*data[i];
  return(sqrt(sum));
}

void SCVector::Normalize(){
  double tmp = 1.0/Norm_l2();
  for(int i=0;i<dimension;i++)
    data[i] = data[i]*tmp;
}


double SCVector::Norm_linf(){
  double maxval = 0.0,tmp;
  
  for(int i=0;i<dimension;i++){
    tmp = fabs(data[i]);
    maxval = (maxval > tmp)?maxval:tmp;
  }
  return(maxval);
}

double SCVector::MaxMod(){
  double maxm = -1.0e+10;

  for(int i=0; i<dimension; i++)
    maxm = (maxm > fabs(data[i]))?maxm:fabs(data[i]);
  
  return maxm;
}

double SCVector::ElementofMaxMod(){
  return(data[MaxModindex()]);
}


int SCVector::MaxModindex(){
  double maxm = -1.0e+10;
  int maxmindex = 0;

  for(int i=0; i<dimension; i++){
    if(maxm<fabs(data[i])){
      maxm = fabs(data[i]);
      maxmindex = i;
    }
  }
  
  return maxmindex;
}

void SCVector::Initialize(double a){
  for(int i=0; i<dimension; i++)
    data[i] = a;
}

void SCVector::Initialize(double *v){
  for(int i=0; i<dimension; i++)
    data[i] = v[i];
}

SCMatrix::SCMatrix(int dim){
  rows = dim;
  columns = dim;
  data = new double*[rows];
  for(int i=0;i<rows;i++){
    data[i] = new double[columns];
    for(int j=0;j<columns;j++)
      data[i][j] = 0.0;
  }
}

  
SCMatrix::SCMatrix(int rows1, int columns1){
  rows = rows1;
  columns = columns1;

  data = new double*[rows];
  for(int i=0;i<rows;i++){
    data[i] = new double[columns];
    for(int j=0;j<columns;j++)
      data[i][j] = 0.0;
  }
}

SCMatrix::SCMatrix(const SCMatrix& m){
  rows = m.rows;
  columns = m.columns;

  data = new double*[rows];

  for(int i=0;i<rows;i++){
    data[i] = new double[columns];
    for(int j=0; j<columns; j++)
      data[i][j] = m.data[i][j];
  }
}

SCMatrix::SCMatrix(int num_SCVectors, const SCVector * q){
  rows = q[0].Dimension();
  columns = num_SCVectors;

  data = new double*[rows];

  for(int i=0;i<rows;i++){
    data[i] = new double[columns];
    for(int j=0; j<columns; j++)
      data[i][j] = q[j](i);
  }
}

SCMatrix::SCMatrix(int rows1, int columns1, double **rowptrs){
  rows = rows1;
  columns = columns1;

  data = new double*[rows];

  for(int i=0;i<rows;i++)
    data[i] = rowptrs[i];
}


SCMatrix::~SCMatrix(){
  for(int i=0;i<rows;i++)
    delete[] data[i];

  rows = 0;
  columns = 0;
  delete[] data;
}

int SCMatrix::Rows() const{
  return(rows);
}  

int SCMatrix::Columns() const{
  return(columns);
}  


double **SCMatrix::GetPointer(){
  return(data);
}

void SCMatrix::GetColumn(int col, SCVector &x){
  x.Initialize(0.0);
  for(int i=0;i<rows;i++)
    x(i) = data[i][col];
}

void SCMatrix::GetColumn(int col, SCVector &x, int rowoffset){
  x.Initialize(0.0);
  for(int i=0;i<rows-rowoffset;i++)
    x(i) = data[i+rowoffset][col];
}

void SCMatrix::PutColumn(int col, const SCVector &x){
  for(int i=0;i<rows;i++)
    data[i][col] = x(i);
}


double SCMatrix::Norm_linf(){
  double maxval = 0.0,sum;
  
  for(int i=0;i<rows;i++){
    sum = 0.0;
    for(int j=0;j<columns;j++)
      sum += fabs(data[i][j]);
    maxval = (maxval > sum)?maxval:sum;
  }
  return(maxval);
}


double SCMatrix::Norm_l1(){
  double maxval = 0.0,sum;

  for(int j=0;j<columns;j++){
    sum = 0.0;
    for(int i=0;i<rows;i++)
      sum += fabs(data[i][j]);
    maxval = (maxval > sum)?maxval:sum;
  }
  return(maxval);
}



SCMatrix& SCMatrix::operator=(const SCMatrix &m){
  if( (rows == m.rows) && (columns == m.columns)){
    for(int i=0; i<rows; i++)
      for(int j=0;j<columns;j++){
	data[i][j] = m.data[i][j];
      }
  }
  else
    cerr << "SCMatrix Error: Cannot equate matrices of different sizes\n";
  return *this;
}

  
double SCMatrix::operator()(const int i, const int j) const {
  if( (i>=0) && (j>=0) && (i<rows) && (j<columns))
    return(data[i][j]);  
  else
    cerr << "SCMatrix Error: Invalid SCMatrix indices (" << i << "," << j << 
      "), for SCMatrix of size " << rows << " X " << columns << endl;
  return((double)0);
}
  

double& SCMatrix::operator()(const int i, const int j) {
  if( (i>=0) && (j>=0) && (i<rows) && (j<columns))
    return(data[i][j]);  
  else
    cerr << "SCMatrix Error: Invalid SCMatrix indices (" << i << "," << j << 
      "), for SCMatrix of size " << rows << " X " << columns << endl;;
  return(data[0][0]);
}


void SCMatrix::Print() const{
  cout << endl;


  cout << "[ ";
  for(int i=0;i<rows;i++){
    cout << data[i][0];
    for(int j=1;j<columns;j++)
      cout << " " << data[i][j];
    if(i!=(rows-1))
      cout << ";\n";
  }
  cout << " ]" << endl;
}


double SCMatrix::MaxModInRow(int row){
  double maxv = -1.0e+10;
  for(int i=0;i<columns;i++)
    maxv = (fabs(data[row][i])>maxv)?fabs(data[row][i]):maxv;

  return maxv;
}

double SCMatrix::MaxModInRow(int row, int starting_column){
  double maxv = -1.0e+10;
  for(int i=starting_column;i<columns;i++)
    maxv = (fabs(data[row][i])>maxv)?fabs(data[row][i]):maxv;

  return maxv;
}

int SCMatrix::MaxModInRowindex(int row){
  int maxvindex = 0;
  double maxv = -1.0e+10;
  
  for(int i=0;i<columns;i++){
    if(maxv < fabs(data[row][i])){
      maxv = fabs(data[row][i]);
      maxvindex = i;
    }
  }

  return maxvindex;
}

int SCMatrix::MaxModInRowindex(int row, int starting_column){
  int maxvindex = 0;
  double maxv = -1.0e+10;

  for(int i=starting_column;i<columns;i++){
    if(maxv < fabs(data[row][i])){
      maxv = fabs(data[row][i]);
      maxvindex = i;
    }
  }
  
  return maxvindex;
}

double SCMatrix::MaxModInColumn(int column){
  double maxv = -1.0e+10;
  for(int i=0;i<rows;i++)
    maxv = (fabs(data[i][column])>maxv)?fabs(data[i][column]):maxv;

  return maxv;
}

double SCMatrix::MaxModInColumn(int column, int starting_row){
  double maxv = -1.0e+10;
  for(int i=starting_row;i<rows;i++)
    maxv = (fabs(data[i][column])>maxv)?fabs(data[i][column]):maxv;

  return maxv;
}

int SCMatrix::MaxModInColumnindex(int column){
  int maxvindex = 0;
  double maxv = -1.0e+10;
  
  for(int i=0;i<rows;i++){
    if(maxv < fabs(data[i][column])){
      maxv = fabs(data[i][column]);
      maxvindex = i;
    }
  }

  return maxvindex;
}

int SCMatrix::MaxModInColumnindex(int column, int starting_column){
  int maxvindex = 0;
  double maxv = -1.0e+10;

  for(int i=starting_column;i<rows;i++){
    if(maxv < fabs(data[i][column])){
      maxv = fabs(data[i][column]);
      maxvindex = i;
    }
  }
  
  return maxvindex;
}

void SCMatrix::RowSwap(int row1, int row2){
  double * tmp = data[row1];
  data[row1] = data[row2];
  data[row2] = tmp;
}



/****************************************************************/
/*                 Operator Definitions                         */
/****************************************************************/


SCVector operator-(const SCVector& v){
  SCVector x(v.Dimension());
  for(int i=0;i<v.Dimension();i++)
    x(i) = -v(i);
  return x;
}


SCVector operator+(const SCVector& v1, const SCVector& v2){
  int min_dim = min_dimension(v1,v2);
  SCVector x(min_dim);
  for(int i=0;i<min_dim;i++)
    x(i) = v1(i) + v2(i);
  return x;
}


SCVector operator-(const SCVector& v1, const SCVector& v2){
  int min_dim = min_dimension(v1,v2);
  SCVector x(min_dim);
  for(int i=0;i<min_dim;i++)
    x(i) = v1(i) - v2(i);
  return x;
}


SCVector operator/(const SCVector& v, const double s) {
  SCVector x(v.Dimension());
  for(int i=0;i<v.Dimension();i++)
    x(i) = v(i)/s;
  return x;
}



SCVector operator*(const double s, const SCVector &v) {
  SCVector x(v.Dimension());
  for(int i=0;i<v.Dimension();i++)
    x(i) = s*v(i);
  return x;
}


SCVector operator*(const SCVector& v, const double s) {
  SCVector x(v.Dimension());
  for(int i=0;i<v.Dimension();i++)
    x(i) = s*v(i);
  return x;
}

SCVector operator*(const SCMatrix& A, const SCVector& x){
  int rows = A.Rows(), columns = A.Columns();
  int dim = x.Dimension();
  SCVector b(dim);
  
  if(columns != dim){
    cerr << "Invalid dimensions given in matrix-vector multiply" << endl;
    return(b);
  }
  
  for(int i=0;i<rows;i++){
    b(i) = 0.0;
    for(int j=0;j<columns;j++){
      b(i) += A(i,j)*x(j);
    }
  }
  
  return b;
}


/****************************************************************/
/*                 Function Definitions                         */
/****************************************************************/

int min_dimension(const SCVector& v1, const SCVector& v2){
  int min_dim = (v1.Dimension()<v2.Dimension())?v1.Dimension():v2.Dimension();
  return(min_dim);
}


double dot(const SCVector& u, const SCVector& v){
  double sum = 0.0;
  int min_dim = min_dimension(u,v);

  for(int i=0;i<min_dim;i++)
    sum += u(i)*v(i);
  
  return sum; 
}


double dot(int N, const SCVector& u, const SCVector& v){
  double sum = 0.0;

  for(int i=0;i<N;i++)
    sum += u(i)*v(i);
  
  return sum;
}


double dot(int N, double *a, double *b){
  double sum = 0.0;
  
  for(int i=0;i<N;i++)
    sum += a[i]*b[i];

  return sum;
}


/*******************************/
/*   Log base 2 of a number    */
/*******************************/

double log2(double x){
  return(log(x)/log(2.0));
}

void Swap(double &a, double &b){
  double tmp = a;
  a = b;
  b = tmp;
}

double Sign(double x){
  double xs;

  xs = (x>=0.0)?1.0:-1.0;

  return xs;
}

//GammaF function valid for x integer, or x (integer+0.5)
double GammaF(double x){
  double gamma = 1.0;
  
  if (x == -0.5) 
    gamma = -2.0*sqrt(M_PI);
  else if (!x) return gamma;
  else if ((x-(int)x) == 0.5){ 
    int n = (int) x;
    double tmp = x;
    
    gamma = sqrt(M_PI);
    while(n--){
      tmp   -= 1.0;
      gamma *= tmp;
    }
  }
  else if ((x-(int)x) == 0.0){
    int n = (int) x;
    double tmp = x;
    
    while(--n){
      tmp   -= 1.0;
      gamma *= tmp;
    }
  }  
  
  return gamma;
}


int Factorial(int n){
  int value=1;
  for(int i=n;i>0;i--)
    value = value*i;

  return value;
}

double ** CreateMatrix(int m, int n){
  double ** mat;
  mat = new double*[m];
  for(int i=0;i<m;i++){
    mat[i] = new double[n];
    for(int j=0;j<m;j++)
      mat[i][j] = 0.0;
  }
  return mat;
}

int ** ICreateMatrix(int m, int n){
  int ** mat;
  mat = new int*[m];
  for(int i=0;i<m;i++){
    mat[i] = new int[n];
    for(int j=0;j<m;j++)
      mat[i][j] = 0;
  }
  return mat;
}

void DestroyMatrix(double ** mat, int m, int n){
  for(int i=0;i<m;i++)
    delete[] mat[i];
  delete[] mat;
}

void IDestroyMatrix(int ** mat, int m, int n){
  for(int i=0;i<m;i++)
    delete[] mat[i];
  delete[] mat;
}







