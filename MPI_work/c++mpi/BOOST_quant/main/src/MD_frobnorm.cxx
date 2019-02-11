#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
// ...
ublas::matrix<double> matrix(3, 3);
int k = 0;
int i;
int j;

void matrixInit() {
for (int i = 0; i < matrix.size1(); i++) { 
  for (int j = 0; j < matrix.size2(); j++) {
    matrix(i, j) = ++k;
  }
 }
}
// ...

double frobeniusNorm(const ublas::matrix<double>& matrix)
{
  double result = 0.0;
  for(unsigned int i = 0; i < matrix.size1(); ++i)
    {
      for(unsigned int j = 0; j < matrix.size2(); ++j)
	{
	  double value = matrix(i, j);
	  result += value * value;
	}
    }
  return sqrt(result);
}

