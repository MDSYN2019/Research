#include <cmath>

double faculty_(int *n)
{
  int i__1;
  double ret_val;
  
  /* Local variables */
  static int j;

  ret_val = 0.0;
  i__1 = *n;
  for (j = 2; j <= i__1; ++j) {
    ret_val += log(j);
  }
  return ret_val;
} 

int main () {
  return 0;
}
