/* serial_dot.c -- compute a dot product on a single processor.
 *
 * Input: 
 *     n: order of vectors
 *     x, y:  the vectors
 *
 * Output:
 *     the dot product of x and y.
 *
 * Note:  Arrays containing vectors are statically allocated.
 *
 * See Chap 5, p. 75 in PPMPI.
 */
#include <stdio.h>

#define MAX_ORDER 100

main() {
    float  x[MAX_ORDER];
    float  y[MAX_ORDER];
    int    n;
    float  dot;

    void Read_vector(char* prompt, float v[], int n);
    float Serial_dot(float x[], float y[], int n);
    
    printf("Enter the order of the vectors\n");
    scanf("%d", &n);
    Read_vector("the first vector", x, n);
    Read_vector("the second vector", y, n);
    dot = Serial_dot(x, y, n);
    printf("The dot product is %f\n", dot);
}  /* main */
   

/*****************************************************************/
void Read_vector(
         char*  prompt  /* in  */,
         float  v[]     /* out */,
         int    n       /* in  */) {
    int i;

    printf("Enter %s\n", prompt);
    for (i = 0; i < n; i++)
        scanf("%f", &v[i]);
}  /* Read_vector */


/*****************************************************************/
float Serial_dot(
          float  x[]  /* in */, 
          float  y[]  /* in */, 
          int    n    /* in */) {

    int    i; 
    float  sum = 0.0;

    for (i = 0; i < n; i++)
        sum = sum + x[i]*y[i];
    return sum;
} /* Serial_dot */
