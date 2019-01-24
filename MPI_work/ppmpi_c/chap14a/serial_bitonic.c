/* serial_bitonic.c -- serial bitonic sort of randomly generated list
 *     of integers
 *
 * Input:
 *     n: the length of the list -- must be a power of 2.
 *
 * Output:
 *     The unsorted list and the sorted list.
 *
 * Notes:
 *     1.  The list is statically allocated -- size specified in MAX.
 *     2.  Keys are in the range 0 -- KEY_MAX-1.
 *
 * See Chap 14, pp. 316 & ff. in PPMPI.
 */

#include <stdio.h>
#include <stdlib.h>
   /* to get rand() */

#define MAX 16384

/* Successive subsequences will switch between
 * increasing and decreasing bitonic splits.
 */
#define INCR 0
#define DECR 1
#define Reverse(ordering) ((ordering) == INCR ? DECR : INCR)

typedef int KEY_T;
#define Swap(a,b) {KEY_T temp; temp = a; a = b; b = temp;}
#define KEY_MAX 32768 

void Generate_list(int n, KEY_T A[]);
void Bitonic_sort_incr(int length, KEY_T B[]);
void Bitonic_sort_decr(int length, KEY_T B[]);
void Print_list(char* title, int n, KEY_T A[]);

/*********************************************************************/
main() {
    int    list_length;
    int    n;
    int    start_index;
    int    ordering;
    KEY_T  A[MAX];

    printf("Enter the list size (a power of 2)\n");
    scanf("%d", &n);

    Generate_list(n, A);

    Print_list("The unsorted list is", n, A);

    for (list_length = 2; list_length <= n; 
         list_length = list_length*2)
        for (start_index = 0, ordering = INCR; 
             start_index < n;
             start_index = start_index + list_length, 
             ordering = Reverse(ordering))
            if (ordering == INCR)
                Bitonic_sort_incr(list_length, 
                    A + start_index);
            else
                Bitonic_sort_decr(list_length, 
                    A + start_index);

    Print_list("The sorted list is", n, A);

}  /* main */


/*********************************************************************/
void Generate_list(int n, KEY_T A[]) {
    int i;

    for (i = 0; i < n; i++)
        A[i] = rand() % KEY_MAX; 
}  /* Generate_list */ 


/*********************************************************************/
void Bitonic_sort_incr(
        int     length  /* in     */, 
        KEY_T   B[]     /* in/out */) {
    int i;
    int half_way;

    /* This is the bitonic split */
    half_way = length/2;
    for (i = 0; i < half_way; i++)
        if (B[i] > B[half_way + i])
            Swap(B[i],B[half_way+i]);

    if (length > 2) {
        Bitonic_sort_incr(length/2, B);
        Bitonic_sort_incr(length/2, B + half_way);
    }
} /* Bitonic_sort_incr */


/*********************************************************************/
void Bitonic_sort_decr(
        int    length   /* in     */, 
        KEY_T  B[]      /* in/out */) {
    int i;
    int half_way;

    /* This is the bitonic split */
    half_way = length/2;
    for (i = 0; i < half_way; i++)
        if (B[i] < B[half_way + i])
            Swap(B[i],B[half_way+i]);

    if (length > 2) {
        Bitonic_sort_decr(length/2, B);
        Bitonic_sort_decr(length/2, B + half_way);
    }
} /* Bitonic_sort_decr */


/*********************************************************************/
void Print_list(char* title, int n, KEY_T A[]) {
    int i;

    printf("%s\n", title);
    for (i = 0; i < n; i++)
        printf("%d ", A[i]);
    printf("\n");
}  /* Print_list */
